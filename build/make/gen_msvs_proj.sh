#!/bin/bash
##
##  Copyright (c) 2010 The VP8 project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##


self=$0
self_basename=${self##*/}
self_dirname=$(dirname "$0")
EOL=$'\n'

show_help() {
    cat <<EOF
Usage: ${self_basename} --name=projname [options] file1 [file2 ...]

This script generates a Visual Studio project file from a list of source
code files.

Options:
    --help                      Print this message
    --exe                       Generate a project for building an Application
    --lib                       Generate a project for creating a static library
    --static-crt                Use the static C runtime (/MT)
    --target=isa-os-cc          Target specifier (required)
    --out=filename              Write output to a file [stdout]
    --name=project_name         Name of the project (required)
    --proj-guid=GUID            GUID to use for the project
    --module-def=filename       File containing export definitions (for DLLs)
    --ver=version               Version (7,8) of visual studio to generate for
    -Ipath/to/include           Additional include directories
    -DFLAG[=value]              Preprocessor macros to define
    -Lpath/to/lib               Additional library search paths
    -llibname                   Library to link against
EOF
    exit 1
}

die() {
    echo "${self_basename}: $@" >&2
    exit 1
}

die_unknown(){
    echo "Unknown option \"$1\"." >&2
    echo "See ${self_basename} --help for available options." >&2
    exit 1
}

generate_uuid() {
    local hex="0123456789ABCDEF"
    local i
    local uuid=""
    local j
    #93995380-89BD-4b04-88EB-625FBE52EBFB
    for ((i=0; i<32; i++)); do
        (( j = $RANDOM % 16 ))
        uuid="${uuid}${hex:$j:1}"
    done
    echo "${uuid:0:8}-${uuid:8:4}-${uuid:12:4}-${uuid:16:4}-${uuid:20:12}"
}

indent1="    "
indent=""
indent_push() {
    indent="${indent}${indent1}"
}
indent_pop() {
    indent="${indent%${indent1}}"
}

tag_attributes() {
    for opt in "$@"; do
        optval="${opt#*=}"
        [ -n "${optval}" ] ||
            die "Missing attribute value in '$opt' while generating $tag tag"
        echo "${indent}${opt%%=*}=\"${optval}\""
    done
}

open_tag() {
    local tag=$1
    shift
    if [ $# -ne 0 ]; then
        echo "${indent}<${tag}"
        indent_push
        tag_attributes "$@"
        echo "${indent}>"
    else
        echo "${indent}<${tag}>"
        indent_push
    fi
}

close_tag() {
    local tag=$1
    indent_pop
    echo "${indent}</${tag}>"
}

tag() {
    local tag=$1
    shift
    if [ $# -ne 0 ]; then
        echo "${indent}<${tag}"
        indent_push
        tag_attributes "$@"
        indent_pop
        echo "${indent}/>"
    else
        echo "${indent}<${tag}/>"
    fi
}

generate_filter() {
    local var=$1
    local name=$2
    local pats=$3
    local file_list_sz
    local i
    local f
    local saveIFS="$IFS"
    local pack
    echo "generating filter '$name' from ${#file_list[@]} files" >&2
    IFS=*

    open_tag Filter \
        Name=$name \
        Filter=$pats \
        UniqueIdentifier=`generate_uuid`

    file_list_sz=${#file_list[@]}
    for i in ${!file_list[@]}; do
        f=${file_list[i]}
        for pat in ${pats//;/$IFS}; do
            if [ "${f##*.}" == "$pat" ]; then
                unset file_list[i]

                open_tag File RelativePath="./$f"
                if [ "$pat" == "asm" ] && $asm_use_custom_step; then
                    for plat in "${platforms[@]}"; do
                        for cfg in Debug Release; do
                            open_tag  FileConfiguration \
                            Name="${cfg}|${plat}"
                            tag Tool \
                                Name="VCCustomBuildTool" \
                                Description="Assembling \$(InputFileName)" \
                                CommandLine="$(eval echo \$asm_${cfg}_cmdline)"\
                                Outputs="\$(InputName).obj"
                            close_tag FileConfiguration
                        done
                    done
                fi

                if [ "${f##*.}" == "cpp" ]; then
                    for plat in "${platforms[@]}"; do
                        for cfg in Debug Release; do
                        open_tag FileConfiguration \
                            Name="${cfg}|${plat}"
                        tag Tool \
                            Name="VCCLCompilerTool" \
                            CompileAs="2"
                        close_tag FileConfiguration
                        done
                    done
                fi
                close_tag  File

                break
            fi
        done
    done

    close_tag Filter
    IFS="$saveIFS"
}

# Process command line
unset target
for opt in "$@"; do
    optval="${opt#*=}"
    case "$opt" in
    --help|-h) show_help
    ;;
    --target=*) target="${optval}"
    ;;
    --out=*) outfile="$optval"
    ;;
    --name=*) name="${optval}"
    ;;
    --proj-guid=*) guid="${optval}"
    ;;
    --module-def=*)
        link_opts="${link_opts} ModuleDefinitionFile=${optval}"
    ;;
    --exe) proj_kind="exe"
    ;;
    --lib) proj_kind="lib"
    ;;
    --static-crt) use_static_runtime=true
    ;;
    --ver=*) vs_ver="$optval"
             case $optval in
             [78])
             ;;
             *) die Unrecognized Visual Studio Version in $opt
             ;;
             esac
    ;;
    -I*) opt="${opt%/}"
         incs="${incs}${incs:+;}&quot;${opt##-I}&quot;"
         yasmincs="${yasmincs} ${opt}"
    ;;
    -D*) defines="${defines}${defines:+;}${opt##-D}"
    ;;
    -L*) # fudge . to $(OutDir)
         if [ "${opt##-L}" == "." ]; then
             libdirs="${libdirs}${libdirs:+;}&quot;\$(OutDir)&quot;"
         else
             # Also try directories for this platform/configuration
             libdirs="${libdirs}${libdirs:+;}&quot;${opt##-L}&quot;"
             libdirs="${libdirs}${libdirs:+;}&quot;${opt##-L}/\$(PlatformName)/\$(ConfigurationName)&quot;"
             libdirs="${libdirs}${libdirs:+;}&quot;${opt##-L}/\$(PlatformName)&quot;"
         fi
    ;;
    -l*) libs="${libs}${libs:+ }${opt##-l}.lib"
    ;;
    -*) die_unknown $opt
    ;;
    *) file_list[${#file_list[@]}]="$opt"
       case "$opt" in
       *.asm) uses_asm=true;;
       esac
    esac
done
outfile=${outfile:-/dev/stdout}
guid=${guid:-`generate_uuid`}
asm_use_custom_step=false
uses_asm=${uses_asm:-false}
case "${vs_ver:-8}" in
    7) vs_ver_id="7.10"
       asm_use_custom_step=$uses_asm
    ;;
    8) vs_ver_id="8.00"
    ;;
esac

[ -n "$name" ] || die "Project name (--name) must be specified!"
[ -n "$target" ] || die "Target (--target) must be specified!"

if ${use_static_runtime:-false}; then
    release_runtime=0
    debug_runtime=1
    lib_sfx=mt
else
    release_runtime=2
    debug_runtime=3
    lib_sfx=md
fi

# Calculate debug lib names: If a lib ends in ${lib_sfx}.lib, then rename
# it to ${lib_sfx}d.lib. This precludes linking to release libs from a
# debug exe, so this may need to be refactored later.
for lib in ${libs}; do
    if [ "$lib" != "${lib%${lib_sfx}.lib}" ]; then
        lib=${lib%.lib}d.lib
    fi
    debug_libs="${debug_libs}${debug_libs:+ }${lib}"
done


# List Keyword for this target
case "$target" in
    x86*)
        keyword="ManagedCProj"
    ;;
    arm*|iwmmx*)
        keyword="Win32Proj"
    ;;
    *) die "Unsupported target $target!"
esac

# List of all platforms supported for this target
case "$target" in
    x86_64*)
        platforms[0]="x64"
    ;;
    x86*)
        platforms[0]="Win32"
        # these are only used by vs7
        asm_Debug_cmdline="yasm -Xvc -g cv8 -f \$(PlatformName) ${yasmincs} &quot;\$(InputPath)&quot;"
        asm_Release_cmdline="yasm -Xvc -f \$(PlatformName) ${yasmincs} &quot;\$(InputPath)&quot;"
    ;;
    arm*|iwmmx*)
        case "${name}" in
        obj_int_extract) platforms[0]="Win32"
        ;;
        *) platforms[0]="Pocket PC 2003 (ARMV4)"
        ;;
        esac
    ;;
    *) die "Unsupported target $target!"
esac

# List Command-line Arguments for this target
case "$target" in
    arm*|iwmmx*)
        if [ "$name" == "example" ];then
            ARGU="--codec vp6 --flipuv --progress _bnd.vp6"
        fi
        if [ "$name" == "xma" ];then
            ARGU="--codec vp6 -h 240 -w 320 -v"
        fi
    ;;
esac

generate_vcproj() {
    case "$proj_kind" in
    exe) vs_ConfigurationType=1
    ;;
    *)   vs_ConfigurationType=4
    ;;
    esac

    echo "<?xml version=\"1.0\" encoding=\"Windows-1252\"?>"
    open_tag  VisualStudioProject \
                  ProjectType="Visual C++" \
                  Version="${vs_ver_id}" \
                  Name="${name}" \
                  ProjectGUID="{${guid}}" \
                  RootNamespace="${name}" \
                  Keyword="${keyword}"

    open_tag  Platforms
    for plat in "${platforms[@]}"; do
        tag   Platform Name="$plat"
    done
    close_tag Platforms

    open_tag  ToolFiles
    case "$target" in
        x86*) $uses_asm && tag ToolFile RelativePath="$self_dirname/../x86-msvs/yasm.rules"
        ;;
        arm*|iwmmx*)
            if [ "$name" == "vpx_decoder" ];then
            case "$target" in
                armv5*)
                    tag ToolFile RelativePath="$self_dirname/../arm-wince-vs8/armasmv5.rules"
                ;;
                armv6*)
                    tag ToolFile RelativePath="$self_dirname/../arm-wince-vs8/armasmv6.rules"
                ;;
                iwmmxt*)
                    tag ToolFile RelativePath="$self_dirname/../arm-wince-vs8/armasmxscale.rules"
                ;;
            esac
            fi
        ;;
    esac
    close_tag ToolFiles

    open_tag  Configurations
    for plat in "${platforms[@]}"; do
        plat_no_ws=`echo $plat | sed 's/[^A-Za-z0-9_]/_/g'`
        open_tag  Configuration \
                      Name="Debug|$plat" \
                      OutputDirectory="\$(SolutionDir)$plat_no_ws/\$(ConfigurationName)" \
                      IntermediateDirectory="$plat_no_ws/\$(ConfigurationName)/${name}" \
                      ConfigurationType="$vs_ConfigurationType" \
                      CharacterSet="1"

        if [ "$target" == "armv6-wince-vs8" ] || [ "$target" == "armv5te-wince-vs8" ] || [ "$target" == "iwmmxt-wince-vs8" ] || [ "$target" == "iwmmxt2-wince-vs8" ];then
            case "$name" in
                vpx_decoder) tag Tool \
                             Name="VCPreBuildEventTool" \
                             CommandLine="call obj_int_extract.bat \$(ConfigurationName)"
                             tag Tool \
                             Name="VCMIDLTool" \
                             TargetEnvironment="1"
                             tag Tool \
                             Name="VCCLCompilerTool" \
                             ExecutionBucket="7" \
                             Optimization="0" \
                             AdditionalIncludeDirectories="$incs" \
                             PreprocessorDefinitions="_DEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES);WINCE;DEBUG;_LIB;\$(ARCHFAM);\$(_ARCHFAM_);_UNICODE;UNICODE;" \
                             MinimalRebuild="true" \
                             RuntimeLibrary="1" \
                             BufferSecurityCheck="false" \
                             UsePrecompiledHeader="0" \
                             WarningLevel="3" \
                             DebugInformationFormat="1" \
                             CompileAs="1"
                             tag Tool \
                             Name="VCResourceCompilerTool" \
                             PreprocessorDefinitions="_DEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES)" \
                             Culture="1033" \
                             AdditionalIncludeDirectories="\$(IntDir)" \
                ;;
                example|xma) tag Tool \
                             Name="VCCLCompilerTool" \
                             ExecutionBucket="7" \
                             Optimization="0" \
                             AdditionalIncludeDirectories="$incs" \
                             PreprocessorDefinitions="_DEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES);WINCE;DEBUG;_CONSOLE;\$(ARCHFAM);\$(_ARCHFAM_);_UNICODE;UNICODE;" \
                             MinimalRebuild="true" \
                             RuntimeLibrary="1" \
                             BufferSecurityCheck="false" \
                             UsePrecompiledHeader="0" \
                             WarningLevel="3" \
                             DebugInformationFormat="1" \
                             CompileAs="1"
                             tag Tool \
                             Name="VCResourceCompilerTool" \
                             PreprocessorDefinitions="_DEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES)" \
                             Culture="1033" \
                             AdditionalIncludeDirectories="\$(IntDir)" \
                ;;
                obj_int_extract) tag Tool \
                             Name="VCCLCompilerTool" \
                             Optimization="0" \
                             AdditionalIncludeDirectories="$incs" \
                             PreprocessorDefinitions="WIN32;DEBUG;_CONSOLE" \
                             RuntimeLibrary="1" \
                             WarningLevel="3" \
                             DebugInformationFormat="1" \
                ;;
            esac
        fi

        case "$target" in
            x86*) tag Tool \
                Name="VCCLCompilerTool" \
                Optimization="0" \
                AdditionalIncludeDirectories="$incs" \
                PreprocessorDefinitions="WIN32;_DEBUG;_CRT_SECURE_NO_WARNINGS;$defines" \
                RuntimeLibrary="$debug_runtime" \
                UsePrecompiledHeader="0" \
                WarningLevel="3" \
                DebugInformationFormat="1" \
                Detect64BitPortabilityProblems="true" \

                $uses_asm && tag Tool Name="YASM"  IncludePaths="$incs" Debug="1"
            ;;
        esac

        case "$proj_kind" in
            exe)
                case "$target" in
                    x86*) tag Tool \
                          Name="VCLinkerTool" \
                          AdditionalDependencies="$debug_libs \$(NoInherit)" \
                          AdditionalLibraryDirectories="$libdirs" \
                          GenerateDebugInformation="true" \
                          ProgramDatabaseFile="\$(OutDir)/${name}.pdb" \

                    ;;
                    arm*|iwmmx*)
                        case "$name" in
                            obj_int_extract) tag Tool \
                                Name="VCLinkerTool" \
                                OutputFile="${name}.exe" \
                                GenerateDebugInformation="true"
                            ;;
                            *) tag Tool \
                                Name="VCLinkerTool" \
                                AdditionalDependencies="$debug_libs" \
                                OutputFile="\$(OutDir)/${name}.exe" \
                                LinkIncremental="2" \
                                AdditionalLibraryDirectories="${libdirs};&quot;..\lib/$plat_no_ws&quot;" \
                                DelayLoadDLLs="\$(NOINHERIT)" \
                                GenerateDebugInformation="true" \
                                ProgramDatabaseFile="\$(OutDir)/${name}.pdb" \
                                SubSystem="9" \
                                StackReserveSize="65536" \
                                StackCommitSize="4096" \
                                EntryPointSymbol="mainWCRTStartup" \
                                TargetMachine="3"
                            ;;
                        esac
                     ;;
                 esac
            ;;
            lib)
                case "$target" in
                      arm*|iwmmx*) tag Tool \
                                    Name="VCLibrarianTool" \
                                    AdditionalOptions=" /subsystem:windowsce,4.20 /machine:ARM" \
                                    OutputFile="\$(OutDir)/${name}.lib" \
                                ;;
                                *) tag Tool \
                                    Name="VCLibrarianTool" \
                                    OutputFile="\$(OutDir)/${name}${lib_sfx}d.lib" \
                                ;;
                esac
            ;;
            dll) tag Tool \
                 Name="VCLinkerTool" \
                AdditionalDependencies="\$(NoInherit)" \
                LinkIncremental="2" \
                GenerateDebugInformation="true" \
                AssemblyDebug="1" \
                TargetMachine="1" \
                      $link_opts
        esac

        if [ "$target" == "armv6-wince-vs8" ] || [ "$target" == "armv5te-wince-vs8" ] || [ "$target" == "iwmmxt-wince-vs8" ] || [ "$target" == "iwmmxt2-wince-vs8" ];then
            case "$name" in
                vpx_decoder) tag DeploymentTool \
                             ForceDirty="-1" \
                             RegisterOutput="0"
                                ;;
                example|xma) tag DeploymentTool \
                             ForceDirty="-1" \
                             RegisterOutput="0"
                             tag DebuggerTool \
                             Arguments="${ARGU}"
                                ;;
            esac
        fi
        close_tag Configuration

        open_tag  Configuration \
                      Name="Release|$plat" \
                      OutputDirectory="\$(SolutionDir)$plat_no_ws/\$(ConfigurationName)" \
                      IntermediateDirectory="$plat_no_ws/\$(ConfigurationName)/${name}" \
                      ConfigurationType="$vs_ConfigurationType" \
                      CharacterSet="1" \
                      WholeProgramOptimization="0"

        if [ "$target" == "armv6-wince-vs8" ] || [ "$target" == "armv5te-wince-vs8" ] || [ "$target" == "iwmmxt-wince-vs8" ] || [ "$target" == "iwmmxt2-wince-vs8" ];then
            case "$name" in
                vpx_decoder) tag Tool \
                                     Name="VCPreBuildEventTool" \
                                     CommandLine="call obj_int_extract.bat \$(ConfigurationName)"
                             tag Tool \
                                     Name="VCMIDLTool" \
                                     TargetEnvironment="1"
                             tag Tool \
                                             Name="VCCLCompilerTool" \
                                             ExecutionBucket="7" \
                                             Optimization="2" \
                                             FavorSizeOrSpeed="1" \
                                             AdditionalIncludeDirectories="$incs" \
                                             PreprocessorDefinitions="NDEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES);WINCE;_LIB;\$(ARCHFAM);\$(_ARCHFAM_);_UNICODE;UNICODE;" \
                                             RuntimeLibrary="0" \
                                             BufferSecurityCheck="false" \
                                             UsePrecompiledHeader="0" \
                                             WarningLevel="3" \
                                             DebugInformationFormat="0" \
                                             CompileAs="1"
                             tag Tool \
                                             Name="VCResourceCompilerTool" \
                                             PreprocessorDefinitions="NDEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES)" \
                                             Culture="1033" \
                                             AdditionalIncludeDirectories="\$(IntDir)" \
                ;;
                example|xma) tag Tool \
                             Name="VCCLCompilerTool" \
                             ExecutionBucket="7" \
                             Optimization="2" \
                             FavorSizeOrSpeed="1" \
                             AdditionalIncludeDirectories="$incs" \
                             PreprocessorDefinitions="NDEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES);WINCE;_CONSOLE;\$(ARCHFAM);\$(_ARCHFAM_);_UNICODE;UNICODE;" \
                             RuntimeLibrary="0" \
                             BufferSecurityCheck="false" \
                             UsePrecompiledHeader="0" \
                             WarningLevel="3" \
                             DebugInformationFormat="0" \
                             CompileAs="1"
                             tag Tool \
                             Name="VCResourceCompilerTool" \
                             PreprocessorDefinitions="NDEBUG;_WIN32_WCE=\$(CEVER);UNDER_CE;\$(PLATFORMDEFINES)" \
                             Culture="1033" \
                             AdditionalIncludeDirectories="\$(IntDir)" \
                ;;
                obj_int_extract) tag Tool \
                             Name="VCCLCompilerTool" \
                             AdditionalIncludeDirectories="$incs" \
                             PreprocessorDefinitions="WIN32;NDEBUG;_CONSOLE" \
                             RuntimeLibrary="0" \
                             UsePrecompiledHeader="0" \
                             WarningLevel="3" \
                             Detect64BitPortabilityProblems="true" \
                             DebugInformationFormat="0" \
                ;;
            esac
        fi

    case "$target" in
        x86*) tag       Tool \
                      Name="VCCLCompilerTool" \
                      AdditionalIncludeDirectories="$incs" \
                      PreprocessorDefinitions="WIN32;NDEBUG;_CRT_SECURE_NO_WARNINGS;$defines" \
                      RuntimeLibrary="$release_runtime" \
                      UsePrecompiledHeader="0" \
                      WarningLevel="3" \
                      DebugInformationFormat="0" \
                      Detect64BitPortabilityProblems="true"

                $uses_asm && tag Tool Name="YASM"  IncludePaths="$incs"
                ;;
                esac

        case "$proj_kind" in
            exe)
                case "$target" in
                    x86*) tag Tool \
                                  Name="VCLinkerTool" \
                                  AdditionalDependencies="$libs \$(NoInherit)" \
                                  AdditionalLibraryDirectories="$libdirs" \
                    ;;
                    arm*|iwmmx*)
                        case "$name" in
                            obj_int_extract) tag Tool \
                                Name="VCLinkerTool" \
                                OutputFile="${name}.exe" \
                                LinkIncremental="1" \
                                GenerateDebugInformation="false" \
                                SubSystem="0" \
                                OptimizeReferences="0" \
                                EnableCOMDATFolding="0" \
                                TargetMachine="0"
                            ;;
                            *) tag Tool \
                                Name="VCLinkerTool" \
                                AdditionalDependencies="$libs" \
                                OutputFile="\$(OutDir)/${name}.exe" \
                                LinkIncremental="1" \
                                AdditionalLibraryDirectories="${libdirs};&quot;..\lib/$plat_no_ws&quot;" \
                                DelayLoadDLLs="\$(NOINHERIT)" \
                                GenerateDebugInformation="true" \
                                ProgramDatabaseFile="\$(OutDir)/${name}.pdb" \
                                SubSystem="9" \
                                StackReserveSize="65536" \
                                StackCommitSize="4096" \
                                OptimizeReferences="2" \
                                EnableCOMDATFolding="2" \
                                EntryPointSymbol="mainWCRTStartup" \
                                TargetMachine="3"
                            ;;
                        esac
                     ;;
                 esac
            ;;
        lib)
                case "$target" in
                      arm*|iwmmx*) tag Tool \
                                    Name="VCLibrarianTool" \
                                    AdditionalOptions=" /subsystem:windowsce,4.20 /machine:ARM" \
                                    OutputFile="\$(OutDir)/${name}.lib" \
                                ;;
                                *) tag Tool \
                                    Name="VCLibrarianTool" \
                                    OutputFile="\$(OutDir)/${name}${lib_sfx}.lib" \
                                ;;
                esac
        ;;
        dll) # note differences to debug version: LinkIncremental, AssemblyDebug
             tag Tool \
                      Name="VCLinkerTool" \
                      AdditionalDependencies="\$(NoInherit)" \
                      LinkIncremental="1" \
                      GenerateDebugInformation="true" \
                      TargetMachine="1" \
                      $link_opts
        esac

        if [ "$target" == "armv6-wince-vs8" ] || [ "$target" == "armv5te-wince-vs8" ] || [ "$target" == "iwmmxt-wince-vs8" ] || [ "$target" == "iwmmxt2-wince-vs8" ];then
            case "$name" in
                vpx_decoder) tag DeploymentTool \
                             ForceDirty="-1" \
                             RegisterOutput="0"
                ;;
                example|xma) tag DeploymentTool \
                             ForceDirty="-1" \
                             RegisterOutput="0"
                             tag DebuggerTool \
                             Arguments="${ARGU}"
                                ;;
            esac
        fi

        close_tag Configuration
    done
    close_tag Configurations

    open_tag  Files
    generate_filter srcs   "Source Files"   "cpp;c;cc;cxx;def;odl;idl;hpj;bat;asm;asmx"
    generate_filter hdrs   "Header Files"   "h;hpp;hxx;hm;inl;inc;xsd"
    generate_filter resrcs "Resource Files" "rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe;resx;tiff;tif;png;wav"
    generate_filter resrcs "Build Files"    "mk"
    close_tag Files

    tag       Globals
    close_tag VisualStudioProject

    # This must be done from within the {} subshell
    echo "Ignored files list (${#file_list[@]} items) is:" >&2
    for f in "${file_list[@]}"; do
        echo "    $f" >&2
    done
}

generate_vcproj |
    sed  -e '/"/s;\([^ "]\)/;\1\\;g' > ${outfile}

exit
<!--
TODO: Add any files not captured by filters.
                <File
                        RelativePath=".\ReadMe.txt"
                        >
                </File>
-->
