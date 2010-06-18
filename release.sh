#!/bin/sh
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

for opt; do
    case $opt in
        --clean) clean=yes;;
        -j*) jopt=$opt;;
        *) echo "Unsupported option $opt"; exit 1;;
    esac
done

TAB="$(printf '\t')"
cat > release.mk << EOF
%\$(BUILD_SFX).tar.bz2: %/.done
${TAB}@echo "\$(subst .tar.bz2,,\$@): tarball"
${TAB}@cd \$(dir \$<); tar -cf - \$(subst .tar.bz2,,\$@) | bzip2 > ../\$@

%\$(BUILD_SFX).zip: %/.done
${TAB}@echo "\$(subst .zip,,\$@): zip"
${TAB}@rm -f \$@; cd \$(dir \$<); zip -rq ../\$@ \$(subst .zip,,\$@)

logs/%\$(BUILD_SFX).log.bz2: %/.done
${TAB}@echo "\$(subst .log.bz2,,\$(notdir \$@)): tarlog"
${TAB}@mkdir -p logs
${TAB}@cat \$< | bzip2 > \$@

%/.done:
${TAB}@mkdir -p \$(dir \$@)
${TAB}@echo "\$(dir \$@): configure \$(CONFIG_OPTS) \$(EXTRA_PATH)"
${TAB}@cd \$(dir \$@); export PATH=\$\$PATH\$(EXTRA_PATH); ../\$(SRC_ROOT)/configure \$(CONFIG_OPTS) >makelog.txt 2>&1
${TAB}@echo "\$(dir \$@): make"
${TAB}@cd \$(dir \$@); PATH=\$\$PATH\$(EXTRA_PATH) \$(MAKE) >>makelog.txt 2>&1
${TAB}@echo "\$(dir \$@): test install"
${TAB}@cd \$(dir \$@); PATH=\$\$PATH\$(EXTRA_PATH) \$(MAKE) install >>makelog.txt 2>&1
${TAB}@cd \$(dir \$@)/dist/build; PATH=\$\$PATH\$(EXTRA_PATH) \$(MAKE) >>makelog.txt 2>&1
${TAB}@echo "\$(dir \$@): install"
${TAB}@cd \$(dir \$@); PATH=\$\$PATH\$(EXTRA_PATH) \$(MAKE) install DIST_DIR=\$(TGT) >>makelog.txt 2>&1
${TAB}@touch \$@

#include release-deps.mk
EOF

#[ -f release-deps.mk ] || \
#    find ${self%/*} -name .git -prune -o -type f -print0 \
#    | xargs -0 -n1 echo \
#    | sed -e 's; ;\\ ;g' | awk '{print "$(TGT)/.done: "$0}' > release-deps.mk

build_config_list() {
    for codec in $CODEC_LIST; do
        for arch in $ARCH_LIST; do
            if [ -n "$OS_LIST" ]; then
                for os in $OS_LIST; do
                    CONFIGS="$CONFIGS vpx-${codec}-${arch}-${os}"
                done
            else
                CONFIGS="$CONFIGS vpx-${codec}-${arch}"
            fi
        done
    done
}

CODEC_LIST="vp8 vp8cx vp8dx"
case `uname` in
    Linux*)
        ARCH_LIST="x86 x86_64"
        OS_LIST="linux"
        build_config_list
        ARCH_LIST="armv5te armv6 armv7"
        OS_LIST="linux-gcc"

        ;;
    CYGWIN*)
        TAR_SFX=.zip
        for vs in vs7 vs8; do
            for arch in x86-win32 x86_64-win64; do
                for msvcrt in md mt; do
                    case $vs,$arch in
                        vs7,x86_64-win64) continue ;;
                    esac
                    ARCH_LIST="$ARCH_LIST ${arch}${msvcrt}-${vs}"
                done
            done
        done
        ;;
    Darwin*)
        ARCH_LIST="universal"
        OS_LIST="darwin8 darwin9"
        ;;
    sun_os*)
        ARCH_LIST="x86 x86_64"
        OS_LIST="solaris"
        ;;
esac
build_config_list

TAR_SFX=${TAR_SFX:-.tar.bz2}
ARM_TOOLCHAIN=/usr/local/google/csl-2009q3-67
for cfg in $CONFIGS; do
    full_cfg=$cfg
    cfg=${cfg#vpx-}
    opts=
    rm -f makelog.txt

    case $cfg in
        src-*)  opts="$opts --enable-codec-srcs"
                cfg=${cfg#src-}
                ;;
        eval-*) opts="$opts --enable-eval-limit"
                cfg=${cfg#src-}
                ;;
    esac

    case $cfg in
        #
        # Linux
        #
        *x86-linux)
            opts="$opts --target=x86-linux-gcc" ;;
        *x86_64-linux)
            opts="$opts --target=x86_64-linux-gcc" ;;
        *arm*-linux-gcc)
            armv=${cfg##*armv}
            armv=${armv%%-*}
            opts="$opts --target=armv${armv}-linux-gcc" ;;
        *arm*-linux-rvct)
            armv=${cfg##*armv}
            armv=${armv%%-*}
            opts="$opts --target=armv${armv}-linux-rvct"
            opts="$opts --libc=${ARM_TOOLCHAIN}/arm-none-linux-gnueabi/libc" ;;


        #
        # Windows
        #
        # need --enable-debug-libs for now until we're smarter about
        # building the debug/release from the customer installed
        # environment
        *-x86-win32*-vs*)
            opts="$opts --target=x86-win32-vs${cfg##*-vs} --enable-debug-libs";;
        *-x86_64-win64*-vs8)
            opts="$opts --target=x86_64-win64-vs8 --enable-debug-libs" ;;

        #
        # Darwin
        #
        *-universal-darwin*)
            opts="$opts --target=universal-darwin${cfg##*-darwin}-gcc" ;;

        #
        # Solaris
        #
        *x86-solaris)
            opts="$opts --target=x86-solaris-gcc" ;;
        *x86_64-solaris)
            opts="$opts --target=x86_64-solaris-gcc" ;;
    esac

    case $cfg in
        *x86-linux | *x86-solaris) opts="$opts --enable-pic" ;;
    esac

    case $cfg in
        *-win[36][24]mt*)  opts="$opts --enable-static-msvcrt" ;;
        *-win[36][24]md*)  opts="$opts --disable-static-msvcrt" ;;
    esac

    opts="$opts --disable-codecs"
    case $cfg in
        vp8*) opts="$opts --enable-vp8" ;;
    esac
    case $cfg in
        *cx-*) opts="${opts}-encoder" ;;
        *dx-*) opts="${opts}-decoder" ;;
    esac
    opts="$opts --enable-postproc"

    [ "x${clean}" = "xyes" ] \
        && rm -rf ${full_cfg}${BUILD_SFX}${TAR_SFX} \
        && rm -rf logs/${full_cfg}${BUILD_SFX}.log.bz2

    TGT=${full_cfg}${BUILD_SFX}
    BUILD_TARGETS="logs/${TGT}.log.bz2 ${TGT}${TAR_SFX}"
    echo "${BUILD_TARGETS}: CONFIG_OPTS=$opts" >>release.mk
    echo "${BUILD_TARGETS}: TGT=${TGT}" >>release.mk
    case $cfg in
        *-arm*-linux-*)
            echo "${BUILD_TARGETS}: EXTRA_PATH=:${ARM_TOOLCHAIN}/bin/" >>release.mk ;;
        *-vs7)
            echo "${BUILD_TARGETS}: EXTRA_PATH=:/cygdrive/c/Program\ Files/Microsoft\ Visual\ Studio\ .NET\ 2003/Common7/IDE" >>release.mk ;;
        *-vs8)
            echo "${BUILD_TARGETS}: EXTRA_PATH=:/cygdrive/c/Program\ Files/Microsoft\ Visual\ Studio\ 8/Common7/IDE" >>release.mk ;;
    esac
    MAKE_TGTS="$MAKE_TGTS ${TGT}${TAR_SFX} logs/${TGT}.log.bz2"
done


${MAKE:-make} ${jopt:--j3} -f release.mk  \
    SRC_ROOT=${self%/*} BUILD_SFX=${BUILD_SFX} ${MAKE_TGTS}
