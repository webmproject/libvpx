#!/bin/sh
self=$0

usage() {
  cat <<EOF >&2
Usage: $self [options] FILE

Reads the Run Time CPU Detections definitions from FILE and generates a
C header file on stdout.

Options:
  --arch=ARCH   Architecture to generate defs for (required)
  --disable-EXT Disable support for EXT extensions
  --sym=SYMBOL  Unique symbol to use for RTCD initialization function

EOF
  exit 1
}

die() {
  echo "$@" >&2
  exit 1
}

die_argument_required() {
  die "Option $opt requires argument"
}

for opt; do
  optval="${opt#*=}"
  case "$opt" in
    --arch) die_argument_required;;
    --arch=*) arch=${optval};;
    --disable-*) eval "disable_${opt#--disable-}=true";;
    --sym) die_argument_required;;
    --sym=*) symbol=${optval};;
    --rtcd=*) CONFIG_RUNTIME_CPU_DETECT=${optval};;
    -h|--help)
      usage
      ;;
    -*)
      die "Unrecognized option: ${opt%%=*}"
      ;;
    *)
      defs_file="$defs_file $opt"
      ;;
  esac
  shift
done
for f in $defs_file; do [ -f "$f" ] || usage; done
[ -n "$arch" ] || usage

#
# Routines for the RTCD DSL to call
#
prototype() {
  local rtyp="$1"
  local fn="$2"
  local args="$3"

  eval "${2}_rtyp='$1'"
  eval "${2}_args='$3'"
  ALL_FUNCS="$ALL_FUNCS $fn"
}

specialize() {
  local fn="$1"
  shift
  for opt in c "$@"; do
    eval "${fn}_${opt}=${fn}_${opt}"
  done
}

require() {
  for fn in $ALL_FUNCS; do
    for opt in "$@"; do
      local ofn=$(eval "echo \$${fn}_${opt}")
      [ -z "$ofn" ] && continue

      # if we already have a default, then we can undefine it, as we know
      # we can do better.
      local best=$(eval "echo \$${fn}_default")
      [ -n "$best" ] && eval "unset $best"
      eval "${fn}_default=${fn}_${opt}"
    done
  done
}

forward_decls() {
  ALL_FORWARD_DECLS="$ALL_FORWARD_DECLS $1"
}

#
# Include the user's directives
#
for f in $defs_file; do
  . $f
done

#
# Process the directives according to the command line
#
process_forward_decls() {
  for fn in $ALL_FORWARD_DECLS; do
    eval $fn
  done
}

determine_indirection() {
  [ "$CONFIG_RUNTIME_CPU_DETECT" = "yes" ] || require $ALL_ARCHS
  for fn in $ALL_FUNCS; do
    local n=""
    local rtyp=$(eval "echo \$${fn}_rtyp")
    local args=$(eval "echo \$${fn}_args")
    local dfn=$(eval "echo \$${fn}_default")
    dfn=$(eval "echo \$${dfn}")
    for opt in "$@"; do
      local ofn=$(eval "echo \$${fn}_${opt}")
      [ -z "$ofn" ] && continue
      n="${n}x"
    done
    if [ "$n" = "x" ]; then
      eval "${fn}_indirect=false"
    else
      eval "${fn}_indirect=true"
    fi
    echo
  done
}

declare_function_pointers() {
  for fn in $ALL_FUNCS; do
    local n=""
    local rtyp=$(eval "echo \$${fn}_rtyp")
    local args=$(eval "echo \$${fn}_args")
    local dfn=$(eval "echo \$${fn}_default")
    dfn=$(eval "echo \$${dfn}")
    for opt in "$@"; do
      local ofn=$(eval "echo \$${fn}_${opt}")
      [ -z "$ofn" ] && continue
      n="${n}x"
      echo "$rtyp ${ofn}($args);"
    done
    if [ "$n" = "x" ]; then
      echo "#define ${fn} ${dfn}"
    else
      echo "RTCD_EXTERN $rtyp (*${fn})($args);"
    fi
    echo
  done
}

set_function_pointers() {
  for fn in $ALL_FUNCS; do
    local n=""
    local rtyp=$(eval "echo \$${fn}_rtyp")
    local args=$(eval "echo \$${fn}_args")
    local dfn=$(eval "echo \$${fn}_default")
    dfn=$(eval "echo \$${dfn}")
    if $(eval "echo \$${fn}_indirect"); then
      echo "    $fn = $dfn;"
      for opt in "$@"; do
        local ofn=$(eval "echo \$${fn}_${opt}")
        [ -z "$ofn" ] && continue
        [ "$ofn" = "$dfn" ] && continue;
        echo "    if (have_${opt}) $fn = $ofn;"
      done
      echo
    fi
  done
}

filter() {
  local filtered
  for opt in "$@"; do
    [ -z $(eval "echo \$disable_${opt}") ] && filtered="$filtered $opt"
  done
  echo $filtered
}

#
# Helper functions for generating the arch specific RTCD files
#
common_top() {
  local outfile_basename=$(basename ${outfile:-rtcd.h})
  local include_guard=$(echo -n $outfile_basename | tr [a-z] [A-Z] | tr -c [A-Z] _)
  cat <<EOF
#ifndef ${include_guard}
#define ${include_guard}

#ifdef RTCD_C
#define RTCD_EXTERN
#else
#define RTCD_EXTERN extern
#endif

$(process_forward_decls)

$(declare_function_pointers c $ALL_ARCHS)
EOF
}

common_bottom() {
  cat <<EOF
#endif
EOF
}

x86() {
  determine_indirection c $ALL_ARCHS
  cat <<EOF
$(common_top)
void ${symbol:-rtcd}(void);

#ifdef RTCD_C
#include "vpx_ports/x86.h"
void ${symbol:-rtcd}(void)
{
    int flags = x86_simd_caps();
EOF

  # Write out the helper variable for each enabled extension
  for opt in $ALL_ARCHS; do
    local uc=$(echo -n $opt | tr [a-z] [A-Z])
    echo "    int have_${opt} = flags & HAS_${uc};"
  done
  cat <<EOF

$(set_function_pointers c $ALL_ARCHS)
}
#endif
$(common_bottom)
EOF
}

arm() {
  determine_indirection c $ALL_ARCHS
  cat <<EOF
$(common_top)
#include "vpx_config.h"
#include "vp8/decoder/onyxd_int.h"

void ${symbol:-rtcd}(VP8D_COMP *pbi);

#ifdef RTCD_C
void ${symbol:-rtcd}(VP8D_COMP *pbi)
{
#if CONFIG_RUNTIME_CPU_DETECT
    int flags = pbi->common.rtcd.flags;

    int have_v5te = flags & HAS_EDSP;
    int have_v6 = flags & HAS_MEDIA;
    int have_neon = flags & HAS_NEON;
#endif

$(set_function_pointers c $ALL_ARCHS)
}
#endif
$(common_bottom)
EOF
}

#
# Main Driver
#
require c
case $arch in
  x86)
    ALL_ARCHS=$(filter mmx sse sse2 sse3 sse4_1)
    x86
    ;;
  x86_64)
    ALL_ARCHS=$(filter mmx sse sse2 sse3 sse4_1)
    require $(filter mmx sse sse2)
    x86
    ;;
  armv5te)
    ALL_ARCHS=$(filter v5te)
    arm
    ;;
  armv6)
    ALL_ARCHS=$(filter v5te v6)
    arm
    ;;
  armv7)
    ALL_ARCHS=$(filter v5te v6 neon)
    arm
    ;;
  *)
    die "Unrecognized architecture: $arch"
esac
