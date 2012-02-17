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
  --require-EXT Require support for EXT extensions
  --sym=SYMBOL  Unique symbol to use for RTCD initialization function
  --config=FILE File with CONFIG_FOO=yes lines to parse
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
    --require-*) REQUIRES="${REQUIRES}${opt#--require-} ";;
    --sym) die_argument_required;;
    --sym=*) symbol=${optval};;
    --config=*) config_file=${optval};;
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

# Import the configuration
[ -f "$config_file" ] && eval $(grep CONFIG_ "$config_file")

#
# Routines for the RTCD DSL to call
#
prototype() {
  local rtyp
  case "$1" in
    unsigned) rtyp="$1 "; shift;;
  esac
  rtyp="${rtyp}$1"
  local fn="$2"
  local args="$3"

  eval "${2}_rtyp='$rtyp'"
  eval "${2}_args='$3'"
  ALL_FUNCS="$ALL_FUNCS $fn"
  specialize $fn c
}

specialize() {
  local fn="$1"
  shift
  for opt in "$@"; do
    eval "${fn}_${opt}=${fn}_${opt}"
  done
}

require() {
  for fn in $ALL_FUNCS; do
    for opt in "$@"; do
      local ofn=$(eval "echo \$${fn}_${opt}")
      [ -z "$ofn" ] && continue

      # if we already have a default, then we can disable it, as we know
      # we can do better.
      local best=$(eval "echo \$${fn}_default")
      local best_ofn=$(eval "echo \$${best}")
      [ -n "$best" ] && [ "$best_ofn" != "$ofn" ] && eval "${best}_link=false"
      eval "${fn}_default=${fn}_${opt}"
      eval "${fn}_${opt}_link=true"
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
    local rtyp="$(eval "echo \$${fn}_rtyp")"
    local args="$(eval "echo \"\$${fn}_args\"")"
    local dfn="$(eval "echo \$${fn}_default")"
    dfn=$(eval "echo \$${dfn}")
    for opt in "$@"; do
      local ofn=$(eval "echo \$${fn}_${opt}")
      [ -z "$ofn" ] && continue
      local link=$(eval "echo \$${fn}_${opt}_link")
      [ "$link" = "false" ] && continue
      n="${n}x"
    done
    if [ "$n" = "x" ]; then
      eval "${fn}_indirect=false"
    else
      eval "${fn}_indirect=true"
    fi
  done
}

declare_function_pointers() {
  for fn in $ALL_FUNCS; do
    local rtyp="$(eval "echo \$${fn}_rtyp")"
    local args="$(eval "echo \"\$${fn}_args\"")"
    local dfn="$(eval "echo \$${fn}_default")"
    dfn=$(eval "echo \$${dfn}")
    for opt in "$@"; do
      local ofn=$(eval "echo \$${fn}_${opt}")
      [ -z "$ofn" ] && continue
      echo "$rtyp ${ofn}($args);"
    done
    if [ "$(eval "echo \$${fn}_indirect")" = "false" ]; then
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
    local rtyp="$(eval "echo \$${fn}_rtyp")"
    local args="$(eval "echo \"\$${fn}_args\"")"
    local dfn="$(eval "echo \$${fn}_default")"
    dfn=$(eval "echo \$${dfn}")
    if $(eval "echo \$${fn}_indirect"); then
      echo "    $fn = $dfn;"
      for opt in "$@"; do
        local ofn=$(eval "echo \$${fn}_${opt}")
        [ -z "$ofn" ] && continue
        [ "$ofn" = "$dfn" ] && continue;
        local link=$(eval "echo \$${fn}_${opt}_link")
        [ "$link" = "false" ] && continue
        local cond="$(eval "echo \$have_${opt}")"
        echo "    if (${cond}) $fn = $ofn;"
      done
    fi
    echo
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
  local outfile_basename=$(basename ${symbol:-rtcd.h})
  local include_guard=$(echo $outfile_basename | tr '[a-z]' '[A-Z]' | tr -c '[A-Z]' _)
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

  # Assign the helper variable for each enabled extension
  for opt in $ALL_ARCHS; do
    local uc=$(echo $opt | tr '[a-z]' '[A-Z]')
    eval "have_${opt}=\"flags & HAS_${uc}\""
  done

  cat <<EOF
$(common_top)
void ${symbol:-rtcd}(void);

#ifdef RTCD_C
#include "vpx_ports/x86.h"
void ${symbol:-rtcd}(void)
{
    int flags = x86_simd_caps();

    (void)flags;

$(set_function_pointers c $ALL_ARCHS)
}
#endif
$(common_bottom)
EOF
}

arm() {
  determine_indirection c $ALL_ARCHS

  # Assign the helper variable for each enabled extension
  for opt in $ALL_ARCHS; do
    local uc=$(echo $opt | tr '[a-z]' '[A-Z]')
    eval "have_${opt}=\"flags & HAS_${uc}\""
  done

  cat <<EOF
$(common_top)
#include "vpx_config.h"

void ${symbol:-rtcd}(void);

#ifdef RTCD_C
#include "vpx_ports/arm.h"
void ${symbol:-rtcd}(void)
{
    int flags = arm_cpu_caps();

    (void)flags;

$(set_function_pointers c $ALL_ARCHS)
}
#endif
$(common_bottom)
EOF
}


unoptimized() {
  determine_indirection c
  cat <<EOF
$(common_top)
#include "vpx_config.h"

void ${symbol:-rtcd}(void);

#ifdef RTCD_C
void ${symbol:-rtcd}(void)
{
$(set_function_pointers c)
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
    ALL_ARCHS=$(filter mmx sse sse2 sse3 ssse3 sse4_1)
    x86
    ;;
  x86_64)
    ALL_ARCHS=$(filter mmx sse sse2 sse3 ssse3 sse4_1)
    REQUIRES=${REQUIRES:-mmx sse sse2}
    require $(filter $REQUIRES)
    x86
    ;;
  armv5te)
    ALL_ARCHS=$(filter edsp)
    arm
    ;;
  armv6)
    ALL_ARCHS=$(filter edsp media)
    arm
    ;;
  armv7)
    ALL_ARCHS=$(filter edsp media neon)
    arm
    ;;
  *)
    unoptimized
    ;;
esac
