#!/bin/sh
##
##  Copyright (c) 2014 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##
##
## This script generates 'VPX.framework'. An iOS app can encode and decode VPx
## video by including 'VPX.framework'.
##
## Run iosbuild.sh to create 'VPX.framework' in the current directory.
##
set -e
devnull='> /dev/null 2>&1'

BUILD_ROOT="_iosbuild"
DIST_DIR="_dist"
FRAMEWORK_DIR="VPX.framework"
HEADER_DIR="${FRAMEWORK_DIR}/Headers/vpx"
MAKE_JOBS=1
LIBVPX_SOURCE_DIR=$(dirname "$0")
LIPO=$(xcrun -sdk iphoneos${SDK} -find lipo)
ORIG_PWD="$(pwd)"
TARGETS="armv6-darwin-gcc
         armv7-darwin-gcc
         armv7s-darwin-gcc
         x86-iphonesimulator-gcc
         x86_64-iphonesimulator-gcc"

# This variable is set to the last dist dir used with make dist, and reused when
# populating the framework directory to get the path to the most recent
# includes.
TARGET_DIST_DIR=""

# List of library files passed to lipo.
LIBS=""

build_target() {
  local target="$1"
  local old_pwd="$(pwd)"

  vlog "***Building target: ${target}***"

  mkdir "${target}"
  cd "${target}"
  eval "../../${LIBVPX_SOURCE_DIR}/configure" --target="${target}" \
      --disable-docs ${devnull}
  export DIST_DIR
  eval make -j ${MAKE_JOBS} dist ${devnull}
  cd "${old_pwd}"

  vlog "***Done building target: ${target}***"
}

build_targets() {
  local targets="$1"
  local target

  # Clean up from previous build(s).
  rm -rf "${BUILD_ROOT}" "${FRAMEWORK_DIR}"

  # Create output dirs.
  mkdir -p "${BUILD_ROOT}"
  mkdir -p "${HEADER_DIR}"

  cd "${BUILD_ROOT}"

  for target in ${targets}; do
    build_target "${target}"
    TARGET_DIST_DIR="${BUILD_ROOT}/${target}/${DIST_DIR}"
    LIBS="${LIBS} ${TARGET_DIST_DIR}/lib/libvpx.a"
  done

  cd "${ORIG_PWD}"
}

cleanup() {
  cd "${ORIG_PWD}"

  if [ "${PRESERVE_BUILD_OUTPUT}" != "yes" ]; then
    rm -rf "${BUILD_ROOT}"
  fi
}

iosbuild_usage() {
cat << EOF
  Usage: ${0##*/} [arguments]
    --help: Display this message and exit.
    --jobs: Number of make jobs.
    --preserve-build-output: Do not delete the build directory.
    --show-build-output: Show output from each library build.
    --verbose: Output information about the environment and each stage of the
               build.
EOF
}

vlog() {
  if [ "${VERBOSE}" = "yes" ]; then
    echo "$@"
  fi
}

trap cleanup EXIT

# Parse the command line.
while [ -n "$1" ]; do
  case "$1" in
    --help)
      iosbuild_usage
      exit
      ;;
    --jobs)
      MAKE_JOBS="$2"
      shift
      ;;
    --preserve-build-output)
      PRESERVE_BUILD_OUTPUT=yes
      ;;
    --show-build-output)
      devnull=
      ;;
    --verbose)
      VERBOSE=yes
      ;;
    *)
      iosbuild_usage
      exit 1
      ;;
  esac
  shift
done

if [ "${VERBOSE}" = "yes" ]; then
cat << EOF
  BUILD_ROOT=${BUILD_ROOT}
  DIST_DIR=${DIST_DIR}
  FRAMEWORK_DIR=${FRAMEWORK_DIR}
  HEADER_DIR=${HEADER_DIR}
  MAKE_JOBS=${MAKE_JOBS}
  PRESERVE_BUILD_OUTPUT=${PRESERVE_BUILD_OUTPUT}
  LIBVPX_SOURCE_DIR=${LIBVPX_SOURCE_DIR}
  LIPO=${LIPO}
  ORIG_PWD=${ORIG_PWD}
  TARGETS="${TARGETS}"
EOF
fi

build_targets "${TARGETS}"

# Includes are identical for all platforms, and according to dist target
# behavior vpx_config.h and vpx_version.h aren't actually necessary for user
# apps built with libvpx. So, just copy the includes from the last target built.
# TODO(tomfinegan): The above is a lame excuse. Build common config/version
# includes that use the preprocessor to include the correct file.
cp -p "${TARGET_DIST_DIR}"/include/vpx/* "${HEADER_DIR}"
${LIPO} -create ${LIBS} -output ${FRAMEWORK_DIR}/VPX

vlog "Created fat library ${FRAMEWORK_DIR}/VPX containing:"
for lib in ${LIBS}; do
  vlog "  $(echo ${lib} | awk -F / '{print $2, $NF}')"
done
