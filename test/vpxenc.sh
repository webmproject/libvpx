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
##  This file tests vpxenc using hantro_collage_w352h288.yuv as input. To add
##  new tests to this file, do the following:
##    1. Write a shell function (this is your test).
##    2. Add the function to vpxenc_tests (on a new line).
##
. $(dirname $0)/tools_common.sh

readonly TEST_FRAMES=10

# Environment check: Make sure input is available.
vpxenc_verify_environment() {
  if [ ! -e "${YUV_RAW_INPUT}" ]; then
    echo "The file ${YUV_RAW_INPUT##*/} must exist in LIBVPX_TEST_DATA_PATH."
    return 1
  fi
  if [ -z "$(vpx_tool_path vpxenc)" ]; then
    elog "vpxenc not found. It must exist in LIBVPX_BIN_PATH or its parent."
    return 1
  fi
}

vpxenc_can_encode_vp8() {
  if [ "$(vp8_encode_available)" = "yes" ]; then
    echo yes
  fi
}

vpxenc_can_encode_vp9() {
  if [ "$(vp9_encode_available)" = "yes" ]; then
    echo yes
  fi
}

# Wrapper function for running vpxenc with pipe input. Requires that
# LIBVPX_BIN_PATH points to the directory containing vpxenc. $1 is used as the
# input file path and shifted away. All remaining parameters are passed through
# to vpxenc.
vpxenc_pipe() {
  local readonly encoder="$(vpx_tool_path vpxenc)"
  local readonly input="$1"
  shift
  cat "${input}" | eval "${VPX_TEST_PREFIX}" "${encoder}" - "$@" ${devnull}
}

# Wrapper function for running vpxenc. Requires that LIBVPX_BIN_PATH points to
# the directory containing vpxenc. $1 one is used as the input file path and
# shifted away. All remaining parameters are passed through to vpxenc.
vpxenc() {
  local readonly encoder="$(vpx_tool_path vpxenc)"
  local readonly input="${1}"
  shift
  eval "${VPX_TEST_PREFIX}" "${encoder}" "$input" "$@" ${devnull}
}

vpxenc_vp8_ivf() {
  if [ "$(vpxenc_can_encode_vp8)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp8.ivf"
    vpxenc --codec=vp8 \
      --width="${YUV_RAW_INPUT_WIDTH}" \
      --height="${YUV_RAW_INPUT_HEIGHT}" \
      --limit="${TEST_FRAMES}" \
      --ivf \
      --output="${output}" \
      "${YUV_RAW_INPUT}"

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_vp8_ivf_piped_input() {
  if [ "$(vpxenc_can_encode_vp8)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp8_piped_input.ivf"
    cat "${YUV_RAW_INPUT}" \
      | vpxenc --codec=vp8 \
        --width="${YUV_RAW_INPUT_WIDTH}" \
        --height="${YUV_RAW_INPUT_HEIGHT}" \
        --limit="${TEST_FRAMES}" \
        --ivf \
        --output="${output}" \
        -

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_vp8_webm() {
  if [ "$(vpxenc_can_encode_vp8)" = "yes" ] && \
     [ "$(webm_io_available)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp8.webm"
    vpxenc --codec=vp8 \
      --width="${YUV_RAW_INPUT_WIDTH}" \
      --height="${YUV_RAW_INPUT_HEIGHT}" \
      --limit="${TEST_FRAMES}" \
      --output="${output}" \
      "${YUV_RAW_INPUT}"

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_vp9_ivf() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp9.ivf"
    vpxenc --codec=vp9 \
      --width="${YUV_RAW_INPUT_WIDTH}" \
      --height="${YUV_RAW_INPUT_HEIGHT}" \
      --limit="${TEST_FRAMES}" \
      --ivf \
      --test-decode=fatal \
      --output="${output}" \
      "${YUV_RAW_INPUT}"

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_vp9_webm() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ] && \
     [ "$(webm_io_available)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp9.webm"
    vpxenc --codec=vp9 \
      --width="${YUV_RAW_INPUT_WIDTH}" \
      --height="${YUV_RAW_INPUT_HEIGHT}" \
      --limit="${TEST_FRAMES}" \
      --test-decode=fatal \
      --output="${output}" \
      "${YUV_RAW_INPUT}"

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_vp9_ivf_lossless() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp9_lossless.ivf"
    vpxenc --codec=vp9 \
      --width="${YUV_RAW_INPUT_WIDTH}" \
      --height="${YUV_RAW_INPUT_HEIGHT}" \
      --limit="${TEST_FRAMES}" \
      --ivf \
      --output="${output}" \
      --lossless=1 \
      --test-decode=fatal \
      "${YUV_RAW_INPUT}"

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_vp9_ivf_minq0_maxq0() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ]; then
    local readonly output="${VPX_TEST_OUTPUT_DIR}/vp9_lossless_minq0_maxq0.ivf"
    vpxenc --codec=vp9 \
      --width="${YUV_RAW_INPUT_WIDTH}" \
      --height="${YUV_RAW_INPUT_HEIGHT}" \
      --limit="${TEST_FRAMES}" \
      --ivf \
      --output="${output}" \
      --min-q=0 \
      --max-q=0 \
      --test-decode=fatal \
      "${YUV_RAW_INPUT}"

    if [ ! -e "${output}" ]; then
      elog "Output file does not exist."
      return 1
    fi
  fi
}

vpxenc_tests="vpxenc_vp8_ivf
              vpxenc_vp8_webm
              vpxenc_vp8_ivf_piped_input
              vpxenc_vp9_ivf
              vpxenc_vp9_webm
              vpxenc_vp9_ivf_lossless
              vpxenc_vp9_ivf_minq0_maxq0"

run_tests vpxenc_verify_environment "${vpxenc_tests}"
