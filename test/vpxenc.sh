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

TEST_FRAMES=10

# Environment check: Make sure input is available.
vpxenc_verify_environment() {
  if [ ! -e "${YUV_RAW_INPUT}" ]; then
    echo "The file ${YUV_RAW_INPUT##*/} must exist in LIBVPX_TEST_DATA_PATH."
    return 1
  fi
}

vpxenc_can_encode_vp8() {
  if [ "$(vpxenc_available)" = "yes" ] && \
     [ "$(vp8_encode_available)" = "yes" ]; then
    echo yes
  fi
}

vpxenc_can_encode_vp9() {
  if [ "$(vpxenc_available)" = "yes" ] && \
     [ "$(vp9_encode_available)" = "yes" ]; then
    echo yes
  fi
}

# Echoes yes to stdout when vpxenc exists according to vpx_tool_available().
vpxenc_available() {
  [ -n $(vpx_tool_available vpxenc) ] && echo yes
}

# Wrapper function for running vpxenc. Positional parameters are interpreted as
# follows:
#   1 - codec name
#   2 - input width
#   3 - input height
#   4 - number of frames to encode
#   5 - path to input file
#   6 - path to output file
#       Note: The output file path must end in .ivf to output an IVF file.
#   7 - extra flags
#       Note: Extra flags currently supports a special case: when set to "-"
#             input is piped to vpxenc via cat.
vpxenc() {
  local encoder="${LIBVPX_BIN_PATH}/vpxenc${VPX_TEST_EXE_SUFFIX}"
  local codec="${1}"
  local width=${2}
  local height=${3}
  local frames=${4}
  local input=${5}
  local output="${VPX_TEST_OUTPUT_DIR}/${6}"
  local extra_flags=${7}

  # Because --ivf must be within the command line to get IVF from vpxenc.
  if echo "${output}" | egrep -q 'ivf$'; then
    use_ivf=--ivf
  else
    unset use_ivf
  fi

  if [ "${extra_flags}" = "-" ]; then
    pipe_input=yes
    extra_flags=${8}
  else
    unset pipe_input
  fi

  if [ -z "${pipe_input}" ]; then
    eval "${VPX_TEST_PREFIX}" "${encoder}" --codec=${codec} --width=${width} \
        --height=${height} --limit=${frames} ${use_ivf} ${extra_flags} \
        --output="${output}" "${input}" ${devnull}
  else
    cat "${input}" \
        | eval "${VPX_TEST_PREFIX}" "${encoder}" --codec=${codec} \
            --width=${width} --height=${height} --limit=${frames} ${use_ivf} \
            ${extra_flags} --output="${output}" - ${devnull}
  fi

  if [ ! -e "${output}" ]; then
    # Return non-zero exit status: output file doesn't exist, so something
    # definitely went wrong.
    return 1
  fi
}

vpxenc_vp8_ivf() {
  if [ "$(vpxenc_can_encode_vp8)" = "yes" ]; then
    vpxenc vp8 ${YUV_RAW_INPUT_WIDTH} ${YUV_RAW_INPUT_HEIGHT} ${TEST_FRAMES} \
        "${YUV_RAW_INPUT}" vp8.ivf
  fi
}

vpxenc_vp8_ivf_pipe_input() {
  if [ "$(vpxenc_can_encode_vp8)" = "yes" ]; then
    vpxenc vp8 ${YUV_RAW_INPUT_WIDTH} ${YUV_RAW_INPUT_HEIGHT} ${TEST_FRAMES} \
        "${YUV_RAW_INPUT}" vp8.ivf -
  fi
}

vpxenc_vp8_webm() {
  if [ "$(vpxenc_can_encode_vp8)" = "yes" ] &&
     [ "$(webm_io_available)" = "yes" ] ; then
    vpxenc vp8 ${YUV_RAW_INPUT_WIDTH} ${YUV_RAW_INPUT_HEIGHT} ${TEST_FRAMES} \
        "${YUV_RAW_INPUT}" vp8.webm
  fi
}

vpxenc_vp9_ivf() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ]; then
    vpxenc vp9 ${YUV_RAW_INPUT_WIDTH} ${YUV_RAW_INPUT_HEIGHT} ${TEST_FRAMES} \
        "${YUV_RAW_INPUT}" vp9.ivf
  fi
}

vpxenc_vp9_webm() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ] &&
     [ "$(webm_io_available)" = "yes" ] ; then
    vpxenc vp9 ${YUV_RAW_INPUT_WIDTH} ${YUV_RAW_INPUT_HEIGHT} ${TEST_FRAMES} \
        "${YUV_RAW_INPUT}" vp9.webm
  fi
}

vpxenc_vp9_ivf_lossless() {
  if [ "$(vpxenc_can_encode_vp9)" = "yes" ]; then
    vpxenc vp9 ${YUV_RAW_INPUT_WIDTH} ${YUV_RAW_INPUT_HEIGHT} ${TEST_FRAMES} \
        "${YUV_RAW_INPUT}" vp9_lossless.ivf --lossless=1
  fi
}

vpxenc_tests="vpxenc_vp8_ivf
              vpxenc_vp8_webm
              vpxenc_vp8_ivf_pipe_input
              vpxenc_vp9_ivf
              vpxenc_vp9_webm
              vpxenc_vp9_ivf_lossless"

run_tests vpxenc_verify_environment "${vpxenc_tests}"
