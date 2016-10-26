##
## Copyright (c) 2016, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
cmake_minimum_required(VERSION 3.2)

include("${AOM_ROOT}/build/cmake/compiler_flags.cmake")

include(FindGit)

# Defaults for every libaom configuration variable.
set(RESTRICT)
set(INLINE)
set(ARCH_ARM 0)
set(ARCH_MIPS 0)
set(ARCH_X86 0)
set(ARCH_X86_64 0)
set(HAVE_EDSP 0)
set(HAVE_MEDIA 0)
set(HAVE_NEON 0)
set(HAVE_NEON_ASM 0)
set(HAVE_MIPS32 0)
set(HAVE_DSPR2 0)
set(HAVE_MSA 0)
set(HAVE_MIPS64 0)
set(HAVE_MMX 0)
set(HAVE_SSE 0)
set(HAVE_SSE2 0)
set(HAVE_SSE3 0)
set(HAVE_SSSE3 0)
set(HAVE_SSE4_1 0)
set(HAVE_AVX 0)
set(HAVE_AVX2 0)
set(HAVE_AOM_PORTS 0)
set(HAVE_PTHREAD_H 0)
set(HAVE_UNISTD_H 0)
set(CONFIG_DEPENDENCY_TRACKING 1)
set(CONFIG_EXTERNAL_BUILD 0)
set(CONFIG_INSTALL_DOCS 0)
set(CONFIG_INSTALL_BINS 0)
set(CONFIG_INSTALL_LIBS 0)
set(CONFIG_INSTALL_SRCS 0)
set(CONFIG_USE_X86INC 0)
set(CONFIG_DEBUG 0)
set(CONFIG_GPROF 0)
set(CONFIG_GCOV 0)
set(CONFIG_RVCT 0)
set(CONFIG_GCC 0)
set(CONFIG_MSVS 0)
set(CONFIG_PIC 0)
set(CONFIG_BIG_ENDIAN 0)
set(CONFIG_CODEC_SRCS 0)
set(CONFIG_DEBUG_LIBS 0)
set(CONFIG_DEQUANT_TOKENS 0)
set(CONFIG_DC_RECON 0)
set(CONFIG_RUNTIME_CPU_DETECT 0)
set(CONFIG_MULTITHREAD 0)
set(CONFIG_INTERNAL_STATS 0)
set(CONFIG_AV1_ENCODER 1)
set(CONFIG_AV1_DECODER 1)
set(CONFIG_AV1 1)
set(CONFIG_ENCODERS 1)
set(CONFIG_DECODERS 1)
set(CONFIG_STATIC_MSVCRT 0)
set(CONFIG_SPATIAL_RESAMPLING 1)
set(CONFIG_REALTIME_ONLY 0)
set(CONFIG_ONTHEFLY_BITPACKING 0)
set(CONFIG_ERROR_CONCEALMENT 0)
set(CONFIG_SHARED 0)
set(CONFIG_STATIC 1)
set(CONFIG_SMALL 0)
set(CONFIG_OS_SUPPORT 0)
set(CONFIG_UNIT_TESTS 0)
set(CONFIG_WEBM_IO 0)
set(CONFIG_LIBYUV 0)
set(CONFIG_ACCOUNTING 0)
set(CONFIG_DECODE_PERF_TESTS 0)
set(CONFIG_ENCODE_PERF_TESTS 0)
set(CONFIG_MULTI_RES_ENCODING 0)
set(CONFIG_TEMPORAL_DENOISING 1)
set(CONFIG_COEFFICIENT_RANGE_CHECKING 0)
set(CONFIG_AOM_HIGHBITDEPTH 0)
set(CONFIG_EXPERIMENTAL 0)
set(CONFIG_SIZE_LIMIT 0)
set(CONFIG_AOM_QM 0)
set(CONFIG_SPATIAL_SVC 0)
set(CONFIG_FP_MB_STATS 0)
set(CONFIG_EMULATE_HARDWARE 0)
set(CONFIG_CLPF 0)
set(CONFIG_DERING 0)
set(CONFIG_REF_MV 0)
set(CONFIG_SUB8X8_MC 0)
set(CONFIG_EXT_INTRA 0)
set(CONFIG_EXT_INTERP 0)
set(CONFIG_EXT_TX 0)
set(CONFIG_MOTION_VAR 0)
set(CONFIG_EXT_REFS 0)
set(CONFIG_EXT_COMPOUND 0)
set(CONFIG_SUPERTX 0)
set(CONFIG_ANS 0)
set(CONFIG_EC_MULTISYMBOL 0)
set(CONFIG_DAALA_EC 0)
set(CONFIG_PARALLEL_DEBLOCKING 0)
set(CONFIG_CB4X4 0)
set(CONFIG_PALETTE 0)
set(CONFIG_FRAME_SIZE 0)
set(CONFIG_FILTER_7BIT 0)
set(CONFIG_DELTA_Q 0)
set(CONFIG_ADAPT_SCAN 0)
set(CONFIG_BITSTREAM_DEBUG 0)
set(CONFIG_TILE_GROUPS 0)
set(CONFIG_EC_ADAPT 0)

# TODO(tomfinegan): consume trailing whitespace after configure_file().
configure_file("${AOM_ROOT}/build/cmake/aom_config.h.cmake"
               "${CMAKE_CURRENT_BINARY_DIR}/aom_config.h")

# Read the current git hash.
find_package(Git)
if (GIT_FOUND)
  # TODO(tomfinegan): Make this smart enough to write a proper version string
  # when in a repo that is on a label and clean.
  # TODO(tomfinegan): In addition to the one above, also make this a custom
  # build rule so users don't have to re-run cmake to create accurately
  # versioned cmake builds.
  execute_process(COMMAND ${GIT_EXECUTABLE}
                  --git-dir=${AOM_ROOT}/.git rev-parse HEAD
                  OUTPUT_VARIABLE AOM_GIT_HASH)
  # Consume the newline at the end of the git output.
  string(STRIP ${AOM_GIT_HASH} AOM_GIT_HASH)
else ()
  set(AOM_GIT_HASH)
endIf ()

# TODO(tomfinegan): An alternative to dumping the configure command line to
# aom_config.c is needed in cmake. Normal cmake generation runs do not make the
# command line available in the cmake script. For now, we just set the variable
# to the following. The configure_file() command will expand the message in
# aom_config.c.
# Note: This message isn't strictly true. When cmake is run in script mode (with
# the -P argument), CMAKE_ARGC and CMAKE_ARGVn are defined (n = 0 through
# n = CMAKE_ARGC become valid). Normal cmake generation runs do not make the
# information available.
set(AOM_CMAKE_CONFIG "cmake")
configure_file("${AOM_ROOT}/build/cmake/aom_config.c.cmake"
               "${CMAKE_CURRENT_BINARY_DIR}/aom_config.c")
