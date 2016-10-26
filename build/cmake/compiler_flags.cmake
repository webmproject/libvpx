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

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

function (add_c_flag_if_supported c_flag)
  unset(C_FLAG_SUPPORTED CACHE)
  message("Checking C compiler flag support for: " ${c_flag})
  check_c_compiler_flag("${c_flag}" C_FLAG_SUPPORTED)
  if (C_FLAG_SUPPORTED)
    set(CMAKE_C_FLAGS "${c_flag} ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
  endif ()
endfunction ()

function (add_cxx_flag_if_supported cxx_flag)
  unset(CXX_FLAG_SUPPORTED CACHE)
  message("Checking CXX compiler flag support for: " ${cxx_flag})
  check_cxx_compiler_flag("${cxx_flag}" CXX_FLAG_SUPPORTED)
  if (CXX_FLAG_SUPPORTED)
    set(CMAKE_CXX_FLAGS "${cxx_flag} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
  endif ()
endfunction ()

function (add_compiler_flag_if_supported flag)
  add_c_flag_if_supported(${flag})
  add_cxx_flag_if_supported(${flag})
endfunction ()

# Set warning levels.
if (MSVC)
  add_compiler_flag_if_supported("/W3")
  # Disable MSVC warnings that suggest making code non-portable.
  add_compiler_flag_if_supported("/wd4996")
  if (ENABLE_WERROR)
    add_compiler_flag_if_supported("/WX")
  endif ()
else ()
  add_compiler_flag_if_supported("-Wall")
  add_compiler_flag_if_supported("-Wextra")
  add_compiler_flag_if_supported("-Wno-deprecated")
  add_compiler_flag_if_supported("-Wshorten-64-to-32")
  add_compiler_flag_if_supported("-Wnarrowing")
  if (ENABLE_WERROR)
    add_compiler_flag_if_supported("-Werror")
  endif ()
endif ()
