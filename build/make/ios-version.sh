#!/bin/sh
##
##  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##

if [ "$1" = "--enable-shared" ]; then
  # Shared library framework builds are only possible on iOS 8 and later.
  echo "8.0"
else
  echo "6.0"
fi
