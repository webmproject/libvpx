#!/bin/sh

platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom

bitdepth=--enable-highbitdepth
#bitdepth=
#debug=--enable-debug
debug=

if [ "$#" -eq 1 ]; then
  tool=$1
  ../$libsrc/configure --target=$platform $codec $bitdepth $debug --enable-experimental --enable-$tool --disable-unit-tests > /dev/null
else
  tool_1=$1
  tool_2=$2
  ../$libsrc/configure --target=$platform $codec $bitdepth --enable-experimental --enable-$tool_1 --enable-$tool_2 --disable-unit-tests --disable-daala-ec > /dev/null
fi
