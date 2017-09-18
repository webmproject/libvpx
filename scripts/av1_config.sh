#!/bin/sh
#set -x

if [ "$2" == "dbg" ]; then
  debug="--enable-debug --disable-optimizations"
else
  debug=
fi
tool=--enable-$1

platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom
test_dir=~/Dev/info
script_dir=~/Dev/sandbox/libvpx/scripts


common="--disable-docs --enable-experimental"

bd_config="--enable-lowbitdepth"
#--disable-highbitdepth"

. $script_dir/disabled_list.sh

../$libsrc/configure $debug $bd_config $common $disabled $tool > /dev/null
if [ $? -ne 0 ]; then
  echo "Error: configure fails!" > $test_dir/error_config.txt
  exit 1
fi
