#!/bin/sh

platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom
test_dir=~/Dev/nightly
script_dir=~/Dev/sandbox/libvpx/scripts

tool=--enable-$1
common="--disable-unit-tests --disable-docs --enable-experimental"
debug=

. $script_dir/disabled_list.sh

../$libsrc/configure $debug $common $disabled $tool > /dev/null
if [ $? -ne 0 ]; then
  echo "Error: configure fails!" > $test_dir/error_config.txt
  exit 1
fi
