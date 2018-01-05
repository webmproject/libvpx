#!/bin/sh
#set -x

platform=x86_64-linux-gcc
codec=--enable-av1
libsrc=aom
test_dir=~/Dev/nightly
script_dir=~/Dev/sandbox/libvpx/scripts
tool=

common="--disable-unit-tests --disable-docs"
debug=

#. $script_dir/disabled_list.sh
disabled=

echo ../$libsrc/configure $common $debug $disabled $tool

../$libsrc/configure $common $debug $disabled $tool > /dev/null
if [ $? -ne 0 ]; then
  echo "Error: configure fails!" > $test_dir/aom_error_config.txt
  exit 1
fi
