#!/bin/sh

root_dir=$1
bitdepthflag=$2

build_dir=$root_dir/release

if [ "$bitdepthflag" == "highbitdepth" ]; then
  build_flag=--enable-vp9-highbitdepth
else
  build_flag=
fi

cd $build_dir
make clean > /dev/null

common_flag="--disable-unit-tests --disable-docs"

../libvpx/configure $common_flag $build_flag > /dev/null

make -j > /dev/null
if [ $? -ne 0 ]; then
  echo "VP9 build failed"
  exit 1
fi

test_dir=~/Dev/nightly

cp -f ./vpxenc $test_dir/.
cp -f ./vpxdec $test_dir/.
