#!/bin/sh

root_dir=$1
build_dir=$root_dir/release
script_dir=~/Dev/sandbox/libvpx/scripts
exp_tool=experimental

cd $build_dir
make clean > /dev/null
$script_dir/nightly_config.sh $exp_tool
make -j > /dev/null
if [ $? -ne 0 ]; then
  echo "AV1 build failed on experiment: " $exp_tool
  exit 1
fi

test_dir=~/Dev/nightly

cp -f ./aomenc $test_dir/.
cp -f ./aomdec $test_dir/.
