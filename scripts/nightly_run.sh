#!/bin/sh

script_path=~/Dev/sandbox/libvpx/scripts

av1_code=~/Dev/av1d
vp9_code=~/Dev/vp9d

log_path=~/Dev/log
#date_str=`date +%H:%M_%b_%d_%Y`
date_str=`date +%b_%d_%Y`
log_file=report_$date_str.txt

test_dir=~/Dev/nightly
rm $test_dir/*

# AV1
$script_path/sync_codebase.sh $av1_code/aom > $log_path/$log_file 2>&1

$script_path/av1_conf_build.sh $av1_code >> $log_path/$log_file 2>&1
$script_path/nightly_speed.sh $av1_code 0 >> $log_path/$log_file 2>&1
$script_path/nightly_speed.sh $av1_code 2 >> $log_path/$log_file 2>&1

# VP9
$script_path/sync_codebase.sh $vp9_code/libvpx >> $log_path/$log_file 2>&1

$script_path/vp9_conf_build.sh $vp9_code lowbitdepth >> $log_path/$log_file 2>&1
$script_path/vp9_nightly_speed.sh $vp9_code 0 >> $log_path/$log_file 2>&1

$script_path/vp9_conf_build.sh $vp9_code highbitdepth >> $log_path/$log_file 2>&1
$script_path/vp9_nightly_speed.sh $vp9_code 2 >> $log_path/$log_file 2>&1

users=luoyi
host_name=`hostname`
sender=luoyi
#cc_list="--cc=yunqingwang"

sendgmr --to=$users $cc_list --subject="Codec Daily Report" --from=$sender --reply_to=$sender < $log_path/$log_file
