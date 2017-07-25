#!/bin/sh

script_path=~/Dev/sandbox/libvpx/scripts

av1_code=~/Dev/av1w

log_path=~/Dev/log
date_str=`date +%b_%d_%Y`
log_file=weekly_report_$date_str.txt

test_dir=~/Dev/weekly
rm $test_dir/*

$script_path/sync_codebase.sh $av1_code/aom > $log_path/$log_file 2>&1
$script_path/weekly_speed.sh $av1_code 0 >> $log_path/$log_file 2>&1

users=luoyi
host_name=`hostname`
sender=luoyi
#cc_list="--cc=yunqingwang"

sendgmr --to=$users $cc_list --subject="Codec Weekly Report" --from=$sender --reply_to=$sender < $log_path/$log_file
