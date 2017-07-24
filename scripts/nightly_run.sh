#!/bin/sh

script_path=~/Dev/sandbox/libvpx/scripts
codebase=~/Dev/av1d
log_path=~/Dev/log
date_str=`date +%H:%M_%b_%d_%Y`
log_file=rep_$date_str.txt

$script_path/nightly_speed.sh $codebase > $log_path/$log_file 2>&1

users=luoyi
host_name=`hostname`
sender=luoyi
cc_list=yunqingwang

sendgmr --to=$users --cc=$cc_list --subject="AV1 Build/Running Report" --from=$sender --reply_to=$sender < $log_path/$log_file
