#!/bin/sh

script_path=~/Dev/sandbox/libvpx/scripts

av1_code=~/Dev/av1d
vp9_code=~/Dev/vp9d

log_path=~/Dev/log
date_str=`date +%H:%M_%b_%d_%Y`
log_file=rep_$date_str.txt

$script_path/nightly_speed.sh $av1_code > $log_path/$log_file 2>&1
$script_path/vp9_nightly_speed.sh $vp9_code >> $log_path/$log_file 2>&1

users=luoyi
host_name=`hostname`
sender=luoyi
cc_list=yunqingwang
#--cc=$cc_list

sendgmr --to=$users --cc=$cc_list --subject="Codec Daily Report" --from=$sender --reply_to=$sender < $log_path/$log_file
