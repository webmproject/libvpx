#!/bin/sh

script_path=~/Dev/sandbox/libvpx/scripts
log_path=~/Dev/log
log_file=special_log_08_19.txt

$script_path/list_exp_speed.sh ~/Dev/av1w > $log_path/$log_file 2>&1

users=luoyi
sender=luoyi
cc_list="--cc=yunqingwang"

sendgmr --to=$users $cc_list --subject="Weekly Codec Report" --from=$sender --reply_to=$sender --body_file=$log_path/$log_file
