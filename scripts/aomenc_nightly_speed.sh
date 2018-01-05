#!/bin/sh
#set -x

script_path=~/Dev/sandbox/libvpx/scripts

av1_code=~/Dev/av1d
log_path=~/Dev/log

date_str=`date -d tomorrow +%b_%d_%Y`

log_file=aomenc_$date_str.txt
html_log_file=aomenc_$date_str.html

prev_date_str=`date +%b_%d_%Y`
prev_log_file=aomenc_$prev_date_str.txt

test_dir=~/Dev/nightly
rm $test_dir/*

$script_path/gen_html_header.sh > $log_path/$html_log_file

echo "<p>" >> $log_path/$html_log_file
$script_path/sync_codebase.sh $av1_code/aom >> $log_path/$html_log_file 2>&1
echo "</p>" >> $log_path/$html_log_file

echo "<p>" >> $log_path/$html_log_file
$script_path/aom_conf_build.sh $av1_code >> $log_path/$html_log_file 2>&1
echo "</p>" >> $log_path/$html_log_file

pdfps=`cat $log_path/$prev_log_file | grep e_ok | awk '{print $2}' | awk 'NR==1 {print $1}'`
petime=`cat $log_path/$prev_log_file | grep e_ok | awk '{print $1}' | awk 'NR==1 {print $1}'`
$script_path/aom_nightly_speed.sh $av1_code $pdfps $petime $html_log_file >> $log_path/$log_file 2>&1

# Send an email to coworkers
users=luoyi
host_name=`hostname`
sender=luoyi
cc_list="--cc=yunqingwang,vpx-eng"

$script_path/gen_html_footer.sh >> $log_path/$html_log_file

sendgmr --to=$users $cc_list --subject="AV1 Encoder Nightly Report" --from=$sender --reply_to=$sender --html_file=/usr/local/google/home/luoyi/Dev/log/$html_log_file --body_file=/usr/local/google/home/luoyi/Dev/log/$html_log_file
