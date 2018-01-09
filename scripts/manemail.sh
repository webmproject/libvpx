#!/bin/sh
# Note:
# - With any path, be very careful to type the report HTML file name.
# - Use prodcertstatus --show_expiration_time to see when ticket expires.
#

report_html=$1
report_path=/usr/local/google/home/luoyi/Dev/log
report_send=$report_path/$report_html

sendgmr --to=luoyi --cc=yunqingwang,vpx-eng --from=luoyi --reply_to=luoyi --subject="AV1 Nightly Speed Report" --html_file=$report_send --body_file=$report_send
