#!/bin/sh
# Note:
# - Be very careful to type the report HTML file name (without path)
# - Use prodcertstatus --show_expiration_time to see when ticket expires.
#

report_file=$1

sendgmr --to=luoyi --cc=yunqingwang,vpx-eng --subject="Codec Daily Report" --from=luoyi --reply_to=luoyi --html_file=/usr/local/google/home/luoyi/Dev/log/$report_file --body_file=/usr/local/google/home/luoyi/Dev/log/$report_file
