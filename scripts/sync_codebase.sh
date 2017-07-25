#!/bin/sh

code_dir=$1

cd $code_dir
git checkout -q master
git pull -q
git log -1 --oneline
