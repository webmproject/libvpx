#!/usr/bin/python

import subprocess
import sys

def RunCommand(command):
  run = subprocess.Popen(command, shell=True)
  output = run.communicate()
  if run.returncode:
    print "Non-zero return code: " + str(run.returncode) + " => exiting!"
    sys.exit(1)

def list_of_experiments():
  experiments = []
  configure_file = open("configure")
  list_start = False
  for line in configure_file.read().split("\n"):
    if line == 'EXPERIMENT_LIST="':
      list_start = True
    elif line == '"':
      list_start = False
    elif list_start:
      currently_broken = ["csm"]
      experiment = line[4:]
      if experiment not in currently_broken:
        experiments.append(experiment)
  return experiments

def main():
  base_command = "./configure --enable-internal-stats"
  test_build(base_command)
  for experiment_name in list_of_experiments():
    test_build("%s --enable-experimental --enable-%s" % (base_command,
      experiment_name))

def test_build(configure_command):
  print "\033[34m\033[47mTesting %s\033[0m" % (configure_command)
  RunCommand(configure_command)
  RunCommand("make clean")
  RunCommand("make")

if __name__ == "__main__":
  main()
