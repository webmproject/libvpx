# Copyright (c) 2014 The WebM project authors. All Rights Reserved.
#
# Use of this source code is governed by a BSD-style license
# that can be found in the LICENSE file in the root of the source
# tree. An additional intellectual property rights grant can be found
# in the file PATENTS.  All contributing project authors may
# be found in the AUTHORS file in the root of the source tree.

"""Standalone script which parses a gtest log for json.

Json is returned returns as an array.  This script is used by the libvpx
waterfall to gather json results mixed in with gtest logs.  This is
dubious software engineering.
"""

import json
import re
import sys


def main():
  blob = sys.stdin.read()
  json_string = '[' + ','.join('{' + x + '}' for x in
                               re.findall(r'{([^}]*.?)}', blob)) + ']'
  print json.dumps(json.loads(json_string), indent=4, sort_keys=True)

if __name__ == '__main__':
  sys.exit(main())
