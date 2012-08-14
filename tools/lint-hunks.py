#!/usr/bin/python
##  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##
"""Performs style checking on each diff hunk."""
import getopt
import os
import StringIO
import subprocess
import sys

import diff


SHORT_OPTIONS = "h"
LONG_OPTIONS = ["help"]

TOPLEVEL_CMD = ["git", "rev-parse", "--show-toplevel"]
DIFF_CMD = ["git", "diff-index", "-u", "--cached", "HEAD", "--"]
SHOW_CMD = ["git", "show"]
CPPLINT_FILTERS = ["-readability/casting", "-runtime/int"]


class Usage(Exception):
    pass


class SubprocessException(Exception):
    def __init__(self, args):
        msg = "Failed to execute '%s'"%(" ".join(args))
        super(SubprocessException, self).__init__(msg)


class Subprocess(subprocess.Popen):
    """Adds the notion of an expected returncode to Popen."""

    def __init__(self, args, expected_returncode=0, **kwargs):
        self._args = args
        self._expected_returncode = expected_returncode
        super(Subprocess, self).__init__(args, **kwargs)

    def communicate(self, *args, **kwargs):
        result = super(Subprocess, self).communicate(*args, **kwargs)
        if self._expected_returncode is not None:
            try:
                ok = self.returncode in self._expected_returncode
            except TypeError:
                ok = self.returncode == self._expected_returncode
            if not ok:
                raise SubprocessException(self._args)
        return result


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, _ = getopt.getopt(argv[1:], SHORT_OPTIONS, LONG_OPTIONS)
        except getopt.error, msg:
            raise Usage(msg)

        # process options
        for o, _ in opts:
            if o in ("-h", "--help"):
                print __doc__
                sys.exit(0)

        # Find the fully qualified path to the root of the tree
        tl = Subprocess(TOPLEVEL_CMD, stdout=subprocess.PIPE)
        tl = tl.communicate()[0].strip()

        # Build the command line to execute cpplint
        cpplint_cmd = [os.path.join(tl, "tools", "cpplint.py"),
                       "--filter=" + ",".join(CPPLINT_FILTERS),
                       "-"]

        # Get a list of all affected lines
        file_affected_line_map = {}
        p = Subprocess(DIFF_CMD, stdout=subprocess.PIPE)
        stdout = p.communicate()[0]
        for hunk in diff.ParseDiffHunks(StringIO.StringIO(stdout)):
            filename = hunk.right.filename[2:]
            if filename not in file_affected_line_map:
                file_affected_line_map[filename] = set()
            file_affected_line_map[filename].update(hunk.right.delta_line_nums)

        # Run each affected file through cpplint
        lint_failed = False
        for filename, affected_lines in file_affected_line_map.iteritems():
            if filename.split(".")[-1] not in ("c", "h", "cc"):
                continue
            show_cmd = SHOW_CMD + [":" + filename]
            show = Subprocess(show_cmd, stdout=subprocess.PIPE)
            lint = Subprocess(cpplint_cmd, expected_returncode=(0, 1),
                              stdin=show.stdout, stderr=subprocess.PIPE)
            lint_out = lint.communicate()[1]
            if lint.returncode == 1:
                lint_failed = True

            for line in lint_out.split("\n"):
                fields = line.split(":")
                if fields[0] != "-":
                    continue
                warning_line_num = int(fields[1])
                if warning_line_num in affected_lines:
                    print "%s:%d:%s"%(filename, warning_line_num,
                                      ":".join(fields[2:]))

        # Propagate lint's exit status
        if lint_failed:
            return 1

    except Usage, err:
        print >>sys.stderr, err
        print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())
