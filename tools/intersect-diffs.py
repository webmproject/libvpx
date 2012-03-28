#!/usr/bin/env python
##  Copyright (c) 2012 The WebM project authors. All Rights Reserved.
##
##  Use of this source code is governed by a BSD-style license
##  that can be found in the LICENSE file in the root of the source
##  tree. An additional intellectual property rights grant can be found
##  in the file PATENTS.  All contributing project authors may
##  be found in the AUTHORS file in the root of the source tree.
##
"""Calculates the "intersection" of two unified diffs.

Given two diffs, A and B, it finds all hunks in B that had non-context lines
in A and prints them to stdout. This is useful to determine the hunks in B that
are relevant to A. The resulting file can be applied with patch(1) on top of A.
"""

__author__ = "jkoleszar@google.com"

import re
import sys


class DiffLines(object):
    """A container for one half of a diff."""

    def __init__(self, filename, offset, length):
        self.filename = filename
        self.offset = offset
        self.length = length
        self.lines = []
        self.delta_line_nums = []

    def Append(self, line):
        l = len(self.lines)
        if line[0] != " ":
            self.delta_line_nums.append(self.offset + l)
        self.lines.append(line[1:])
        assert l+1 <= self.length

    def Complete(self):
        return len(self.lines) == self.length

    def __contains__(self, item):
        return item >= self.offset and item <= self.offset + self.length - 1


class DiffHunk(object):
    """A container for one diff hunk, consisting of two DiffLines."""

    def __init__(self, header, file_a, file_b, start_a, len_a, start_b, len_b):
        self.header = header
        self.left = DiffLines(file_a, start_a, len_a)
        self.right = DiffLines(file_b, start_b, len_b)
        self.lines = []

    def Append(self, line):
        """Adds a line to the DiffHunk and its DiffLines children."""
        if line[0] == "-":
            self.left.Append(line)
        elif line[0] == "+":
            self.right.Append(line)
        elif line[0] == " ":
            self.left.Append(line)
            self.right.Append(line)
        else:
            assert False, ("Unrecognized character at start of diff line "
                           "%r" % line[0])
        self.lines.append(line)

    def Complete(self):
        return self.left.Complete() and self.right.Complete()

    def __repr__(self):
        return "DiffHunk(%s, %s, len %d)" % (
            self.left.filename, self.right.filename,
            max(self.left.length, self.right.length))


def ParseDiffHunks(stream):
    """Walk a file-like object, yielding DiffHunks as they're parsed."""

    file_regex = re.compile(r"(\+\+\+|---) (\S+)")
    range_regex = re.compile(r"@@ -(\d+)(,(\d+))? \+(\d+)(,(\d+))?")
    hunk = None
    while True:
        line = stream.readline()
        if not line:
            break

        if hunk is None:
            # Parse file names
            diff_file = file_regex.match(line)
            if diff_file:
              if line.startswith("---"):
                  a_line = line
                  a = diff_file.group(2)
                  continue
              if line.startswith("+++"):
                  b_line = line
                  b = diff_file.group(2)
                  continue

            # Parse offset/lengths
            diffrange = range_regex.match(line)
            if diffrange:
                if diffrange.group(2):
                    start_a = int(diffrange.group(1))
                    len_a = int(diffrange.group(3))
                else:
                    start_a = 1
                    len_a = int(diffrange.group(1))

                if diffrange.group(5):
                    start_b = int(diffrange.group(4))
                    len_b = int(diffrange.group(6))
                else:
                    start_b = 1
                    len_b = int(diffrange.group(4))

                header = [a_line, b_line, line]
                hunk = DiffHunk(header, a, b, start_a, len_a, start_b, len_b)
        else:
            # Add the current line to the hunk
            hunk.Append(line)

            # See if the whole hunk has been parsed. If so, yield it and prepare
            # for the next hunk.
            if hunk.Complete():
                yield hunk
                hunk = None

    # Partial hunks are a parse error
    assert hunk is None


def FormatDiffHunks(hunks):
    """Re-serialize a list of DiffHunks."""
    r = []
    last_header = None
    for hunk in hunks:
        this_header = hunk.header[0:2]
        if last_header != this_header:
            r.extend(hunk.header)
            last_header = this_header
        else:
            r.extend(hunk.header[2])
        r.extend(hunk.lines)
        r.append("\n")
    return "".join(r)


def ZipHunks(rhs_hunks, lhs_hunks):
    """Join two hunk lists on filename."""
    for rhs_hunk in rhs_hunks:
        rhs_file = rhs_hunk.right.filename.split("/")[1:]

        for lhs_hunk in lhs_hunks:
            lhs_file = lhs_hunk.left.filename.split("/")[1:]
            if lhs_file != rhs_file:
                continue
            yield (rhs_hunk, lhs_hunk)


def main():
    old_hunks = [x for x in ParseDiffHunks(open(sys.argv[1], "r"))]
    new_hunks = [x for x in ParseDiffHunks(open(sys.argv[2], "r"))]
    out_hunks = []

    # Join the right hand side of the older diff with the left hand side of the
    # newer diff.
    for old_hunk, new_hunk in ZipHunks(old_hunks, new_hunks):
        if new_hunk in out_hunks:
            continue
        old_lines = old_hunk.right
        new_lines = new_hunk.left

        # Determine if this hunk overlaps any non-context line from the other
        for i in old_lines.delta_line_nums:
            if i in new_lines:
                out_hunks.append(new_hunk)
                break

    if out_hunks:
        print FormatDiffHunks(out_hunks)
        sys.exit(1)

if __name__ == "__main__":
    main()
