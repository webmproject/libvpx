#!/bin/sh
self="$0"

usage() {
  cat <<EOF >&2
Usage: $self [option]

This script applies a whitespace transformation to the commit at HEAD. If no
options are given, then the modified files are left in the working tree.

Options:
  -n, --dry-run  Shows a diff of the changes to be made.
  --amend        Squashes the changes into the commit at HEAD
  --commit       Creates a new commit containing only the whitespace changes
EOF
  rm -f ${CLEAN_FILES}
  exit 1
}


log() {
  echo "${self##*/}: $@" >&2
}


vpx_style() {
  astyle --style=bsd --min-conditional-indent=0 --break-blocks \
         --pad-oper --pad-header --unpad-paren \
         --align-pointer=name \
         --indent-preprocessor --convert-tabs --indent-labels \
         --suffix=none --quiet "$@"
  sed -i 's/[[:space:]]\{1,\},/,/g' "$@"
}


apply() {
  patch -p1 < "$1"
}


commit() {
  LAST_CHANGEID=$(git show | awk '/Change-Id:/{print $2}')
  if [ -z "$LAST_CHANGEID" ]; then
    log "HEAD doesn't have a Change-Id, unable to generate a new commit"
    exit 1
  fi

  # Build a deterministic Change-Id from the parent's
  NEW_CHANGEID=${LAST_CHANGEID}-styled
  NEW_CHANGEID=I$(echo $NEW_CHANGEID | git hash-object --stdin)

  # Commit, preserving authorship from the parent commit.
  git commit -a -C HEAD > /dev/null
  git commit --amend -F- << EOF
Cosmetic: Fix whitespace in change ${LAST_CHANGEID:0:9}

Change-Id: ${NEW_CHANGEID}
EOF
}


amend() {
  git commit -a --amend -C HEAD
}


# Temporary files
ORIG_DIFF=orig.diff.$$
MODIFIED_DIFF=modified.diff.$$
FINAL_DIFF=final.diff.$$
CLEAN_FILES="${ORIG_DIFF} ${MODIFIED_DIFF} ${FINAL_DIFF}"

# Preconditions
[ $# -lt 2 ] || usage

if ! git diff --quiet HEAD; then
  log "Working tree is dirty, commit your changes first"
  exit 1
fi

# Need to be in the root
cd "$(git rev-parse --show-toplevel)"

# Collect the original diff
git show > "${ORIG_DIFF}"

# Apply the style guide on the modified files and collect its diff
for f in $(git diff HEAD^ --name-only | grep '\.[ch]$'); do
  case "$f" in
    third_party/*) continue;;
    nestegg/*) continue;;
  esac
  vpx_style "$f"
done
git diff --no-color --no-ext-diff > "${MODIFIED_DIFF}"

# Intersect the two diffs
$(dirname ${self})/intersect-diffs.py \
    "${ORIG_DIFF}" "${MODIFIED_DIFF}" > "${FINAL_DIFF}"
INTERSECT_RESULT=$?
git reset --hard >/dev/null

if [ $INTERSECT_RESULT -eq 0 ]; then
  # Handle options
  if [ -n "$1" ]; then
    case "$1" in
      -h|--help) usage;;
      -n|--dry-run) cat "${FINAL_DIFF}";;
      --commit) apply "${FINAL_DIFF}"; commit;;
      --amend) apply "${FINAL_DIFF}"; amend;;
      *) usage;;
    esac
  else
    apply "${FINAL_DIFF}"
    if ! git diff --quiet; then
      log "Formatting changes applied, verify and commit."
      log "See also: http://www.webmproject.org/code/contribute/conventions/"
      git diff --stat
    fi
  fi
fi

rm -f ${CLEAN_FILES}
