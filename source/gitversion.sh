#!/usr/bin/env bash
#
# Label the version of Cloudy to be shown on its output, using git metadata.
# The branch name will be used, along with the most recent tag, if found, or the
# SHA1 string of the most recent commit.  If the code has been modified, a
# string indicating so is also appended.  The label is dash-delimited.
#
#
# Created: Dec 05, 2020
# Author: M. Chatzikos
#
# Updated: Dec 11, 2020
# Author: M. Chatzikos
# Comment: Added support for git tags.
#

branch=`git branch | grep '^\*' | awk '{ print $2 }'`
tag=`git describe --tags --abbrev=0 2> /dev/null`
if [ -z "$tag" ]; then
	tag=`git log --oneline | head -n 1 | awk '{print $1}'`
fi
[[ -z "`git status -s -uno`" ]] && modified="" || modified="-modified"
echo $branch-$tag$modified
