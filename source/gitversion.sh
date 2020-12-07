#!/usr/bin/env bash
#
# Set the version of the branch to be used on Cloudy's output.
#
# Created: Dec 05, 2020
# Author: M. Chatzikos
#

branch=`git branch | grep '^\*' | gawk  -e '{ print $2 }'`
gitrev=`git log --oneline | head -n 1 | gawk -e '{print $1}'`
[[ -z "`git status -s -uno`" ]] && modified="" || modified="-modified"
echo $branch-$gitrev$modified
