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
# Updated: Oct 18, 2022
# Author: M. Chatzikos
# Comment: Added support for compilation sans git (e.g., release tarballs).
#
# Updated: Jan 26, 2024
# Author: M. Chatzikos
# Comment: Implement Christophe Morisset's request to emit bare numbers for
#	   release tags.  Simplify logic along the way.
#
# Updated: Dec 27, 2024
# Author: P.A.M. van Hoof
# Comment: Handle working copies in a detached state correctly (see PR #497)
#

is_repo=`git rev-parse --is-inside-work-tree 2>&1 | grep true`
if [[ $is_repo != 'true' ]]; then
	#
	# This branch is for tarball releases.
	# The release number is set internally (version.cpp).
	#
	echo
	exit 1
fi

tag=`git describe --tags --abbrev=0 2> /dev/null`
if ! [ -z "$tag" ]; then
	#
	# Cloudy version number, sans the initial 'c'.
	# This branch is for official releases.
	#
	tag=`echo $tag | sed -E 's/^(c|C)//'`
	echo $tag
else
	#
	# This branch is for development.
	#
	sha1=`git log --oneline | head -n 1 | awk '{print $1}'`
	branch=`git branch | grep '^\*' | sed 's/(HEAD//' | awk '{ print $2 }'`
	branch=`echo $branch | sed 's/(no//'`
	[[ -z "`git status -s -uno`" ]] && modified="" || modified="-modified"
	
	if [ -z "$branch" ];
	then
		echo $sha1$modified
	else
		echo $branch-$sha1$modified
	fi
fi
