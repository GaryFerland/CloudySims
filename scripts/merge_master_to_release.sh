#!/bin/bash

set -e

echo "Pulling latest from origin..."
git fetch origin

echo "Checking out master branch..."
git checkout master

echo "🧹 Resetting master to match origin/master (fresh copy)..."
git reset --hard origin/master

echo "Updating copyright years..."
CURRENT_YEAR=$(date +%Y)
find ./ -type f -exec sed -i -e "s/1978-2023/1978-${CURRENT_YEAR}/g" {} \;

echo "Commit any files that have changed."
read -p "Press ENTER when done committing and pushing to master..."

echo "Checking out release branch..."
git switch release

echo "Pulling latest release branch from origin..."
git pull origin release

echo "Merging master into release (using 'theirs' strategy to favor master)..."
git merge -X theirs master -m "Merge master into release (resolved in favor of master)"

echo "Pushing updated release branch..."
git push origin release

echo "Merge process complete!"