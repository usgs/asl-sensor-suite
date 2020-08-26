#!/bin/bash
# Move to scripts directory
cd "$(dirname "$0")"

rm *.jar
git fetch --all
git checkout -f origin/staging

gradle build
