#!/bin/bash
# Move to scripts directory
cd "$( dirname "${BASH_SOURCE[0]}" )"

rm *.jar
git fetch --all
git checkout -f origin/staging

# Set version number to #.#.#_bYYJJJ_CommitHash
sed -i -E "s/(version\s*=\s*'[0-9]*\.[0-9]*\.[0-9]*)(')/\1_b$(date +%y%j)_$(git rev-parse --short HEAD)\2/" build.gradle
gradle copyJar
gradle copyServerJar
# Undo our change
git checkout build.gradle