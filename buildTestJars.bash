#!/bin/bash
git fetch --all
git checkout origin/devel
sed -i -E "s/(version\s*=\s*'[0-9]\.[0-9]\.[0-9])(')/\1_b$(date +%y%j)_$(git rev-parse --short HEAD)\2/" build.gradle
gradle copyJar
gradle copyServerJar
git checkout build.gradle
