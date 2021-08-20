#!/bin/bash

path=$1

for f in $path/*.png; do feh $f; done
