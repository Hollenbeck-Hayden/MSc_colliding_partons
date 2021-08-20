#!/bin/bash

./reweight "$@"
./chi2 "$@"
python plot.py "$@"
