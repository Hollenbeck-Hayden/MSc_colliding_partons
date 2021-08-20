#!/bin/bash

name=$1
cutoff=$2
expfile=$3
nreplicas=200

path=${name}_${cutoff}

rm -rf $path
mkdir -p $path

echo -e "$nreplicas\n$cutoff\n" > $path/experiments.txt
for exp in $(cat $expfile); do echo -e "$exp\t1_3_3_1.0_1" >> $path/experiments.txt; done


plotexps=( "ALICE_2760a_PI0" "ALICE_2760b_PI0" "ALICE_7000_PI0" "ALICE_8000_PI0" "PHENIX_200b_PIsum" "PHENIX_200_PI0" "PHENIX_510_PI0" "STAR_200a_PI0" "STAR_200b_PI0" )

echo -e "$nreplicas\n$cutoff\n" > $path/plot_experiments.txt
for exp in ${plotexps[*]}; do echo -e "$exp\t1_3_3_1.0_1" >> $path/plot_experiments.txt; done
