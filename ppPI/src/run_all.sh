#!/bin/bash

log_file=parallel_log.txt
output_file=parallel_output.txt
exp_file=../data/successful_experiments.txt

rm -rf ${log_file} ${output_file}
touch ${log_file} ${output_file}
for exp in $(cat ${exp_file}); do
	parallel_run.sh $exp >> ${output_file}
done
