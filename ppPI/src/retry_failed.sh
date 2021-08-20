#!/bin/bash

log_file=parallel_log.txt
failed_file=failed_replicas.txt

echo "=== FAILED RUNS ===" >> $log_file

i=0
exp_name=""
pdf_rep=""
ff_rep=""

for line in $(cat ${failed_file}); do

	if [[ $i -eq 0 ]]; then
		exp_name=$line
		i=1
	elif [[ $i -eq 1 ]]; then
		pdf_rep=$line
		i=2
	elif [[ $i -eq 2 ]]; then
		ff_rep=$line
		i=0

		echo $exp_name >> $log_file
		./run.sh $exp_name "${pdf_rep} ${ff_rep}"
	fi
done

echo "=== DONE ===" >> $log_file

echo "----- DONE! -----"
