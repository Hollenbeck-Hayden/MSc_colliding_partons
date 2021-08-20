#!/bin/bash

log_file=parallel_log.txt
exp_name=$1


echo "=== ${exp_name} ===" >> $log_file


# Build the replica file
replica_file=../theory/$exp_name/pdf_replicas.txt
touch "${replica_file}"
echo 0 0 > ${replica_file}
for i in {1..200}; do
	pdf_replica=$[$RANDOM % 100 + 1]
	echo "${i} ${pdf_replica}" >> "${replica_file}"
done

# Run in parallel
time parallel -j+0 --eta --link ./run.sh ::: $exp_name :::: ${replica_file}

echo "--- DONE ---" >> $log_file
echo "" >> $log_file
