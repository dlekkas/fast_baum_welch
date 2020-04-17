#!/usr/bin/env bash

echo > out

output_dir=../output/
echo "Running on traces:"

for i in $(seq 0 20)
do
    obs="../input/observations_${i}.txt"
    s="../input/sample_${i}.txt"
    o="output_cpp_${i}.txt"
    ../bin/baum_welch $s $obs > $output_dir/$o
    echo "$i"
done
