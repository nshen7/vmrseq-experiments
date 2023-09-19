#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/code/case_studies/luo2017mice_full

declare -a chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)
chrslength=${#chrs[@]}

for (( i=0; i<${chrslength}; i++ )); do
    echo $chrs[$i]
    my_command="qsub -v chr=${chrs[$i]} 2-3format_data_smallwood.pbs"
    eval $my_command
done