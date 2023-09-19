#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/code/sim_studies/benchmark_sim_chr

declare -a N=2000
echo $N

declare -a NPs=(12 20 2 3 4 5 8)
# declare -a NPs=(12)
npslength=${#NPs[@]}

declare -a sparseLevels=(1 2 3)
# declare -a sparseLevels=(1)
slslength=${#sparseLevels[@]}

declare -a alpha=($(seq 0.001 0.001 0.005) $(seq 0.01 0.01 0.1) 0.12 0.15 0.2 0.3 0.4)
# declare -a alpha=(0.05)
alplength=${#alpha[@]}

for (( i=0; i<${npslength}; i++ )); do
  for (( j=0; j<${slslength}; j++ )); do
    for (( k=0; k<${alplength}; k++ )); do
      echo ${NPs[$i]}
      echo ${sparseLevels[$j]}
      echo ${alpha[$k]}
      my_command="qsub -v N=${N},NP=${NPs[$i]},sparseLevel=${sparseLevels[$j]},alpha=${alpha[$k]} 3-0run_vmrseq.pbs"
      eval $my_command
    done
  done
done