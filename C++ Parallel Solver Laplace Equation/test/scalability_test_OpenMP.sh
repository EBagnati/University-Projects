#!/bin/bash


cd ../src/

echo "============================================================"
echo "Scalability Test"
echo "============================================================"

echo -e "\n1 PROCESS - 4 THREADS\n"
mpiexec -np 1 -x OMP_NUM_THREADS=4 ./main

echo -e "\n2 PROCESSES - 4 THREADS\n"
mpiexec -np 2 -x OMP_NUM_THREADS=4 ./main

echo -e "\n4 PROCESSES - 4 THREADS\n"
mpiexec -np 4 -x OMP_NUM_THREADS=4 ./main


