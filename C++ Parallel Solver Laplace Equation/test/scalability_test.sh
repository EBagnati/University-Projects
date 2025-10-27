#!/bin/bash


cd ../src/

echo "============================================================"
echo "Scalability Test"
echo "============================================================"

echo -e "\n1 PROCESS\n"
mpiexec -np 1 ./main

echo -e "\n2 PROCESSES\n"
mpiexec -np 2 ./main

echo -e "\n4 PROCESSES\n"
mpiexec -np 4 ./main


