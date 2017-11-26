#!/bin/bash

clear
echo "Compiling..."
/opt/openmpi-2.1.1/bin/mpic++ mpi/hea.cpp -o bin/HEA_mpi.o -std=c++11
echo "done"
