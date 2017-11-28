#!/bin/bash

echo "Compiling..."
mpic++ mpi_openmp/hea.cpp -o bin/HEA_mpi_openmp.o -fopenmp -std=c++11 -O3
echo "done"
echo ""
