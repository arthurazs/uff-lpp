#!/bin/bash

echo "Compiling..."
c++ openmp/hea.cpp -o bin/HEA_openmp.o -fopenmp -std=c++11 -O3
echo "done"
echo ""
