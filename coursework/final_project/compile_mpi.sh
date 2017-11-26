#!/bin/bash

echo "Compiling..."
mpic++ mpi/hea.cpp -o bin/HEA_mpi.o -std=c++11
echo "done"
echo ""
