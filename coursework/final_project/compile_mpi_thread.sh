#!/bin/bash

echo "Compiling..."
mpic++ mpi_thread/hea.cpp -o bin/HEA_mpi_thread.o -std=c++11 -lpthread -O3
echo "done"
echo ""
