#!/bin/bash

clear
echo "Compiling..."
/opt/openmpi-2.1.1/bin/mpic++ mpi_thread/hea.cpp -o bin/HEA_mpi_thread.o -std=c++11 -lpthread
echo "done"
