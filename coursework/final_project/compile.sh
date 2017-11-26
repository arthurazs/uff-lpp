#!/bin/bash

echo ""

echo "Compiling Thread"
date
sh compile_thread.sh

# echo "Compiling OpenMP"
# date
# sh compile_openmp.sh

echo "Compiling MPI"
date
sh compile_mpi.sh

echo "Compiling MPI + Thread"
date
sh compile_mpi_thread.sh

# echo "Compiling MPI + OpenMP"
# date
# sh compile_mpi_openmp.sh

echo "Script finished"
date
