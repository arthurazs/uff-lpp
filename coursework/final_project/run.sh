#!/bin/bash

echo ""

echo "Starting Thread"
date
sh run_thread.sh

# echo "Starting OpenMP"
# date
# sh run_openmp.sh

echo "Starting MPI"
date
sh run_mpi.sh

echo "Starting MPI + Thread"
date
sh run_mpi_thread.sh

# echo "Starting MPI + OpenMP"
# date
# sh run_mpi_openmp.sh

echo "Script finished"
date
