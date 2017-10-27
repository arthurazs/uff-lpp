#!/bin/bash

mpirun -n 4 bin/HEA.o -c input/cluster.vcl -w input/miniworkflow.dag
# ./bin/HEA -c input/cluster.vcl -w input/miniworkflow.dag -v
