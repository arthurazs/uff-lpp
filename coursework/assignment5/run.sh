#!/bin/bash

mpirun -n 4 bin/HEA.o -c input/cluster.vcl -w input/miniworkflow.dag
