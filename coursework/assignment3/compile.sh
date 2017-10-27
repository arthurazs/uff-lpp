#!/bin/bash

clear
echo "Compiling..."
mpic++ source/hea.cpp -std=c++11 -o bin/HEA.o
# c++ source/hea.cpp -std=c++11 -o bin/HEA.o
echo "done"
