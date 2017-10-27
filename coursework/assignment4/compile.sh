#!/bin/bash

clear
echo "Compiling..."
c++ source/hea.cpp -o bin/HEA.o -lpthread -std=c++11
echo "done"
