#!/bin/bash

echo "Compiling..."
c++ thread/hea.cpp -o bin/HEA_thread.o -lpthread -std=c++11 -O3
echo "done"
echo ""
