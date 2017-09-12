# Parallel Programming Lab (UFF)

This repository contains all the assignments implemented for the
Parallel Programming Lab at the Universidade Federal Fluminense (Masters
Degree).

**Contents**

- [openmpi](#openmpi)
    - [Installation](#installation)
    - [Implementation](#implementation)
- [Coursework](#coursework)
    - [Matrix Multiplication (1st Assignment)](#matrix-multiplication)

## openmpi

### Installation

    sudo apt install openmpi-bin openssh-client openssh-server libopenmpi-dev

### Implementation

Check [teste.c](teste.c) for a simple example:

    $ mpicc teste.c -o teste
    $ mpirun -n 4 teste
    a Process 0 on inspiron out of 4
    a Process 1 on inspiron out of 4
    a Process 2 on inspiron out of 4
    a Process 3 on inspiron out of 4

## Coursework

- [Matrix Multiplication (1st Assignment)](#matrix-multiplication)

### Matrix Multiplication

**TODO**
See [matrix_multiplication.c](coursework/assignment1/matrix_multiplication.c).
