# Parallel Programming Lab (UFF)

This repository contains all the assignments implemented for the
Parallel Programming Lab at the Universidade Federal Fluminense (Master's
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

    $ cd uff-lpp
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

Check the [slide](https://gitpitch.com/arthurazs/uff-lpp).

See the source code [matrix_mpi.c](coursework/assignment1/matrix_mpi.c):

    $ cd uff-lpp/coursework/assignment1
    $ mpicc matrix_mpi.c -o matrix
    $ mpirun -n 4 matrix

    P(2): Multiplying row 2 with column 2
    P(1): Multiplying row 1 with column 1
    P(3): Waiting multiplication results
    P(3): Received message from P(0)
    P(3): Received message from P(1)
    P(0): Multiplying row 0 with column 0
    P(3): Received message from P(2)

    Matrix A [3][2]
    9 0
    5 6
    1 2

    Matrix B [2][3]
    2 4 3
    7 8 9

    Matrix X [3][3]
    18 36 27
    52 68 69
    16 20 21

### ?
TODO
