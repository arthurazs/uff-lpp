# Parallel Programming Lab (UFF)

This repository contains all the assignments implemented for the
Parallel Programming Lab at the Universidade Federal Fluminense
(Master's Degree).

**Contents**

- [openmpi](#openmpi)
    - [Installation](#installation)
    - [Implementation](#implementation)
- [Coursework](#coursework)
    - [Matrix Multiplication (1st Assignment)](#matrix-multiplication)
    - [Definite Integral (2nd Assignment)](#definite-integral)
    - [Final Project (HEA)](#final-project)
        - [PThreads](#pthreads)
        - [OpenMP](#openmp)
        - [MPI](#mpi)
            - [MPI + PThreads](#mpi-pthreads)
            - [MPI + OpenMP](#mpi-openmp)

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
- [Definite Integral (2nd Assignment)](#definite-integral)
- [Final Project](#final-project)

### Matrix Multiplication

Check the [slide](https://gitpitch.com/arthurazs/uff-lpp/master?p=coursework/assignment1).

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


### Definite Integral

Check the [slide](https://gitpitch.com/arthurazs/uff-lpp/master?p=coursework/assignment2).

Numerical Integration of a Definite Integral using MPI.
Integral (1 to 6) of **1 / x**:

![Graph](coursework/assignment2/graph.png)

![Math](coursework/assignment2/math.png)

The graph was generated using
[symbolab.com](https://www.symbolab.com/solver/definite-integral-calculator/%5Cint_%7B1%7D%5E%7B6%7D%20%5Cfrac%7B1%7D%7Bx%7D%20dx),
check it out!

Running [integral_mpi.c](coursework/assignment2/integral_mpi.c):

    $ cd uff-lpp/coursework/assignment2
    $ mpicc integral_mpi.c -o integral
    $ mpirun -n 4 integral

    P(3): Estimated 0.233402 for a(4.75) to b(6.00)
    P(0): Estimated 0.812278 for a(1.00) to b(2.25)
    P(2): Estimated 0.305017 for a(3.50) to b(4.75)
    P(1): Estimated 0.441064 for a(2.25) to b(3.50)

    P(0): MPI_Reduce
    Integral (1 to 6) of 1/x

    Expected  result  = 1.79175946923
    Estimated result  = 1.79176079291
                        -------------
    Difference        = 0.00000132368


### Final Project

[HEA](https://github.com/luanteylo/HEA) is a Hybrid Evolutionary
Algorithm by [Luan Teylo](https://github.com/luanteylo) developed to
solve the problem of task scheduling and data allocations for scientific
workflows.
HEA is a linear code with many heavy computational algorithms, meaning
it may be slow to solve big problems. This Final Project aims to study
HEA's linear code, implement several parallelism techniques and compare
its performances

You can check the presentation's slide SOON.

HEA's [linear code](https://github.com/luanteylo/HEA/blob/master/source/hea.cpp) performance was compared with PThreads, OpenMP, MPI and a combination of MPI + PThreads and MPI + OpenMP.

- [PThreads](#pthreads)
- [OpenMP](#openmp)
- [MPI](#mpi)
    - [MPI + PThreads](#mpi-openmp)
    - [MPI + OpenMP](#mpi-pthreads)

#### PThreads

    TODO

#### OpenMP

    TODO

#### MPI

    TODO

#### MPI + PThreads

    TODO

#### MPI + OpenMP

    TODO

## Authors

- Matrix Multiplication & Definite Integral
    - [Arthur Zopellaro](https://github.com/arthurazs)
- Final Project
    - [Arthur Zopellaro](https://github.com/arthurazs)
    - Marcos Melos
