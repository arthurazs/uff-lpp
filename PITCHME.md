## Parallel Programming
### HEA

Arthur Zopellaro

---

## Linear Code
#### Running

+++

```

$ c++ hea.cpp -o HEA.o -std=c++11
$ ./HEA.o -w miniworkflow.dag -c cluster.vcl

Best fitness: 13.9833(min) Runtime: 1.3814(sec)

```

---

## MPI
#### Code

+++?code=coursework/assignment3/source/hea.cpp&lang=cpp

@[28-37](EXPLAIN)

+++

### Running

+++

```
$ cd uff-lpp/coursework/assignment3
$ mpic++ source/hea.cpp -std=c++11 -o bin/HEA.o
$ mpirun -n 4 bin/HEA.o -c input/cluster.vcl -w input/miniworkflow.dag

Best fitness found by rank: 0
Best fitness: 13.9833(min) Runtime: 0.319574(sec)

```
@[6](Runtime has decreased by 1.06 seconds: 76.86% Improvement)

---

## PThreads
#### Code

+++?code=coursework/assignment4/source/hea.cpp&lang=cpp

@[28-37](EXPLAIN)

+++

### Running

+++

```
$ cd uff-lpp/coursework/assignment4
$ c++ source/hea.cpp -o bin/HEA.o -lpthread -std=c++11
$ ./bin/HEA.o -c input/cluster.vcl -w input/miniworkflow.dag

Best fitness: 13.9833(min) Runtime: 1.68889(sec)

```

@[5](Runtime has increased by 0.30 seconds: 22.26% Decline)

---

## MPI + PThreads
#### Code

+++?code=coursework/assignment5/source/hea.cpp&lang=cpp

@[28-37](EXPLAIN)

+++

### Running

+++

```
$ cd uff-lpp/coursework/assignment5
$ c++ source/hea.cpp -o bin/HEA.o -lpthread -std=c++11
$ ./bin/HEA.o -c input/cluster.vcl -w input/miniworkflow.dag

Best fitness found by rank: 0
Best fitness: 13.9833(min) Runtime: 0.931743(sec)

```

@[6](Runtime has decreased by 0.45 seconds: 32.55% Improvement)

---

## Thank you!
#### All codes at https://github.com/arthurazs/uff-lpp
