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

@[525-529](Initialize MPI)
@[539-544](Calculate local chromosomes)
@[549-553](Each process send it's best to the master)
@[556-563](Master receives all bests...)
@[565-572](and decides which process had the best result)
@[574-582](The best process will receive TRUE and print its result)
@[589-602](The best process prints its result)

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

@[6-13](Declare variables)
@[546-552](Initialize PThreads)
@[553-558](Calculate local generations and starts each thread)
@[380-382](Change function call to suit pthreads)
@[439-443](Locks the variable to initialize the first run)
@[462-465](Locks when a improvement was detected)

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

+++

### Combination of both codes shown before

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

## Conclusion

<br>

| Method | Improved | Runtime |
| ------ | :------: | :-----: |
| MPI | True | -1.06s |
| PThreads | False | +0.30s |
| Both | True | -0.45s |

---

## Thank you!
#### https://github.com/arthurazs/uff-lpp
