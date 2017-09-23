## Definite Integral
### openmpi

Arthur Zopellaro

---

## Integral (1 to 6) of 1/x
[symbolab.com](https://www.symbolab.com/solver/definite-integral-calculator/)

+++

`$\int_{1}^{6}\frac{1}{x}dx$`

`$ln(6)\approx1.79175946923...$`

+++

![Graph](coursework/assignment2/graph.png)

---

## Code in openmpi
#### Using MPI_Reduce

+++?code=coursework/assignment2/integral_mpi.c&lang=c

@[29-36](Initializing variables)
@[37-49](Allocating remaining processes)
@[51-52](`Calls calculate()`)
@[5-19](Integrationg of f(x`)` = 1/x)
@[53-61](Prints a feedback and SUM the integration using MPI_Reduce)
@[63-76](Prints the results)

---

## Running the code

+++

```
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

```

@[1-3](Compiling and running the code)
@[5-8](Feedback [MPI_Recv])
@[10-16](Comparison between expected and obtained results)

---

## Thank you!
#### Check the code at https://github.com/arthurazs/uff-lpp
