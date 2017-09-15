## Matrix  Multiplication
### openmpi

Arthur Zopellaro

---

## Code in C
#### Linear Programming

+++?code=coursework/assignment1/matrix_c.c&lang=c

@[28-37](Matrices A[3][2] and B[2][3])
@[41-49](Matrix Multiplication)

---

## Code in openmpi
#### Parallel Programming
$ mpirun **-n 4**

+++?code=coursework/assignment1/matrix_mpi.c&lang=c

@[39-48](Same matrices)
@[40-42, 46-47](We'll use each process to multiply one line with each column)
@[40](Process 0)
@[41](Process 1)
@[42](Process 2)
@[40-42, 46-47](Which leaves us with one extra process)
@[50](Which will be used to join the results of each process into the Matrix X)

---

## Adapting the C code
#### for openmpi

+++

```c
//for (int row = 0; row < A_ROWS; row++) {
if(rank < 3) {
    int line[X_COLS];
    for (int col = 0; col < B_COLS; col++) {
        int sum = 0;
        for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
            //sum += matrix_a[row][ctrl] * matrix_b[ctrl][col];
            sum += matrix_a[rank][ctrl] * matrix_b[ctrl][col];
        //matrix_x[row][col] = sum;
        line[col] = sum
    }
}
```

@[1-2](Removed the 1st `for`)
@[3, 9, 10](Use array instead of Matrix X)
@[4-6](Remains the same)
@[7-8](Use rank instead of row)

---

## Final code in openmpi
#### Sending and receiving messages

+++?code=coursework/assignment1/matrix_mpi.c&lang=c

@[52-63](Processes 0, 1 and 2)
@[54](Prints a feedback)
@[55-60](The multiplication)
@[61](Sends a message with the result to the master (Process 3))
@[64-77](Process 3)
@[65](Prints a feedback)
@[66-72](Calls MPI_Recv for each origin (Processes 0, 1 and 2))
@[68](Although MPI_Recv will block the reception buffer)
@[69-70](Copy the values from the array of the origin to the Matrix X)
@[71](Prints a feedback)
@[73-75](After receiving all messages, it prints Matrices A, B and X)


---

## Running the code

+++

```
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
```

@[1-3](Compiling and running the code)
@[5-11](Feedback (MPI_Recv))
@[13-25](Matrices)

---

## Future works
#### Ways to improve this code

+++?code=coursework/assignment1/matrix_mpi.c&lang=c

@[37](Stop limiting the code to specifically X_ROWS + 1 processes)
@[79-85](Stop limiting the code to specifically X_ROWS + 1 processes)
@[66-72](Implement MPI_Irecv)

---

## Thank you!
#### All the codes are at https://github.com/arthurazs/uff-lpp
