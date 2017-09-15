## Matrix  Multiplication
### openmpi

Arthur Zopellaro

---

## Code in C

+++

```c
for (int row = 0; row < A_ROWS; row++) {
    for (int col = 0; col < B_COLS; col++) {
        int sum = 0;
        for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
            sum += matrix_a[row][ctrl] * matrix_b[ctrl][col];
        matrix_x[row][col] = sum;
    }
}

```
Linear Programming
---

## Code in mpi

+++

The matrices A[3][2] and B[2][3] will be the following:

```c
int matrix_a [A_ROWS][A_COLS] = {
    {9, 0},
    {5, 6},
    {1, 2}
};

int matrix_b [B_ROWS][B_COLS] = {
    {2, 4, 3},
    {7, 8, 9}
};
```

$ mpirun **-n 4**

+++

```c
for (int row = 0; row < A_ROWS; row++) {
    for (int col = 0; col < B_COLS; col++) {
        int sum = 0;
        for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
            sum += matrix_a[row][ctrl] * matrix_b[ctrl][col];
        matrix_x[row][col] = sum;
    }
}

```

@[1](Remove the 1st `for`)
@[2-6](Run the rest in process 0, 1, 2)

+++

```c
//for (int row = 0; row < A_ROWS; row++) {
    for (int col = 0; col < B_COLS; col++) {
        int sum = 0;
        for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
            sum += matrix_a[row][ctrl] * matrix_b[ctrl][col];
        matrix_x[row][col] = sum;
    }
//}
```

@[1, 8](Removed the 1st `for`)
@[5-6](But now we don't have `row`, so...)

+++

```c
if(rank < 3) {
    int line[X_COLS];
    for (int col = 0; col < B_COLS; col++) {
        int sum = 0;
        for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
            sum += matrix_a[row][ctrl] * matrix_b[ctrl][col];
        // matrix_x[row][col] = sum;
        line[]
    }
}
```

@[1, 10](Since we are using `openmpi`, we will put our code to be used for ranks from 0 to 2)
@[5-6](But now we don't have `row`, so...)
