// mat_x num of rows equals to
// mat_a num of rows
#define A_ROWS 3
#define X_ROWS 3
// mat_x num of cols equals to
// mat_b num of cols
#define B_COLS 3
#define X_COLS 3
// mat_a num of cols should be equals to
// mat_b num of rows
#define A_COLS 2
#define B_ROWS 2

#include <stdio.h>
// #include <mpi.h>

void print_matrix(char* name, int rows, int cols, int matrix[rows][cols]) {

    printf("\n%s [%d][%d]\n", name, rows, cols);
    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++)
            printf("%d ", matrix[row][col]);
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

    int matrix_a [A_ROWS][A_COLS] = {
        {9, 0},
        {5, 6},
        {1, 2}
    };

    int matrix_b [B_ROWS][B_COLS] = {
        {2, 4, 3},
        {7, 8, 9}
    };

    int matrix_x [X_ROWS][X_COLS];

    // multipy matrices a and b
    for (int row = 0; row < A_ROWS; row++) {
        for (int col = 0; col < B_COLS; col++) {
            int sum = 0;
            for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
                sum = sum + matrix_a[row][ctrl] * matrix_b[ctrl][col];
            matrix_x[row][col] = sum;
        }
    }

    print_matrix("Matrix A", A_ROWS, A_COLS, matrix_a);
    print_matrix("Matrix B", B_ROWS, B_COLS, matrix_b);
    print_matrix("Matrix X", X_ROWS, X_COLS, matrix_x);

    printf("\n");

    // TODO

    // int numprocs, rank, namelen;
    // char processor_name[MPI_MAX_PROCESSOR_NAME];

    // MPI_Init(&argc, &argv);
    // MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Get_processor_name(processor_name, &namelen);

    // MPI_Finalize();

    return 0;
}
