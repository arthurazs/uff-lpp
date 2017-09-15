// mat_x num of rows equals to
// mat_a num of rows
#define A_ROWS 3
#define X_ROWS A_ROWS
// mat_x num of cols equals to
// mat_b num of cols
#define B_COLS 3
#define X_COLS B_COLS
// mat_a num of cols should be equals to
// mat_b num of rows
#define A_COLS 2
#define B_ROWS A_COLS

#include <stdio.h>
#include <mpi.h>

void print_matrix(char* name, int rows, int cols, int matrix[rows][cols]) {

    printf("\n%s [%d][%d]\n", name, rows, cols);
    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++)
            printf("%d ", matrix[row][col]);
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

    int numprocs, rank, tag = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int master = X_ROWS;

    if (numprocs == X_ROWS + 1) {

        int matrix_a [A_ROWS][A_COLS] = {
            {9, 0},
            {5, 6},
            {1, 2}
        };

        int matrix_b [B_ROWS][B_COLS] = {
            {2, 4, 3},
            {7, 8, 9}
        };

        int matrix_x [X_ROWS][X_COLS] = {0};

        if (rank < X_ROWS) {
            int line[X_COLS];
            for (int col = 0; col < B_COLS; col++) {
                printf("P(%1$d): Multiplying row %1$d with column %2$d\n", rank, col);
                line[col] = 0;
                for (int ctrl = 0; ctrl < B_ROWS; ctrl++)
                    line[col] +=
                        matrix_a[rank][ctrl] * matrix_b[ctrl][col];
            }
            MPI_Send(
                line, X_COLS, MPI_INT, master, tag, MPI_COMM_WORLD);
        }
        else {
            printf("P(%d): Waiting multiplication results\n", rank);
            for (int origin = 0; origin < master; origin++) {
                int line[X_COLS];
                MPI_Recv(line, X_COLS, MPI_INT, origin, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int col = 0; col < X_COLS; col++)
                    matrix_x[origin][col] = line[col];
                printf("P(%d): Received message from P(%d)\n", rank, origin);
            }
            print_matrix("Matrix A", A_ROWS, A_COLS, matrix_a);
            print_matrix("Matrix B", B_ROWS, B_COLS, matrix_b);
            print_matrix("Matrix X", X_ROWS, X_COLS, matrix_x);
            printf("\n");
        }
    }
    else {
        if (rank == 0){
            printf("\nERROR\n");
            printf("This code needs exactly %d process slots (running %d right now)\n", X_ROWS + 1, numprocs);
            printf("Try: mpirun -n %d\n\n", X_ROWS + 1);
        }
    }

    MPI_Finalize();

    return 0;
}
