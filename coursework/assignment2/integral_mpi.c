#include <stdio.h>
#include <mpi.h>
#include <math.h>

float f(float x) {
    return 1/x;
}

float calculate(float local_a, float local_b, int local_n, float h) {

    float x = local_a;
    float integral = (f(local_a) + f(local_b)) / 2.0;

    for (int i = 1; i <= local_n - 1; i++) {
        x = x + h;
        integral = integral + f(x);
    }

    integral = integral * h;

    return integral;
}

int main(int argc, char *argv[]) {

    int numprocs, rank, tag = 0, master = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Integral (1 to 6) of 1/x
    float a = 1, b = 6, h;
    int trapezoids = 2000;
    h = (b-a)/trapezoids;

    // Local Variables
    //TODO Fix local_n (it HAS to be = 2000, not < 2000)
    float local_n = trapezoids / numprocs;
    float local_a = a + rank * local_n * h;
    float local_b = local_a + local_n * h;

    // Calculate
    float integral_local = calculate(local_a, local_b, local_n, h);
    float integral;

    MPI_Reduce(&integral_local, &integral, 1, MPI_FLOAT, MPI_SUM, master, MPI_COMM_WORLD);

    if (rank == master) {
        float expected = 1.79175946923;

        printf("Integral (%.0f to %.0f) of 1/x \n\n", a, b);
        printf("Expected  result  = %.11f\n", expected);
        printf("Estimated result  = %.11f\n", integral);
        printf("                    -------------\n");

        float difference = expected - integral;
        difference = fabsf(difference); // absolute value (no minus sign)
        printf("Difference        = %.11f\n", difference);
    }

    MPI_Finalize();

    return 0;
}
