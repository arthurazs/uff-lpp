#include <stdio.h>
#include <math.h>
#include <mpi.h>

double f(double x) { return 1 / x; }

double calculate(double local_a, double local_b, int local_n, double h){

    double x = local_a;
    double integral = (f(local_a) + f(local_b)) / 2.0;

    for (int i = 1; i < local_n; i++) {
        x = x + h;
        integral = integral + f(x);
    }

    integral = integral * h;
    return integral;
}

int main(int argc, char *argv[]) {

    int numprocs, rank, master = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Integral (1 to 6) of 1/x
    int a = 1, b = 6;
    double h;
    int trapezoids = 1237;
    h = (double)(b-a)/trapezoids;

    // Local Variables
    int local_n = trapezoids / numprocs;  //if trapezoids not divisible
    int extra = trapezoids % numprocs;    //by numprocs, split the rest
    if (extra > rank) {                   //to other processes.
        local_n++;                        //first few processes get the
    }                                     //extra elements

    double local_a = a + rank * local_n * h;
    if (extra != 0 && extra <= rank) {
        //recalculate local_a
        //if there are processes
        //with extra elements
        local_a = local_a + (extra * h);
    }
    double local_b = local_a + local_n * h;

    // Calculate
    double integral_local = calculate(local_a, local_b, local_n, h);
    printf("P(%d): Estimated %f for a(%.2f) to b(%.2f)\n", rank, integral_local, local_a, local_b);
    double integral;

    MPI_Reduce(
        &integral_local, &integral,
        1, MPI_DOUBLE, MPI_SUM,
        master, MPI_COMM_WORLD);

    if (rank == master) {
        double expected = 1.79175946923;

        printf("\nP(%d): MPI_Reduce\nIntegral (%d to %d) of 1/x\n\n", rank, a, b);

        printf("Expected  result  = %.11f\n", expected);
        printf("Estimated result  = %.11f\n", integral);
        printf("                    -------------\n");

        double difference = expected - integral;
        difference = fabsf(difference); //absolute value (no minus sign)
        printf("Difference        = %.11f\n", difference);
    }

    MPI_Finalize();

    return 0;
}
