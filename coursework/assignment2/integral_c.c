#include <stdio.h>
#include <math.h>

double f(double x) { return 1 / x; }

int main(int argc, char *argv[]) {

    // Integral (1 to 6) of 1/x
    double integral;
    int a = 1, b = 6;
    int n = 2000;
    double h;
    double x;

    h = (double)(b-a)/n;

    integral = (f(a) + f(b)) / 2.0;
    x = a;

    for (int i = 1; i < n; i++) {
        x = x + h;
        integral = integral + f(x);
    }

    integral = integral * h;

    double expected = 1.79175946923;
    printf("Integral (%d to %d) of 1/x \n\n", a, b);
    printf("Expected  result  = %.11f\n", expected);
    printf("Estimated result  = %.11f\n", integral);
    printf("                    -------------\n");

    double difference = expected - integral;
    difference = fabsf(difference); // absolute value (no minus sign)
    printf("Difference        = %.11f\n", difference);

    return 0;
}
