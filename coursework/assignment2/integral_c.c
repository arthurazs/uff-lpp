#include <stdio.h>
#include <math.h>

float f(float x) {
    return 1/x;
}

int main(int argc, char *argv[]) {

    // Integral (1 to 6) of 1/x
    float integral;
    float a = 1, b = 6;
    int n = 2000;
    float h;
    float x;

    h = (b-a)/n;

    integral = (f(a) + f(b)) / 2.0;
    x = a;

    for (int i = 1; i <= n - 1; i++) {
        x = x + h;
        integral = integral + f(x);
    }

    integral = integral * h;

    float expected = 1.79175946923;
    printf("Integral (%.0f to %.0f) of 1/x \n\n", a, b);
    printf("Expected  result  = %.11f\n", expected);
    printf("Estimated result  = %.11f\n", integral);
    printf("                    -------------\n");

    float difference = expected - integral;
    difference = fabsf(difference); // absolute value (no minus sign)
    printf("Difference        = %.11f\n", difference);

    return 0;
}
