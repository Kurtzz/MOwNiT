#include <stdio.h>
#include <stdlib.h>

int main() {
    int i, count = 50;

    float f_x, f_tmp1, f_tmp2;
    double d_x, d_tmp1, d_tmp2;

    f_x = 0.01;
    d_x = 0.01;

    for(i = 1; i<=count; i++) {
        printf("%d: float: %f \t double: %f\n", i, f_x, d_x);
        f_tmp1 = 3.0 * f_x;
        f_tmp2 = 1 - f_x;
        f_x = f_x + f_tmp1 * f_tmp2;

        d_tmp1 = 3.0 * d_x;
        d_tmp2 = 1 - d_x;
        d_x = d_x + d_tmp1 * d_tmp2;
    }

    return 0;
}
