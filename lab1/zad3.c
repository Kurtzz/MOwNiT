#include <stdio.h>

 int main()  {
    double epsilon = 1.0;

    printf( "epsilon;  1 + epsilon\n" );

    do {
       printf( "%G\t%.20f\n", epsilon, (1.0 + epsilon) );
       epsilon /= 2.0f;
    }
    // If next epsilon yields 1, then break
    while ((1.0 + (epsilon/2.0)) != 1.0); //

    // because current epsilon is the machine epsilon.
    printf( "\nCalculated Machine epsilon: %G\n", epsilon );
    return 0;
 }
