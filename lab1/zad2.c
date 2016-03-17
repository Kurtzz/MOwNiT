#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>

int main (void) {
	int count = 130;

  float f = 1.0/3.0;
  double d = 1.0/3.0;

  double fd = f; /* promote from float to double */

  int i;

  for(i = 0; i<count; i++) {

  	f /= 2.0;
  	d /= 2.0;
  	fd = f;

	  printf(" f="); gsl_ieee_printf_float(&f);
	  printf("\n");

	  // printf("fd="); gsl_ieee_printf_double(&fd);
	  // printf("\n");
    //
	  // printf(" d="); gsl_ieee_printf_double(&d);
	  // printf("\n\n");

	}

  return 0;
}
