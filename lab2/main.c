#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

int N, i,j,k;

int main(int argc, char **argv){
  float a = 100.0;
  float tmp;
  clock_t start, middle, stop;

  if(argc != 2) {
    printf("Niepoprawna liczba argumentów!\n");
    exit(0);
  }

  N = atoi(argv[1]);
  gsl_vector *vector1 = gsl_vector_calloc(N);
  gsl_vector *vector2 = gsl_vector_calloc(N);
  gsl_matrix *matrix = gsl_matrix_calloc(N, N);
  gsl_vector *result = gsl_vector_calloc(N);
  CBLAS_TRANSPOSE_t TransA = CblasNoTrans;

  double res = 0;

  srand(time(NULL));
  for(j = 0; j < N; j++){
    tmp = ((float)rand()/(float)(RAND_MAX)) * a;
    gsl_vector_set(vector1, j, tmp);

    tmp = ((float)rand()/(float)(RAND_MAX)) * a;
    gsl_vector_set(vector2, j, tmp);

    for(k = 0; k < N; k++){
      tmp = ((float)rand()/(float)(RAND_MAX)) * a;
      gsl_matrix_set(matrix, j, k, tmp);
    }
  }

  FILE *fp;
  fp = fopen("result.txt", "a");

  for(i = 0; i < 10; i++){
    //NEW DATA
    start = clock();
    gsl_blas_ddot(vector1, vector2, &res);
    middle = clock();
    gsl_blas_dgemv(TransA, 1, matrix, vector1, 1, result);
    stop = clock();
    fprintf(fp, "%d,b1,%f\n",N,(double)(middle - start)*1000 / (double)CLOCKS_PER_SEC);
    fprintf(fp, "%d,b2,%f\n",N,(double)(stop - middle)*1000 / (double)CLOCKS_PER_SEC);
  }
  gsl_vector_free(vector1);
  gsl_vector_free(vector2);
  gsl_vector_free(result);
  gsl_matrix_free(matrix);
  fclose(fp);
  return 0;
}
