#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

int i,j,k,N,tmp,l;
clock_t start;
void naive(int**, int**, int**);
void ver2(int **, int **, int **);

int main(int argc, char **argv){
  FILE *fp;
  if(argc != 2) {
    printf("Niepoprawna liczba argumentów!\n");
    exit(0);
  }
  N = atoi(argv[1]);

  int **A, **B, **C;
  A = calloc(N, sizeof(int*));
  B = calloc(N, sizeof(int*));
  C = calloc(N, sizeof(int*));
  for(i = 0; i<N; i++){
    A[i] = calloc(N, sizeof(int));
    B[i] = calloc(N, sizeof(int));
    C[i] = calloc(N, sizeof(int));
  }

  gsl_matrix *matrix1 = gsl_matrix_calloc(N, N);
  gsl_matrix *matrix2 = gsl_matrix_calloc(N, N);
  gsl_matrix *result = gsl_matrix_calloc(N, N);
  CBLAS_TRANSPOSE_t TransA = CblasNoTrans;

  srand(time(NULL));
  fp = fopen("result.txt", "a");

  //fill matrix
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      tmp = rand() % 100;
      A[i][j] = tmp;
      gsl_matrix_set(matrix1, j, k, tmp);

      tmp = rand() % 100;
      B[i][j] = tmp;
      gsl_matrix_set(matrix2, j, k, tmp);
    }
  }

  for(l = 0; l < 20; l++){
    //algorithms
    start = clock();
    naive(A,B,C);
    fprintf(fp, "%d,alg1,%f\n", N, (double)(clock() - start)/CLOCKS_PER_SEC);
    //printf("Algorytm 1: %g [s]\n", (double)(clock() - start)/CLOCKS_PER_SEC);

    start = clock();
    ver2(A,B,C);
    fprintf(fp, "%d,alg2,%f\n", N, (double)(clock() - start)/CLOCKS_PER_SEC);
    //printf("Algorytm 2: %g [s]\n", (double)(clock() - start)/CLOCKS_PER_SEC);

    start = clock();
    gsl_blas_dgemm (TransA, TransA, 1, matrix1, matrix2, 1, result);
    fprintf(fp, "%d,blas,%f\n", N, (double)(clock() - start)/CLOCKS_PER_SEC);
    //printf("Algorytm 3: %g [s]\n", (double)(clock() - start)/CLOCKS_PER_SEC);
  }

  return 0;
}

void naive(int **A, int **B, int **C){
  int i,j,k;
  for (i=0;i<N; i++)
    for(j=0;j<N;j++)
      for(k=0;k<N; k++)
        C[i][j]+=A[i][k]*B[k][j];
}

void ver2(int **A, int **B, int **C) {
  int i,j,k;
  for (i=0;i<N; i++)
    for(k=0;k<N;k++)
      for(j=0;j<N; j++)
        C[i][j]+=A[i][k]*B[k][j];
}
