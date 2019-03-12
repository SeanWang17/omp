/* 
 * CS594 openmp homework
 * seq_dgemm.c
 * Qinglei Cao
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include<math.h>
#include <omp.h>

#define A(i, j) A[(i)+(j)*M]
#define B(i, j) B[(i)+(j)*K]
#define C(i, j) C[(i)+(j)*M]

/* 
 * C := alpha*op( A )*op( B ) + beta*C; op (X) = X or op (X) = X';
 * A is M*K matrix, B is K*N matrix, C is M*N matrix
 * In this homework, ignore transpose
 */ 
void xGEMM(char *tA, char *tB, int M, int N, int K, double alpha, 
           double *A, double *B, double beta, double *C){
  int i, j, k;
  for(i = 0; i < M; i++)
    for(j = 0; j < N; j++ )
      for(k = 0; k < K; k++)
        C(i, j) += alpha*A(i, k) * B(k, j);
}

int main(int argc, char *argv[]){
  int i, j, M, N, K;
  char *tran = "N";
  const double alpha = 1.0;
  const double beta = -1.0;
  double time_start, time_end, time;
  double diff =0.0;

  if(!((4 == argc) || (2 == argc))){
      printf("Usage: ./seq_dgemm M N K\n");
      printf("or Usage: ./seq_dgemm N\n");
      exit(1);
  } 

  if((4 == argc)){
      M = atoi(argv[1]);
      N = atoi(argv[2]);
      K = atoi(argv[3]);
   }else{
      M = atoi(argv[1]);
      N = atoi(argv[1]);
      K = atoi(argv[1]);
   }
  double flop = (double)M*(double)N*(double)K; /* inside the loop, one multiply */
  double precision = 2.22e-16*M*N*K;
  
  double *A = (double *) calloc(M*K, sizeof(double));
  double *B = (double *) calloc(K*N, sizeof(double));
  double *C = (double *) calloc(M*N, sizeof(double));
  char *result = (char *)malloc(10*sizeof(char));
  
  /* initialize input */
  for(i=0; i<M; i++)
    for(j=0; j<K; j++)
      A(i, j) = (double)rand()/(double)(RAND_MAX/5); 

  for(i=0; i<K; i++)
    for(j=0; j<N; j++)
      B(i, j) = (double)rand()/(double)(RAND_MAX/5);
  
  /* xGEMM */
  time_start = omp_get_wtime();
  xGEMM(tran, tran, M, N, K, alpha, A, B, beta, C);
  time_end = omp_get_wtime();
  time = time_end - time_start;

  /* BLAS */
  dgemm_(tran, tran, &M, &N, &K, &alpha, A, &N, B, &K, &beta, C, &M);

  /* check the correctness */
  for(i = 0; i < M; i++)
    for(j = 0; j < N; j++)
      diff = (diff > fabs(C(i, j)))? diff: fabs(C(i, j));

  if(diff <= precision)
    result = "correct";
  else
    result = "wrong";

  printf("result: %s: %d, %e, %lf\n", result, M, diff, 1.0e-9*flop/time);

  free(A);
  free(B);
  free(C);
  
  return 0;
}
