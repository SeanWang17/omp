#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include<math.h>
#include <omp.h>
#include <mkl.h>
#define A(i, j) A[(i)+(j)*M]
#define B(i, j) B[(i)+(j)*K]
#define C(i, j) C[(i)+(j)*M]


void xGEMM(char *tA, char *tB, int M, int N, int K, double alpha, 
           double *A, double *B, double beta, double *C){
  int i, j, k;
  #pragma omp parallel for collapse(2) private(i,j,k)
  for(i = 0; i < M; i++)
    for(j = 0; j < N; j++ )
      for(k = 0; k < K; k++)
        C(i, j) += alpha*A(i, k) * B(k, j);
}

  int main(int argc, char *argv[]){
  
  int i, j, M, N, K, tNum;
  tNum = 5;
  char *tran = "N";
  const double alpha = 1.0;
  const double beta = -1.0;
  double time_start, time_end, time;
  double diff =0.0;

  if(!((4 == argc) || (2 == argc) || (3 == argc))){
      printf("Usage: ./seq_dgemm M N K\n");
      printf("or Usage: ./seq_dgemm N\n");
      exit(1);
  } 
  if (3==argc){
      M = atoi(argv[1]);
      N = atoi(argv[1]);
      K = atoi(argv[1]);
      tNum = atoi(argv[2]);
}
  else if((4 == argc)){
      M = atoi(argv[1]);
      N = atoi(argv[2]);
      K = atoi(argv[3]);
   }else{
      M = atoi(argv[1]);
      N = atoi(argv[1]);
      K = atoi(argv[1]);
   }
  omp_set_num_threads(tNum);
  double flop = (double)M*(double)N*(double)K;
  double precision = 2.22e-16*M*N*K;
  
  double *A = (double *) calloc(M*K, sizeof(double));
  double *B = (double *) calloc(K*N, sizeof(double));
  double *C = (double *) calloc(M*N, sizeof(double));
  char *result = (char *)malloc(10*sizeof(char));
  
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

  printf("result: %s, threadNum %d Performance: %lf\n", result, tNum, 1.0e-9*flop/time);
  //printf("result %s size %d performance %lf\n", result, M, 1.0e-9*flop/time);
  free(A);
  free(B);
  free(C);
  
  return 0;
}
