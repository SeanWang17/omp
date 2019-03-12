/* 
 * CS594 openmp homework
 * omp_dtrsm.c
 * Qinglei Cao
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include<math.h>
#include <omp.h>

#define A(i, j) A[(i)*M+(j)]
#define X(i, j) X[(i)*N+(j)]
#define B(i, j) B[(i)*N+(j)]
#define B_mul(i, j) B_mul[(i)*N+(j)]

/* 
 * solve A * X =  B, A is a lower or upper triangular matrix;
 * A is M*M matrix, X is M*N matrix, B is M*N matrix
 * In this homework, ignore transpose
 */ 

/* forward substitution, A is lower triangular matrix */
void forward(int M, int N, double *A, double *X, double *B){
  int i, j, k, tid, num_thread;
  double temp;

  printf("lower triangular matrix\n");
#pragma omp parallel private(i, j, k, tid) 
  {
    tid = omp_get_thread_num();
    num_thread = omp_get_num_threads();
    //printf("tid:%d, num:%d\n", tid, num_thread);
    for(j=0+tid; j<N; j+=num_thread){
       X[j] = B[j]/A[0];
       for(i=1; i<M; i++){
          temp = 0.0;
          for(k=0; k<i; k++)
             temp += A(i, k) * X(k, j);

          X(i, j) = (B(i, j) - temp)/A(i, i);
       }
    }
  }
}

/* backward substitution, A is upper triangular matrix */
void backward(int M, int N, double *A, double *X, double *B){
  int i, j, k, tid, num_thread;
  double temp;

  printf("upper triangular matrix\n");
#pragma omp parallel private(i, j, k, tid) 
  {
    tid = omp_get_thread_num();
    num_thread = omp_get_num_threads();
    //printf("tid:%d, num:%d\n", tid, num_thread);
    for(j=0+tid; j<N; j+=num_thread){
       X(N-1, j) = B(N-1, j)/A(N-1, N-1);

       for(i=M-2; i>=0; i--){
          temp = 0.0;
          for(k=i+1; k<M; k++){
             temp += A(i, k) * X(k, j);
          }
          X(i, j) = (B(i, j)-temp)/A(i, i);
       }
    }
  }
}

void print_m(int M, int N, double *A){
  int i, j;
  for(i=0; i<M; i++){
      for(j=0; j<N; j++)
           printf("%e  ", A[i*N+j]);
      printf("\n");
  }
  printf("\n");
}

int main(int argc, char *argv[]){
  int i, j, k, M, N;
  char *uplo = "L";
  double time_start, time_end, time;

  if(!((4 == argc) || (3 == argc))){
      printf("Usage: ./seq_dgemm M N [U/L]\n");
      printf("or Usage: ./seq_dgemm N [U/L]\n");
      exit(1);
  } 

  if((4 == argc)){
      M = atoi(argv[1]);
      N = atoi(argv[2]);
      uplo = argv[3];
   }else{
      M = atoi(argv[1]);
      N = atoi(argv[1]);
      uplo = argv[2];
   }

  double flop = (double)M*(double)M*(double)N/2.0; /* inside the loop, one multiply */
  
  double *A = (double *) calloc(M*M, sizeof(double));
  double *X = (double *) calloc(M*N, sizeof(double));
  double *B = (double *) calloc(M*N, sizeof(double));
  double *B_mul = (double *) calloc(M*N, sizeof(double));
  char *result = (char *)malloc(10*sizeof(char));

  /* initialize input */
  if(!strcmp(uplo, "L")){
     for(i=0; i<M; i++)
        for(j=0; j<=i; j++)
           A(i, j) = (double)rand()/(double)(RAND_MAX/5); 
   }else if(!strcmp(uplo, "U")){ 
     for(i=0; i<M; i++)
        for(j=i; j<M; j++)
           A(i, j) = (double)rand()/(double)(RAND_MAX/5);
   }else{
 	printf("You should give L or U to set matrix\n");
	exit(1);
   }

   for(i=0; i<M; i++)
     for(j=0; j<N; j++)
        B(i, j) = (double)rand()/(double)(RAND_MAX/5);

  /* my dtrsm */
  time_start = omp_get_wtime();

  if(!strcmp(uplo, "L"))
     forward(M, N, A, X, B);
  else
     backward(M, N, A, X, B);

  time_end = omp_get_wtime();
  time = time_end - time_start;

  /* check the correctness */
  /* get A * X */
  for(i = 0; i < M; i++)
    for(j = 0; j < N; j++ )
      for(k = 0; k < M; k++)
        B_mul(i, j) += A(i, k) * X(k, j);

/*
  puts("A:"); print_m(M, M, A);
  puts("X:"); print_m(M, N, X);
  puts("B:"); print_m(M, N, B);
  puts("B_mul:"); print_m(M, N, B_mul);
*/

  printf("result: %d, %lf\n", M, 1.0e-9*flop/time);

  free(A);
  free(B);
  free(B_mul);
  free(X);
  
  return 0;
}
