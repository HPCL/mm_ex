/*
 * Matrix Multiply.
 *
 * This is a simple matrix multiply program which will compute the product
 *
 *                C  = A * B
 *
 * A ,B and C are both square matrix. They are statically allocated and
 * initialized with constant number, so we can focuse on the parallelism.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "mpi.h"

#ifndef ORDER
#define ORDER 1000   // the order of the matrix
#endif
#define AVAL  3.0    // initial value of A
#define BVAL  5.0    // initial value of B
#define TOL   0.001  // tolerance used to check the result

#define TYPE double

MPI_Status status;

// Initialize the matrices (uniform values to make an easier check)
void matrix_init(TYPE** A, TYPE** B, TYPE** C, int size) {
  int i, j;

  *A = (TYPE*)malloc(size*size*sizeof(TYPE));
  *B = (TYPE*)malloc(size*size*sizeof(TYPE));
  *C = (TYPE*)malloc(size*size*sizeof(TYPE));

  for (j=0; j<size*size; j++) {
    *A[j] = AVAL;
    *B[j] = BVAL;
    *C[j] = 0.0;
  }
 
}

void matrix_free(TYPE* A, TYPE* B, TYPE* C, int size) {
  int i;

  free(A);
  free(B);
  free(C);

}

// The actual mulitplication function, totally naive
double row_multiply(TYPE* A, int row_a, TYPE* B, int col_b, TYPE* out_buffer) {

  int k;
  double start, end;

  for (k=0; k<ORDER; k++){
    out_buffer[row_a*ORDER+col_b] += A[row_a*ORDER+k] * B[k*ORDER+col_b];
  }

}

int print_mat(TYPE* C) {
  
  int i, j;
  double e  = 0.0;
  double ee = 0.0;
  double v  = AVAL * BVAL * ORDER;

  for (i=0; i<ORDER; i++) {
    for (j=0; j<ORDER; j++) {
      printf("%f    ",C[i*ORDER+j]);
    }
    printf("\n\n");
  }

}

// Function to check the result, relies on all values in each initial
// matrix being the same
int check_result(TYPE* C) {
        int i, j;

        double e  = 0.0;
        double ee = 0.0;
        double v  = AVAL * BVAL * ORDER;

        for (i=0; i<ORDER; i++) {
                for (j=0; j<ORDER; j++) {
                        e = C[i*ORDER+j] - v;
                        ee += e * e;
                }
        }
        if (ee > TOL) {
                return 0;
        } else {
                return 1;
        }
}

// main function
int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int correct;
  int err = 0;
  double run_time;
  double mflops;
	
  int num_ranks, rank_id, rows_per_rank;
  int i,j,k,r;
  int master_rows, master_i;

  double start, end;

  TYPE *A, *B, *C, *out_buffer;

  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  rows_per_rank = ORDER/(num_ranks-1);

  if (rank_id == 0) { // master

    //int nt = omp_get_max_threads();
    printf("Available processes =  %d \n", num_ranks);

    // initialize the matrices
    printf("init...\n");
    master_rows =  ORDER%(num_ranks-1);
    master_i = ORDER*(ORDER-master_rows);
    A = (TYPE*)malloc(ORDER*ORDER*sizeof(TYPE));
    B = (TYPE*)malloc(ORDER*ORDER*sizeof(TYPE));
    C = (TYPE*)malloc(ORDER*ORDER*sizeof(TYPE));
    out_buffer = (TYPE*)malloc(master_rows*ORDER*sizeof(TYPE));

    for (j=0; j<ORDER*ORDER; j++) {
      A[j] = AVAL;
      B[j] = BVAL;
      C[j] = 0.0;
    }

    /*
    printf("\n\n");
    printf("rows_per_rank = %d\n", rows_per_rank);
    printf("order         = %d\n", ORDER);
    printf("num_ranks     = %d\n", num_ranks);
    printf("\n\n");
    print_mat(A);
    printf("\n\n");
    print_mat(B);
    printf("\n\n");
    print_mat(C);
    printf("\n\n");
    */

    printf("multiply...\n");
    start = omp_get_wtime();

    //send
    for(r = 1; r < num_ranks; r++){
      MPI_Send(&A[(r-1)*rows_per_rank*ORDER], rows_per_rank*ORDER, MPI_DOUBLE, r, 0, MPI_COMM_WORLD);
      MPI_Send(B, ORDER*ORDER, MPI_DOUBLE, r, 1, MPI_COMM_WORLD);  
    }

    //multiply
    for(i = ORDER-master_rows; i < ORDER; i++){
      for(j = 0; j < ORDER; j++){
        for(k = 0; k < ORDER; k++) {
          C[i*ORDER+j] += A[i*ORDER+k] * B[k*ORDER+j];
        }
      }
    }

    //recieve
    for(r = 1; r < num_ranks; r++){
      MPI_Recv(&C[(r-1)*rows_per_rank*ORDER], rows_per_rank*ORDER, MPI_DOUBLE, r, 2, MPI_COMM_WORLD, &status);
    }

    end = omp_get_wtime();

    // verify that the result is sensible
    printf("check...\n");
    correct  = check_result(C);

    //printf("\n\n");    
    //print_mat(C);
    //printf("\n\n");

    // Compute the number of mega flops
    run_time = end - start;
    mflops = (2.0 * ORDER*ORDER*ORDER) / (1000000.0 * run_time);
    printf("Order %d multiplication in %f seconds \n", ORDER, run_time);
    printf("Order %d multiplication at %f mflops\n", ORDER, mflops);

    // Check results
    if (! correct) {
      fprintf(stderr,"\n Errors in multiplication\n");
      err = 1;
    } else {
      fprintf(stdout,"\n SUCCESS : results match\n");
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    matrix_free(A,B,C,ORDER);
    free(out_buffer);

  } else { // workers
    
    //alloc
    out_buffer = (TYPE*)malloc(rows_per_rank*ORDER*sizeof(TYPE));
    A = (TYPE*)malloc(rows_per_rank*ORDER*sizeof(TYPE));
    B = (TYPE*)malloc(ORDER*ORDER*sizeof(TYPE));
    for(i = 0; i < rows_per_rank*ORDER; i++) out_buffer[i]=0;

    //recieve
    MPI_Recv(A, rows_per_rank*ORDER, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(B, ORDER*ORDER, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

    //multiply
    for(i = 0; i < rows_per_rank; i++){
      for(j = 0; j < ORDER; j++){
        row_multiply(A, i, B, j, out_buffer);
      }
    }

    //send
    MPI_Send(out_buffer, rows_per_rank*ORDER, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    free(out_buffer);
    free(A);
    free(B);

  }

  MPI_Finalize();

}

