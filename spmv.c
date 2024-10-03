#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
//#include <suitesparse/cs.h>
#include <string.h>
#include <unistd.h>
#include "timer.h"
#include "spmv.h"

#define DEFAULT_SIZE 1024
#define DEFAULT_DENSITY 0.25

double *generate_sparse_matrix(unsigned int seed, unsigned int size, double density)
{
  srand(seed);
  double *matrix = (double *) calloc(size * size, sizeof(double));

  for (unsigned int i = 0; i < size * size; i++) {
    if ((rand() % 100) / 100.0 < density) {
      // Get a pseudorandom value between -9.99 e 9.99
      matrix[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }
  }

  return matrix;
}

double *generate_vector(unsigned int seed, unsigned int size)
{
  srand(seed);
  double *vec = (double *) malloc(sizeof(double) * size);

  for (unsigned int i = 0; i < size; i++) {
    vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
  }

  return vec;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return fabs(x - y) <= epsilon * fabs(x);
  // see Knuth section 4.2.2 pages 217-218
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    if (!is_nearly_equal(ref[i], result[i]))
      return 0;
  }

  return 1;
}

int main(int argc, char *argv[])
{
  int size;        // number of rows and cols (size x size matrix)
  double density;  // aprox. ratio of non-zero values

  if (argc < 2) {
    size = DEFAULT_SIZE;
    density = DEFAULT_DENSITY;
  } else if (argc < 3) {
    size = atoi(argv[1]);
    density = DEFAULT_DENSITY;
  } else {
    size = atoi(argv[1]);
    density = atoi(argv[2]);
  }

  double *mat, *vec, *refsol, *mysol;

  mat = generate_sparse_matrix(1, size, density);
  vec = generate_vector(2, size);
  refsol = (double *) malloc(size * sizeof(double));
  mysol = (double *) malloc(size * sizeof(double));

  //
  // Dense computation using LAPACK/OpenBLAS
  //
  timeinfo start, now;
  timestamp(&start);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

  timestamp(&now);
  printf("Time taken by Lapacke/OpenBLAS dense computation: %ld ms\n", diff_milli(&start, &now));


  //
  // Using your own dense implementation
  //
  timestamp(&start);

  my_dense(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  //
  // Let's try now SpmV: Sparse Matrix - Dense Vector computation
  //

  // Convert mat to a sparse format: CSR

  //
  // Sparse computation using CSPARSE
  //
  timestamp(&start);

  CSPARSE(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix - vector product: %ld ms\n", diff_milli(&start, &now));

  //
  // Your own sparse implementation
  //

  // Compare times (and computation correctness!)

  return 0;
}
