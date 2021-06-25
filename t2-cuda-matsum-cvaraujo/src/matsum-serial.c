#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char **argv) {
  int *A, *B, *C;
  int i, j;
  double t;

  // Input
  int rows, cols;
  FILE *input;

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file\n");
    return EXIT_FAILURE;
  }

  if ((input = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error: could not open file\n");
    return EXIT_FAILURE;
  }

  fscanf(input, "%d", &rows);
  fscanf(input, "%d", &cols);

  // Allocate memory on the host
  A = (int *)malloc(sizeof(int) * rows * cols);
  B = (int *)malloc(sizeof(int) * rows * cols);
  C = (int *)malloc(sizeof(int) * rows * cols);

  // Initialize memory
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      A[i * cols + j] = B[i * cols + j] = i + j;
    }
  }

  // Main kernel that sums two matrices element-wise
  t = omp_get_wtime();
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      C[i * cols + j] = A[i * cols + j] + B[i * cols + j];
    }
  }
  t = omp_get_wtime() - t;

  long long int sum = 0;

  // Keep this computation on the CPU
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      sum += C[i * cols + j];
    }
  }

  fprintf(stdout, "%lli\n", sum);
  fprintf(stderr, "%lf\n", t);

  free(A);
  free(B);
  free(C);
}
