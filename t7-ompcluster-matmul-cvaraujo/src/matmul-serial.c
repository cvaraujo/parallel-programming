
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <omp.h>

// Initialize matrices
void initialize_matrices(float *a, float *b, float *c, unsigned size,
                         unsigned seed) {
  srand(seed);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      a[i * size + j] = rand() % 10;
      b[i * size + j] = rand() % 10;
      c[i * size + j] = 0.0f;
    }
  }
}

void multiply(float *a, float *b, float *c, unsigned size) {
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      float sum = 0.0;
      for (int k = 0; k < size; ++k) {
        sum = sum + a[i * size + k] * b[k * size + j];
      }
      c[i * size + j] = sum;
    }
  }
}

// Output matrix to stdout
void print_matrix(float *c, unsigned size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf(" %5.1f", c[i * size + j]);
    }
    printf("\n");
  }
}

int main(int argc, char *argv[]) {
  float *a, *b, *c;
  unsigned seed, size;
  double t;
  FILE *input;

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file\n");
    return 1;
  }

  if ((input = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error: could not open file\n");
    return 1;
  }

  // Read inputs
  fscanf(input, "%u", &size);
  fscanf(input, "%u", &seed);

  // Allocate matrices
  a = (float *)malloc(sizeof(float) * size * size);
  b = (float *)malloc(sizeof(float) * size * size);
  c = (float *)malloc(sizeof(float) * size * size);

  // initialize_matrices with random data
  initialize_matrices(a, b, c, size, seed);

  // Multiply matrices
  t = omp_get_wtime();
  multiply(a, b, c, size);
  t = omp_get_wtime() - t;

  // Show result
  print_matrix(c, size);

  // Output elapsed time
  fprintf(stderr, "%lf\n", t);

  // Release memory
  free(a);
  free(b);
  free(c);

  return EXIT_SUCCESS;
}
