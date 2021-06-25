
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <omp.h>

// Initialize matrices
void initialize_matrices(float *a, float *b, float *c, float *r1, float *r2,
                         unsigned size, unsigned seed) {
  srand(seed);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      a[i * size + j] = rand() % 10;
      b[i * size + j] = rand() % 10;
      c[i * size + j] = rand() % 10;
      r1[i * size + j] = 0.0f;
      r2[i * size + j] = 0.0f;
    }
  }
}

void multiply(float *a, float *b, float *r, unsigned size) {
  float sum = 0.0;
  
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      sum = 0.0;
      for (int k = 0; k < size; ++k) {
        sum = sum + a[i * size + k] * b[k * size + j];
      }
      r[i * size + j] = sum;
    }
  }
}

void sum(float *a, float *b, float *r, unsigned size) {
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      r[i * size + j] = a[i * size + j] + b[i * size + j];
    }
  }
}

// Output matrix to stdout
void print_matrix(float *r, unsigned size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf(" %5.1f", r[i * size + j]);
    }
    printf("\n");
  }
}

int main(int argc, char *argv[]) {
  float *a, *b, *c, *r1, *r2;
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
  r1 = (float *)malloc(sizeof(float) * size * size);
  r2 = (float *)malloc(sizeof(float) * size * size);

  // initialize_matrices with random data
  initialize_matrices(a, b, c, r1, r2, size, seed);

  // Compute R = (A * B) + C
  t = omp_get_wtime();

  {
    // Check for parallelization on this block

    // r1 = a * b
    multiply(a, b, r1, size);

    // r2 = r1 + c
    sum(r1, c, r2, size);
  }

  t = omp_get_wtime() - t;

  // Show result
  print_matrix(r2, size);

  // Output elapsed time
  fprintf(stderr, "%lf\n", t);

  // Release memory
  free(a);
  free(b);
  free(c);
  free(r1);
  free(r2);

  return EXIT_SUCCESS;
}
