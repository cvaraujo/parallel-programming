
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
      b[j * size + i] = rand() % 10;
      c[i * size + j] = rand() % 10;
      r1[i * size + j] = c[i * size + j];
      r2[i * size + j] = 0.0f;
    }
  }
}

void multiply(float *a, float *b, float *c,  float *r, unsigned size) {
  float sum = 0.0;
  int i, j, k, iaux, jaux;

#pragma omp parallel default(none) shared(a, b, c, r, size) reduction(+: sum) private(i, j, k, iaux, jaux) num_threads(4) 
   {
 #pragma omp for schedule(dynamic)
    for (i = 0; i < size; ++i) {
      iaux = i * size;
        
        for (j = 0; j < size; ++j) {
        jaux = j * size; 
        sum = 0.0;
        
        for (k = 0; k < size; ++k) {
            sum = sum + a[iaux + k] * b[jaux + k];
        }
        r[iaux + j] = c[iaux + j] + sum;
       }
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

  // Do not change this line
  omp_set_num_threads(4);

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
    // r1 = (a * b) + c
    multiply(a, b, c, r1, size); 
  }

  t = omp_get_wtime() - t;

  // Show result
  print_matrix(r1, size);

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
