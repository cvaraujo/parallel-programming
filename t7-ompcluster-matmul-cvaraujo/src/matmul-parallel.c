
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <omp.h>

#define BS 500

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

float* getBlockA(float *a, int i, int j, unsigned size){
  int k, l, fc = i*BS*size, q = 0;
  float *ba = (float *)malloc(sizeof(float) * (BS * BS));

  
  for(k = fc; k < (size * size) && k < (i+1)*BS*size; k+=size){
    for(l = k+j*BS; l < (size * size) && l < k+(j+1)*BS; l++){
      ba[q++] = a[l];
    }
  }
  return ba;
}

float* getBlockB(float *b, int i, int j, unsigned size){
  int k, l, fc = i*BS*size, q = 0;
  float *cb = (float *)malloc(sizeof(float) * (BS*BS));

  for(k = fc; k < (size * size) && k < (i+1)*BS*size; k+=size){
    for(l = k+j*BS; l < (size * size) && l < k+(j+1)*BS; l++){
      cb[q++] = b[l];
    }
  }
  return cb;
}

float* getBlockC(float *c, int i, int j, unsigned size){
  int k;
  float *ca = (float *)malloc(sizeof(float) * (BS*BS));
  for(k = 0; k < (BS*BS); k++)
    ca[k] = 0.0;
  return ca;
}


// Parallelize this function using OmpCluster
void multiply(float ***a, float ***b, float ***c, unsigned size) {
#pragma omp parallel
#pragma omp single
  {

    for (int i = 0; i < size/BS; i++) {
      for (int j = 0; j < size/BS; j++) {
	float *bc = c[i][j];
        for (int k = 0; k < size/BS; k++) {
          float *ba = a[i][k];
          float *bb = b[k][j];

#pragma omp target nowait depend(in: ba[0], bb[0]) depend(inout: bc[0]) map(to: ba[:BS*BS], bb[:BS*BS]) map(tofrom: bc[:BS*BS])
	  {
#pragma omp parallel for
	    {
	      for(int ii = 0; ii < BS; ii++) {
		for(int jj = 0; jj < BS; jj++) {
		  for(int kk = 0; kk < BS; kk++) {
		    bc[ii + jj * BS] += ba[kk + ii * BS] * bb[jj + kk * BS];
		  }
		}
	      }
	    }
	  }
        } 
      }
    }
  }
}

// Output matrix to stdout
void print_matrix(float ***c, unsigned size) {
  int k, l;
    for (int i = 0; i < size/BS; i++) {
      for (l = 0; l < BS; l++) {
        for (int j = 0; j < size/BS; j++) {
          for (k = l; k < (BS*BS); k+=BS) {
            printf(" %5.1f", c[i][j][k]);
          }
        }
        printf("\n");
      }
    }
}

int main(int argc, char *argv[]) {
  float *a, *b, *c, ***ba, ***bb, ***bc;
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

  int nb = size/BS;

  ba = (float ***)malloc(sizeof(float**) * nb);
  bb = (float ***)malloc(sizeof(float**) * nb);
  bc = (float ***)malloc(sizeof(float**) * nb);

  for (int i = 0; i < nb; i++) {
    ba[i] = (float **)malloc(sizeof(float*) * nb);
    bb[i] = (float **)malloc(sizeof(float*) * nb);
    bc[i] = (float **)malloc(sizeof(float*) * nb);
  }
  
  
  // initialize_matrices with random data
  initialize_matrices(a, b, c, size, seed);

  // initialize the block matrices
  for(int i = 0; i < nb; i++) {
    for(int j = 0; j < nb; j++) {
      ba[i][j] = getBlockA(a, i, j, size);
      bb[i][j] = getBlockB(b, i, j, size);
      bc[i][j] = getBlockC(c, i, j, size);
    }
  }

  // Multiply matrices
  t = omp_get_wtime();
  multiply(ba, bb, bc, size);
  t = omp_get_wtime() - t;

  // Show result
  print_matrix(bc, size);

  // Output elapsed time
  fprintf(stderr, "%lf\n", t);

  // Release memory
  free(a);
  free(b);
  free(c);

  return EXIT_SUCCESS;
}
