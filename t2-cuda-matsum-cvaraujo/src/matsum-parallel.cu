#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

const int THR = 1024;

__global__ void matrix_sum(int *C, int *A, int *B, int n) {
	   int index = blockIdx.x * blockDim.x + threadIdx.x;
	   if (index < n) C[index] = A[index] + B[index];	   
}

int main(int argc, char **argv) {
  int *h_A, *h_B, *h_C;
  int *d_A, *d_B, *d_C;
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
  h_A = (int *)malloc(sizeof(int) * rows * cols); //host_init(h_A);
  h_B = (int *)malloc(sizeof(int) * rows * cols); //host_init(h_B);
  h_C = (int *)malloc(sizeof(int) * rows * cols); //host_init(h_C);

  // Initialize memory
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      h_A[i * cols + j] = h_B[i * cols + j] = i + j;
    }
  }

  int bytes = sizeof(int) * rows * cols;
  // Copy data to device
  // ...
  cudaMalloc(&d_A, bytes);
  cudaMalloc(&d_B, bytes);
  cudaMalloc(&d_C, bytes);

  cudaMemcpy(d_A, h_A, bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, bytes, cudaMemcpyHostToDevice);

  // Compute matrix sum on device
  // Leave only the kernel and synchronize inside the timing region!
  int N = rows * cols;

  t = omp_get_wtime();
  matrix_sum<<<(N + THR-1)/THR, THR>>>(d_C, d_A, d_B, rows * cols);	
  cudaDeviceSynchronize();
  t = omp_get_wtime() - t;
  

  // Copy data back to host
  // ...
  cudaMemcpy(h_C, d_C, bytes, cudaMemcpyDeviceToHost);

  long long int sum = 0;

  // Keep this computation on the CPU
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      sum += h_C[i * cols + j];
    }
  }

  fprintf(stdout, "%lli\n", sum);
  fprintf(stderr, "%lf\n", t);

  cudaFree(d_A); free(h_A);
  cudaFree(d_B); free(h_B);
  cudaFree(d_C); free(h_C);

  return 0;
}
