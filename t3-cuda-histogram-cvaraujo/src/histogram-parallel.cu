#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define COMMENT "Histogram_GPU"
#define RGB_COMPONENT_COLOR 255

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define BLOCK_SIZE 32

void check_cuda(cudaError_t error, const char *filename, const int line)
{
  if (error != cudaSuccess) {
    fprintf(stderr, "Error: %s:%d: %s: %s\n", filename, line,
                 cudaGetErrorName(error), cudaGetErrorString(error));
    exit(EXIT_FAILURE);
  }
}

#define CUDACHECK(cmd) check_cuda(cmd, __FILE__, __LINE__)

typedef struct {
  unsigned char red, green, blue;
} PPMPixel;

typedef struct {
  int x, y;
  PPMPixel *data;
} PPMImage;

static PPMImage *readPPM(const char *filename) {
  char buff[16];
  PPMImage *img;
  FILE *fp;
  int c, rgb_comp_color;
  fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open file '%s'\n", filename);
    exit(1);
  }

  if (!fgets(buff, sizeof(buff), fp)) {
    perror(filename);
    exit(1);
  }

  if (buff[0] != 'P' || buff[1] != '6') {
    fprintf(stderr, "Invalid image format (must be 'P6')\n");
    exit(1);
  }

  img = (PPMImage *)malloc(sizeof(PPMImage));
  if (!img) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  c = getc(fp);
  while (c == '#') {
    while (getc(fp) != '\n')
      ;
    c = getc(fp);
  }

  ungetc(c, fp);
  if (fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
    fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
    exit(1);
  }

  if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
    fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
    exit(1);
  }

  if (rgb_comp_color != RGB_COMPONENT_COLOR) {
    fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
    exit(1);
  }

  while (fgetc(fp) != '\n')
    ;
  img->data = (PPMPixel *)malloc(img->x * img->y * sizeof(PPMPixel));

  if (!img) {
    fprintf(stderr, "Unable to allocate memory\n");
    exit(1);
  }

  if (fread(img->data, 3 * img->x, img->y, fp) != img->y) {
    fprintf(stderr, "Error loading image '%s'\n", filename);
    exit(1);
  }

  fclose(fp);
  return img;
}

__global__ void histogram_kernel(PPMPixel* data, int rows, int cols, float* h) {
  // Local variables
  long long tid = threadIdx.x + blockDim.x * blockIdx.x;  
  int stride = blockDim.x * gridDim.x;
  int i = threadIdx.x;
  long long n = rows * cols;

  __shared__ int temp[64];
  temp[i%64] = 0;
  
  __syncthreads();

  while (tid < n) {
    atomicAdd(&(temp[int(data[tid].red * 16 + data[tid].green * 4  + data[tid].blue)]), 1);
    tid += stride;
  }
  __syncthreads();
  
  atomicAdd(&(h[i%64]), temp[i%64]);

}

double Histogram(PPMImage *image, float *h_h) {
  float ms;
  cudaEvent_t start, stop;
  PPMPixel *d_pixels;
  float *h_d;
  long long  n = image->y * image->x;
  int i;

  for (i = 0; i < n; i++) {
    image->data[i].red = floor((image->data[i].red * 4) / 256);
    image->data[i].green = floor((image->data[i].green * 4) / 256);
    image->data[i].blue = floor((image->data[i].blue * 4) / 256);
  }

  cudaMalloc(&d_pixels, sizeof(PPMPixel) * n);
  cudaMalloc(&h_d, sizeof(float)*64);

  cudaMemcpy(d_pixels, image->data, n * sizeof(PPMPixel), cudaMemcpyHostToDevice);
  cudaMemcpy(h_d, h_h, 64 * sizeof(float), cudaMemcpyHostToDevice);

  // Create Events
  CUDACHECK(cudaEventCreate(&start));
  CUDACHECK(cudaEventCreate(&stop));
  
  int THR = 64;
  int BLO = (n + (THR-1))/THR;
  
  // Launch kernel and compute kernel runtime.
  // Warning: make sure only the kernel is being profiled, memcpies should be
  // out of this region.
  size_t t_size = 64 * sizeof(float);

  CUDACHECK(cudaEventRecord(start));
  histogram_kernel<<<BLO, THR,t_size>>>(d_pixels, image->x, image->y, h_d);
  CUDACHECK(cudaEventRecord(stop));
  CUDACHECK(cudaEventSynchronize(stop));
  CUDACHECK(cudaEventElapsedTime(&ms, start, stop));

  cudaMemcpy(h_h, h_d, 64 * sizeof(float), cudaMemcpyDeviceToHost);

  // Destroy events
  CUDACHECK(cudaEventDestroy(start));
  CUDACHECK(cudaEventDestroy(stop));
  cudaFree(h_d);

  return ((double)ms) / 1000.0;
}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file\n");
    return 1;
  }

  PPMImage *image = readPPM(argv[1]);
  float *h = (float *)malloc(sizeof(float) * 64);
  long long  n = image->y * image->x;

  // Initialize histogram
  for (int i = 0; i < 64; i++)
    h[i] = 0.0;
  
  // Compute histogram
  double t = Histogram(image, h);

  for (int i = 0; i < 64; i++)
    printf("%0.3f ", h[i]/n);
  printf("\n");

  fprintf(stderr, "%lf\n", t);
  free(h);
}
