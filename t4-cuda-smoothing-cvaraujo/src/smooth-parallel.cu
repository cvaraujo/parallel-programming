#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MASK_WIDTH 15
#define TILE_WIDTH 16

#define COMMENT "Histogram_GPU"
#define RGB_COMPONENT_COLOR 255

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

void writePPM(PPMImage *img) {

  fprintf(stdout, "P6\n");
  fprintf(stdout, "# %s\n", COMMENT);
  fprintf(stdout, "%d %d\n", img->x, img->y);
  fprintf(stdout, "%d\n", RGB_COMPONENT_COLOR);

  fwrite(img->data, 3 * img->x, img->y, stdout);
  fclose(stdout);
}


__global__ void smoothing_kernel(PPMPixel *image, PPMPixel *image_copy, int rows, int cols) {
    long long int col = blockIdx.x * blockDim.x + threadIdx.x;
    long long int row = blockIdx.y * blockDim.y + threadIdx.y;

    int border = (int)((MASK_WIDTH - 1) / 2);
    const int dimension = TILE_WIDTH + (MASK_WIDTH - 1);

    __shared__ PPMPixel shared[dimension+1][dimension+1];

    int size = ceil((float)(dimension * dimension) / (blockDim.x * blockDim.y));

    int block_start_col = blockIdx.x * blockDim.x;
    int block_start_row = blockIdx.y * blockDim.y;

    for(int k = 0; k < size; k++){

    	int index = (threadIdx.y * blockDim.x) + threadIdx.x + (k * blockDim.x * blockDim.y);

        int x = (int)(index / dimension);
        int y = (index % dimension);

        if(x < dimension && y < dimension){
            int sr = block_start_row + x - border;
            int sc = block_start_col + y - border;

            if(sr >= 0 && sc >= 0 && sc < cols && sr < rows) shared[x][y] = image_copy[sr * cols + sc];
            else shared[x][y].red = shared[x][y].green = shared[x][y].blue = 0;
            
        }
    }

    __syncthreads();

    int total_red, total_blue, total_green;
    total_red = total_blue = total_green = 0;
    //Determiando a posicao da matriz
    if(row < rows && col < cols){
        for(int i = threadIdx.y; i < (threadIdx.y + MASK_WIDTH); i++){
            for(int j = threadIdx.x; j < (threadIdx.x + MASK_WIDTH); j++){
                total_red += shared[i][j].red;
                total_green += shared[i][j].green;
                total_blue += shared[i][j].blue;
            } 
        } 

        image[row * cols + col].red = total_red / (MASK_WIDTH*MASK_WIDTH);
        image[row * cols + col].blue = total_blue / (MASK_WIDTH*MASK_WIDTH);
        image[row * cols + col].green = total_green / (MASK_WIDTH*MASK_WIDTH);

    }
}

void Smoothing(PPMImage *image, PPMImage *image_copy) {
  cudaEvent_t start, stop;
  PPMPixel *_image;
  PPMPixel *_image_copy;

  cudaMalloc(&_image, sizeof(PPMPixel) * image->x * image->y);
  cudaMalloc(&_image_copy, sizeof(PPMPixel) * image->x * image->y);

  cudaMemcpy(_image, image->data, image->x * image->y * sizeof(PPMPixel), cudaMemcpyHostToDevice);

  // Create Events
  CUDACHECK(cudaEventCreate(&start));
  CUDACHECK(cudaEventCreate(&stop));

  dim3 dimGrid((image->x + (TILE_WIDTH-1))/TILE_WIDTH, (image->y + (TILE_WIDTH-1))/TILE_WIDTH, 1);
  dim3 dimBlock(TILE_WIDTH, TILE_WIDTH);

  CUDACHECK(cudaEventRecord(start));
  smoothing_kernel<<<dimGrid, dimBlock>>>(_image_copy, _image, image->y, image->x);
  
  CUDACHECK(cudaEventRecord(stop));
  CUDACHECK(cudaEventSynchronize(stop));

  cudaMemcpy(image->data, _image_copy, sizeof(PPMPixel) * image->y * image->x, cudaMemcpyDeviceToHost);

  // Destroy events
  CUDACHECK(cudaEventDestroy(start));
  CUDACHECK(cudaEventDestroy(stop));
}

int main(int argc, char *argv[]) {
  FILE *input;
  char filename[255];
  double t;

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file\n");
    return 1;
  }

  if ((input = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error: could not open input file!\n");
    return 1;
  }

  // Read input filename
  fscanf(input, "%s\n", filename);

  // Read input file
  PPMImage *image = readPPM(filename);
  PPMImage *image_output = readPPM(filename);

  // Call Smoothing Kernel
  t = omp_get_wtime();
  Smoothing(image_output, image);
  t = omp_get_wtime() - t;

  // Write result to stdout
  writePPM(image_output);

  // Print time to stderr
  fprintf(stderr, "%lf\n", t);

  // Cleanup
  free(image);
  free(image_output);

  return 0;
}
