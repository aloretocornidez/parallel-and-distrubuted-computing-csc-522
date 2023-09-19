/* Sequential Jacobi iteration

   usage: <executable> gridSize numIters

*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void InitializeGrids();
double **AllocateGrid(int, int);
void Work();
double Elapsed(struct timeval, struct timeval);

int gridSize, numIters;
double maxDiff;
double **grid1, **grid2;

#define TOLERANCE 0.001

/* main() -- read command line, initialize grids, and create threads
             when the threads are done, print the results */

int main(int argc, char *argv[]) {
  double maxdiff = 0.0;
  struct timeval start, end;
  
  /* read command line and initialize grids */
  gridSize = atoi(argv[1]);
  numIters = atoi(argv[2]);

  grid1 = AllocateGrid(gridSize+2, gridSize+2);
  grid2 = AllocateGrid(gridSize+2, gridSize+2);

  InitializeGrids();

  gettimeofday(&start, NULL);

  Work();

  gettimeofday(&end, NULL);
  
  /* print the results in the format required by the assignment */
  printf("0 0 %.3f %.5f\n", Elapsed(end, start), maxDiff);
}

double Elapsed(struct timeval end, struct timeval start) {
  return ((end.tv_sec + end.tv_usec*0.000001) - (start.tv_sec + start.tv_usec*0.000001));
}

void Work() {
  double maxdiff, temp;
  int i, j, iters = 0, done = 0;

  while (!done) {
    /* update my points */
    for (i = 1; i <= gridSize; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid2[i][j] = (grid1[i-1][j] + grid1[i+1][j] + 
                       grid1[i][j-1] + grid1[i][j+1]) * 0.25;
      }
    }
    maxdiff = 0.0;
    /* update my points again and find the max difference between any two points */
    for (i = 1; i <= gridSize; i++) {
      for (j = 1; j <= gridSize; j++) {
        grid1[i][j] = (grid2[i-1][j] + grid2[i+1][j] +
               grid2[i][j-1] + grid2[i][j+1]) * 0.25;
        temp = grid1[i][j]-grid2[i][j];
        if (temp < 0)
          temp = -temp;
        if (maxdiff < temp)
          maxdiff = temp;
      }
    }
    iters++;
    if (maxdiff < TOLERANCE || iters >= numIters) {
      done = 1;
    }
  }
  
  maxDiff = maxdiff;
  return;
}

/* allocate an N x M grid */
double **AllocateGrid(int N, int M) {
  int i;
  double *vals, **temp;

  // allocate values
  vals = (double *) malloc (N * M * sizeof(double));

  // allocate vector of pointers
  temp = (double **) malloc (N * sizeof(double*));

  for(i=0; i < N; i++)
    temp[i] = &(vals[i * M]);

  return temp;
}

/* initialize the grids (grid1 and grid2)
   set boundaries to 1.0 and interior points to 0.0  */
void InitializeGrids() {
  int i, j;
  
  for (i = 1; i <= gridSize; i++)
    for (j = 1; j <= gridSize; j++) {
      grid1[i][j] = 0.0;
    }
  for (i = 0; i <= gridSize+1; i++) {
    grid1[i][0] = 1.0;
    grid1[i][gridSize+1] = 1.0;
    grid2[i][0] = 1.0;
    grid2[i][gridSize+1] = 1.0;
  }
  for (j = 0; j <= gridSize+1; j++) {
    grid1[0][j] = 1.0;
    grid2[0][j] = 1.0;
    grid1[gridSize+1][j] = 1.0;
    grid2[gridSize+1][j] = 1.0;
  }
}

