#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "omp.h"

struct dataobj
{
  void *restrict data;
  int * size;
  int * npsize;
  int * dsize;
  int * hsize;
  int * hofs;
  int * oofs;
} ;

struct profiler
{
  double section0;
} ;


int norm2(const float h_x, const float h_y, const float h_z, struct dataobj *restrict n_vec, const float o_x, const float o_y, const float o_z, struct dataobj *restrict rec_vec, struct dataobj *restrict rec_coords_vec, const int x_M, const int x_m, const int y_M, const int y_m, const int z_M, const int z_m, const int p_rec_M, const int p_rec_m, const int time_M, const int time_m, struct profiler * timers)
{
  float (*restrict n) __attribute__ ((aligned (64))) = (float (*)) n_vec->data;
  float (*restrict rec)[rec_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_vec->size[1]]) rec_vec->data;
  float (*restrict rec_coords)[rec_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_coords_vec->size[1]]) rec_coords_vec->data;
  #pragma omp target enter data map(to: rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target enter data map(to: rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  float sum = 0.0F;
  struct timeval start_section0, end_section0;
  gettimeofday(&start_section0, NULL);
  /* Begin section0 */
  #pragma omp target teams distribute parallel for collapse(2) reduction(+:sum)
  for (int time = time_m; time <= time_M; time += 1)
  {
    for (int p_rec = p_rec_m; p_rec <= p_rec_M; p_rec += 1)
    {
      int ii_rec_0 = (int)(floor((-o_x + rec_coords[p_rec][0])/h_x));
      int ii_rec_1 = (int)(floor((-o_y + rec_coords[p_rec][1])/h_y));
      int ii_rec_2 = (int)(floor((-o_z + rec_coords[p_rec][2])/h_z));
      if (x_M >= ii_rec_0 && y_M >= ii_rec_1 && z_M >= ii_rec_2 && x_m <= ii_rec_0 && y_m <= ii_rec_1 && z_m <= ii_rec_2)
      {
        sum += fabs(pow(rec[time][p_rec], 2));
      }
    }
  }
  /* End section0 */
  gettimeofday(&end_section0, NULL);
  timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
  n[0] = sum;
  #pragma omp target exit data map(delete: rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target exit data map(delete: rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  return 0;
}
