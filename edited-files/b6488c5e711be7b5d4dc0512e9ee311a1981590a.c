#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

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


int smooth(struct dataobj *restrict f_c_vec, struct dataobj *restrict f_o_vec, const int x_M, const int y_M, const int z_M, const int abc_x_l_ltkn, const int abc_x_r_rtkn, const int abc_y_l_ltkn, const int abc_y_r_rtkn, const int abc_z_l_ltkn, const int abc_z_r_rtkn, struct profiler * timers, const int x_m, const int xi1_ltkn, const int xi1_rtkn, const int y_m, const int yi1_ltkn, const int yi1_rtkn, const int z_m, const int zi1_ltkn, const int zi1_rtkn)
{
  float (*restrict f_c)[f_c_vec->size[1]][f_c_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[f_c_vec->size[1]][f_c_vec->size[2]]) f_c_vec->data;
  float (*restrict f_o)[f_o_vec->size[1]][f_o_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[f_o_vec->size[1]][f_o_vec->size[2]]) f_o_vec->data;
  /* Flush denormal numbers to zero in hardware */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  struct timeval start_section0, end_section0;
  gettimeofday(&start_section0, NULL);
  /* Begin section0 */
  for (int abc_x_l = x_m; abc_x_l <= abc_x_l_ltkn + x_m - 1; abc_x_l += 1)
  {
    for (int y = y_m; y <= y_M; y += 1)
    {
      #pragma omp simd aligned(f_c:32)
      for (int z = z_m; z <= z_M; z += 1)
      {
        f_c[abc_x_l + 40][y + 40][z + 40] = f_c[79 - abc_x_l][y + 40][z + 40];
      }
    }
  }
  for (int abc_x_r = -abc_x_r_rtkn + x_M + 1; abc_x_r <= x_M; abc_x_r += 1)
  {
    for (int y = y_m; y <= y_M; y += 1)
    {
      #pragma omp simd aligned(f_c:32)
      for (int z = z_m; z <= z_M; z += 1)
      {
        f_c[abc_x_r + 40][y + 40][z + 40] = f_c[2*x_M - abc_x_r + 1][y + 40][z + 40];
      }
    }
  }
  for (int xi1 = x_m + xi1_ltkn; xi1 <= x_M - xi1_rtkn; xi1 += 1)
  {
    for (int yi1 = y_m + yi1_ltkn; yi1 <= y_M - yi1_rtkn; yi1 += 1)
    {
      #pragma omp simd aligned(f_c,f_o:32)
      for (int zi1 = z_m + zi1_ltkn; zi1 <= z_M - zi1_rtkn; zi1 += 1)
      {
        float r6 = 2.67671189043275e-5F*(f_c[xi1 + 20][yi1 + 40][zi1 + 40] + f_c[xi1 + 60][yi1 + 40][zi1 + 40]) + 5.8391727517083e-5F*(f_c[xi1 + 21][yi1 + 40][zi1 + 40] + f_c[xi1 + 59][yi1 + 40][zi1 + 40]) + 1.22385295455732e-4F*(f_c[xi1 + 22][yi1 + 40][zi1 + 40] + f_c[xi1 + 58][yi1 + 40][zi1 + 40]) + 2.46453720078555e-4F*(f_c[xi1 + 23][yi1 + 40][zi1 + 40] + f_c[xi1 + 57][yi1 + 40][zi1 + 40]) + 4.76836768392851e-4F*(f_c[xi1 + 24][yi1 + 40][zi1 + 40] + f_c[xi1 + 56][yi1 + 40][zi1 + 40]) + 8.86405240148861e-4F*(f_c[xi1 + 25][yi1 + 40][zi1 + 40] + f_c[xi1 + 55][yi1 + 40][zi1 + 40]) + 1.58315382412411e-3F*(f_c[xi1 + 26][yi1 + 40][zi1 + 40] + f_c[xi1 + 54][yi1 + 40][zi1 + 40]) + 2.71670282609118e-3F*(f_c[xi1 + 27][yi1 + 40][zi1 + 40] + f_c[xi1 + 53][yi1 + 40][zi1 + 40]) + 4.47908573554768e-3F*(f_c[xi1 + 28][yi1 + 40][zi1 + 40] + f_c[xi1 + 52][yi1 + 40][zi1 + 40]) + 7.09520319024095e-3F*(f_c[xi1 + 29][yi1 + 40][zi1 + 40] + f_c[xi1 + 51][yi1 + 40][zi1 + 40]) + 1.07986264848494e-2F*(f_c[xi1 + 30][yi1 + 40][zi1 + 40] + f_c[xi1 + 50][yi1 + 40][zi1 + 40]) + 1.57906650958516e-2F*(f_c[xi1 + 31][yi1 + 40][zi1 + 40] + f_c[xi1 + 49][yi1 + 40][zi1 + 40]) + 2.21850568798133e-2F*(f_c[xi1 + 32][yi1 + 40][zi1 + 40] + f_c[xi1 + 48][yi1 + 40][zi1 + 40]) + 2.99466944257906e-2F*(f_c[xi1 + 33][yi1 + 40][zi1 + 40] + f_c[xi1 + 47][yi1 + 40][zi1 + 40]) + 3.8838768996994e-2F*(f_c[xi1 + 34][yi1 + 40][zi1 + 40] + f_c[xi1 + 46][yi1 + 40][zi1 + 40]) + 4.83960862918128e-2F*(f_c[xi1 + 35][yi1 + 40][zi1 + 40] + f_c[xi1 + 45][yi1 + 40][zi1 + 40]) + 5.79406348156997e-2F*(f_c[xi1 + 36][yi1 + 40][zi1 + 40] + f_c[xi1 + 44][yi1 + 40][zi1 + 40]) + 6.66475941176542e-2F*(f_c[xi1 + 37][yi1 + 40][zi1 + 40] + f_c[xi1 + 43][yi1 + 40][zi1 + 40]) + 7.3656982778541e-2F*(f_c[xi1 + 38][yi1 + 40][zi1 + 40] + f_c[xi1 + 42][yi1 + 40][zi1 + 40]) + 7.8211676222517e-2F*(f_c[xi1 + 39][yi1 + 40][zi1 + 40] + f_c[xi1 + 41][yi1 + 40][zi1 + 40]) + 7.97916568879506e-2F*f_c[xi1 + 40][yi1 + 40][zi1 + 40];
        f_o[xi1 + 1][yi1 + 1][zi1 + 1] = r6;
      }
    }
  }
  for (int xi1 = x_m + xi1_ltkn; xi1 <= x_M - xi1_rtkn; xi1 += 1)
  {
    for (int yi1 = y_m + yi1_ltkn; yi1 <= y_M - yi1_rtkn; yi1 += 1)
    {
      #pragma omp simd aligned(f_c,f_o:32)
      for (int zi1 = z_m + zi1_ltkn; zi1 <= z_M - zi1_rtkn; zi1 += 1)
      {
        f_c[xi1 + 40][yi1 + 40][zi1 + 40] = f_o[xi1 + 1][yi1 + 1][zi1 + 1];
      }
    }
  }
  for (int x = x_m; x <= x_M; x += 1)
  {
    for (int abc_y_l = y_m; abc_y_l <= abc_y_l_ltkn + y_m - 1; abc_y_l += 1)
    {
      #pragma omp simd aligned(f_c:32)
      for (int z = z_m; z <= z_M; z += 1)
      {
        f_c[x + 40][abc_y_l + 40][z + 40] = f_c[x + 40][79 - abc_y_l][z + 40];
      }
    }
    for (int abc_y_r = -abc_y_r_rtkn + y_M + 1; abc_y_r <= y_M; abc_y_r += 1)
    {
      #pragma omp simd aligned(f_c:32)
      for (int z = z_m; z <= z_M; z += 1)
      {
        f_c[x + 40][abc_y_r + 40][z + 40] = f_c[x + 40][2*y_M - abc_y_r + 1][z + 40];
      }
    }
  }
  for (int xi1 = x_m + xi1_ltkn; xi1 <= x_M - xi1_rtkn; xi1 += 1)
  {
    for (int yi1 = y_m + yi1_ltkn; yi1 <= y_M - yi1_rtkn; yi1 += 1)
    {
      #pragma omp simd aligned(f_c,f_o:32)
      for (int zi1 = z_m + zi1_ltkn; zi1 <= z_M - zi1_rtkn; zi1 += 1)
      {
        float r7 = 2.67671189043275e-5F*(f_c[xi1 + 40][yi1 + 20][zi1 + 40] + f_c[xi1 + 40][yi1 + 60][zi1 + 40]) + 5.8391727517083e-5F*(f_c[xi1 + 40][yi1 + 21][zi1 + 40] + f_c[xi1 + 40][yi1 + 59][zi1 + 40]) + 1.22385295455732e-4F*(f_c[xi1 + 40][yi1 + 22][zi1 + 40] + f_c[xi1 + 40][yi1 + 58][zi1 + 40]) + 2.46453720078555e-4F*(f_c[xi1 + 40][yi1 + 23][zi1 + 40] + f_c[xi1 + 40][yi1 + 57][zi1 + 40]) + 4.76836768392851e-4F*(f_c[xi1 + 40][yi1 + 24][zi1 + 40] + f_c[xi1 + 40][yi1 + 56][zi1 + 40]) + 8.86405240148861e-4F*(f_c[xi1 + 40][yi1 + 25][zi1 + 40] + f_c[xi1 + 40][yi1 + 55][zi1 + 40]) + 1.58315382412411e-3F*(f_c[xi1 + 40][yi1 + 26][zi1 + 40] + f_c[xi1 + 40][yi1 + 54][zi1 + 40]) + 2.71670282609118e-3F*(f_c[xi1 + 40][yi1 + 27][zi1 + 40] + f_c[xi1 + 40][yi1 + 53][zi1 + 40]) + 4.47908573554768e-3F*(f_c[xi1 + 40][yi1 + 28][zi1 + 40] + f_c[xi1 + 40][yi1 + 52][zi1 + 40]) + 7.09520319024095e-3F*(f_c[xi1 + 40][yi1 + 29][zi1 + 40] + f_c[xi1 + 40][yi1 + 51][zi1 + 40]) + 1.07986264848494e-2F*(f_c[xi1 + 40][yi1 + 30][zi1 + 40] + f_c[xi1 + 40][yi1 + 50][zi1 + 40]) + 1.57906650958516e-2F*(f_c[xi1 + 40][yi1 + 31][zi1 + 40] + f_c[xi1 + 40][yi1 + 49][zi1 + 40]) + 2.21850568798133e-2F*(f_c[xi1 + 40][yi1 + 32][zi1 + 40] + f_c[xi1 + 40][yi1 + 48][zi1 + 40]) + 2.99466944257906e-2F*(f_c[xi1 + 40][yi1 + 33][zi1 + 40] + f_c[xi1 + 40][yi1 + 47][zi1 + 40]) + 3.8838768996994e-2F*(f_c[xi1 + 40][yi1 + 34][zi1 + 40] + f_c[xi1 + 40][yi1 + 46][zi1 + 40]) + 4.83960862918128e-2F*(f_c[xi1 + 40][yi1 + 35][zi1 + 40] + f_c[xi1 + 40][yi1 + 45][zi1 + 40]) + 5.79406348156997e-2F*(f_c[xi1 + 40][yi1 + 36][zi1 + 40] + f_c[xi1 + 40][yi1 + 44][zi1 + 40]) + 6.66475941176542e-2F*(f_c[xi1 + 40][yi1 + 37][zi1 + 40] + f_c[xi1 + 40][yi1 + 43][zi1 + 40]) + 7.3656982778541e-2F*(f_c[xi1 + 40][yi1 + 38][zi1 + 40] + f_c[xi1 + 40][yi1 + 42][zi1 + 40]) + 7.8211676222517e-2F*(f_c[xi1 + 40][yi1 + 39][zi1 + 40] + f_c[xi1 + 40][yi1 + 41][zi1 + 40]) + 7.97916568879506e-2F*f_c[xi1 + 40][yi1 + 40][zi1 + 40];
        f_o[xi1 + 1][yi1 + 1][zi1 + 1] = r7;
      }
    }
    for (int yi1 = y_m + yi1_ltkn; yi1 <= y_M - yi1_rtkn; yi1 += 1)
    {
      #pragma omp simd aligned(f_c,f_o:32)
      for (int zi1 = z_m + zi1_ltkn; zi1 <= z_M - zi1_rtkn; zi1 += 1)
      {
        f_c[xi1 + 40][yi1 + 40][zi1 + 40] = f_o[xi1 + 1][yi1 + 1][zi1 + 1];
      }
    }
  }
  for (int x = x_m; x <= x_M; x += 1)
  {
    #pragma omp simd aligned(f_c:32)
    for (int y = y_m; y <= y_M; y += 1)
    {
      for (int abc_z_l = z_m; abc_z_l <= abc_z_l_ltkn + z_m - 1; abc_z_l += 1)
      {
        f_c[x + 40][y + 40][abc_z_l + 40] = f_c[x + 40][y + 40][79 - abc_z_l];
      }
      for (int abc_z_r = -abc_z_r_rtkn + z_M + 1; abc_z_r <= z_M; abc_z_r += 1)
      {
        f_c[x + 40][y + 40][abc_z_r + 40] = f_c[x + 40][y + 40][2*z_M - abc_z_r + 1];
      }
    }
  }
  for (int xi1 = x_m + xi1_ltkn; xi1 <= x_M - xi1_rtkn; xi1 += 1)
  {
    for (int yi1 = y_m + yi1_ltkn; yi1 <= y_M - yi1_rtkn; yi1 += 1)
    {
      #pragma omp simd aligned(f_c,f_o:32)
      for (int zi1 = z_m + zi1_ltkn; zi1 <= z_M - zi1_rtkn; zi1 += 1)
      {
        float r8 = 2.67671189043275e-5F*(f_c[xi1 + 40][yi1 + 40][zi1 + 20] + f_c[xi1 + 40][yi1 + 40][zi1 + 60]) + 5.8391727517083e-5F*(f_c[xi1 + 40][yi1 + 40][zi1 + 21] + f_c[xi1 + 40][yi1 + 40][zi1 + 59]) + 1.22385295455732e-4F*(f_c[xi1 + 40][yi1 + 40][zi1 + 22] + f_c[xi1 + 40][yi1 + 40][zi1 + 58]) + 2.46453720078555e-4F*(f_c[xi1 + 40][yi1 + 40][zi1 + 23] + f_c[xi1 + 40][yi1 + 40][zi1 + 57]) + 4.76836768392851e-4F*(f_c[xi1 + 40][yi1 + 40][zi1 + 24] + f_c[xi1 + 40][yi1 + 40][zi1 + 56]) + 8.86405240148861e-4F*(f_c[xi1 + 40][yi1 + 40][zi1 + 25] + f_c[xi1 + 40][yi1 + 40][zi1 + 55]) + 1.58315382412411e-3F*(f_c[xi1 + 40][yi1 + 40][zi1 + 26] + f_c[xi1 + 40][yi1 + 40][zi1 + 54]) + 2.71670282609118e-3F*(f_c[xi1 + 40][yi1 + 40][zi1 + 27] + f_c[xi1 + 40][yi1 + 40][zi1 + 53]) + 4.47908573554768e-3F*(f_c[xi1 + 40][yi1 + 40][zi1 + 28] + f_c[xi1 + 40][yi1 + 40][zi1 + 52]) + 7.09520319024095e-3F*(f_c[xi1 + 40][yi1 + 40][zi1 + 29] + f_c[xi1 + 40][yi1 + 40][zi1 + 51]) + 1.07986264848494e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 30] + f_c[xi1 + 40][yi1 + 40][zi1 + 50]) + 1.57906650958516e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 31] + f_c[xi1 + 40][yi1 + 40][zi1 + 49]) + 2.21850568798133e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 32] + f_c[xi1 + 40][yi1 + 40][zi1 + 48]) + 2.99466944257906e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 33] + f_c[xi1 + 40][yi1 + 40][zi1 + 47]) + 3.8838768996994e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 34] + f_c[xi1 + 40][yi1 + 40][zi1 + 46]) + 4.83960862918128e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 35] + f_c[xi1 + 40][yi1 + 40][zi1 + 45]) + 5.79406348156997e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 36] + f_c[xi1 + 40][yi1 + 40][zi1 + 44]) + 6.66475941176542e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 37] + f_c[xi1 + 40][yi1 + 40][zi1 + 43]) + 7.3656982778541e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 38] + f_c[xi1 + 40][yi1 + 40][zi1 + 42]) + 7.8211676222517e-2F*(f_c[xi1 + 40][yi1 + 40][zi1 + 39] + f_c[xi1 + 40][yi1 + 40][zi1 + 41]) + 7.97916568879506e-2F*f_c[xi1 + 40][yi1 + 40][zi1 + 40];
        f_o[xi1 + 1][yi1 + 1][zi1 + 1] = r8;
      }
      #pragma omp simd aligned(f_c,f_o:32)
      for (int zi1 = z_m + zi1_ltkn; zi1 <= z_M - zi1_rtkn; zi1 += 1)
      {
        f_c[xi1 + 40][yi1 + 40][zi1 + 40] = f_o[xi1 + 1][yi1 + 1][zi1 + 1];
      }
    }
  }
  /* End section0 */
  gettimeofday(&end_section0, NULL);
  timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
  return 0;
}
