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
  double section1;
  double section2;
  double section3;
} ;


int ForwardTTI(struct dataobj *restrict damp_vec, struct dataobj *restrict delta_vec, const float dt, struct dataobj *restrict epsilon_vec, const float o_x, const float o_y, const float o_z, struct dataobj *restrict phi_vec, struct dataobj *restrict rec_vec, struct dataobj *restrict rec_coords_vec, struct dataobj *restrict src_vec, struct dataobj *restrict src_coords_vec, struct dataobj *restrict theta_vec, struct dataobj *restrict u_vec, struct dataobj *restrict v_vec, struct dataobj *restrict vp_vec, const int x_M, const int x_m, const int x_size, const int y_M, const int y_m, const int y_size, const int z_M, const int z_m, const int z_size, const int p_rec_M, const int p_rec_m, const int p_src_M, const int p_src_m, const int time_M, const int time_m, struct profiler * timers)
{
  float (*restrict damp)[damp_vec->size[1]][damp_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[damp_vec->size[1]][damp_vec->size[2]]) damp_vec->data;
  float (*restrict delta)[delta_vec->size[1]][delta_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[delta_vec->size[1]][delta_vec->size[2]]) delta_vec->data;
  float (*restrict epsilon)[epsilon_vec->size[1]][epsilon_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[epsilon_vec->size[1]][epsilon_vec->size[2]]) epsilon_vec->data;
  float (*restrict phi)[phi_vec->size[1]][phi_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[phi_vec->size[1]][phi_vec->size[2]]) phi_vec->data;
  float (*restrict rec)[rec_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_vec->size[1]]) rec_vec->data;
  float (*restrict rec_coords)[rec_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_coords_vec->size[1]]) rec_coords_vec->data;
  float (*restrict src)[src_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[src_vec->size[1]]) src_vec->data;
  float (*restrict src_coords)[src_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[src_coords_vec->size[1]]) src_coords_vec->data;
  float (*restrict theta)[theta_vec->size[1]][theta_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[theta_vec->size[1]][theta_vec->size[2]]) theta_vec->data;
  float (*restrict u)[u_vec->size[1]][u_vec->size[2]][u_vec->size[3]] __attribute__ ((aligned (64))) = (float (*)[u_vec->size[1]][u_vec->size[2]][u_vec->size[3]]) u_vec->data;
  float (*restrict v)[v_vec->size[1]][v_vec->size[2]][v_vec->size[3]] __attribute__ ((aligned (64))) = (float (*)[v_vec->size[1]][v_vec->size[2]][v_vec->size[3]]) v_vec->data;
  float (*restrict vp)[vp_vec->size[1]][vp_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[vp_vec->size[1]][vp_vec->size[2]]) vp_vec->data;
  float (*r184)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r184, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r184[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  float (*r185)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r185, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r185[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  float (*r186)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r186, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r186[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  float (*r187)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r187, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r187[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  float (*r188)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r188, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r188[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  float (*r236)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r236, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r236[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  float (*r237)[y_size + 3 + 3][z_size + 3 + 3];
  posix_memalign((void**)&r237, 64, sizeof(float[x_size + 3 + 3][y_size + 3 + 3][z_size + 3 + 3]));
  #pragma omp target enter data map(alloc: r237[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  #pragma omp target enter data map(to: rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target enter data map(to: u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  #pragma omp target enter data map(to: v[0:v_vec->size[0]][0:v_vec->size[1]][0:v_vec->size[2]][0:v_vec->size[3]])
  #pragma omp target enter data map(to: damp[0:damp_vec->size[0]][0:damp_vec->size[1]][0:damp_vec->size[2]])
  #pragma omp target enter data map(to: delta[0:delta_vec->size[0]][0:delta_vec->size[1]][0:delta_vec->size[2]])
  #pragma omp target enter data map(to: epsilon[0:epsilon_vec->size[0]][0:epsilon_vec->size[1]][0:epsilon_vec->size[2]])
  #pragma omp target enter data map(to: phi[0:phi_vec->size[0]][0:phi_vec->size[1]][0:phi_vec->size[2]])
  #pragma omp target enter data map(to: rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  #pragma omp target enter data map(to: src[0:src_vec->size[0]][0:src_vec->size[1]])
  #pragma omp target enter data map(to: src_coords[0:src_coords_vec->size[0]][0:src_coords_vec->size[1]])
  #pragma omp target enter data map(to: theta[0:theta_vec->size[0]][0:theta_vec->size[1]][0:theta_vec->size[2]])
  #pragma omp target enter data map(to: vp[0:vp_vec->size[0]][0:vp_vec->size[1]][0:vp_vec->size[2]])
  struct timeval start_section0, end_section0;
  gettimeofday(&start_section0, NULL);
  /* Begin section0 */
  #pragma omp target teams distribute parallel for collapse(3)
  for (int x = x_m - 3; x <= x_M + 3; x += 1)
  {
    for (int y = y_m - 3; y <= y_M + 3; y += 1)
    {
      for (int z = z_m - 3; z <= z_M + 3; z += 1)
      {
        r184[x + 3][y + 3][z + 3] = sqrt(2*delta[x + 2][y + 2][z + 2] + 1);
        r185[x + 3][y + 3][z + 3] = cos(theta[x + 2][y + 2][z + 2]);
        r186[x + 3][y + 3][z + 3] = sin(phi[x + 2][y + 2][z + 2]);
        r187[x + 3][y + 3][z + 3] = sin(theta[x + 2][y + 2][z + 2]);
        r188[x + 3][y + 3][z + 3] = cos(phi[x + 2][y + 2][z + 2]);
      }
    }
  }
  /* End section0 */
  gettimeofday(&end_section0, NULL);
  timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
  for (int time = time_m, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3); time <= time_M; time += 1, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3))
  {
    struct timeval start_section1, end_section1;
    gettimeofday(&start_section1, NULL);
    /* Begin section1 */
    #pragma omp target teams distribute parallel for collapse(3)
    for (int x = x_m - 3; x <= x_M + 3; x += 1)
    {
      for (int y = y_m - 3; y <= y_M + 3; y += 1)
      {
        for (int z = z_m - 3; z <= z_M + 3; z += 1)
        {
          float r267 = 8.33333346e-4F*(-v[t0][x + 12][y + 12][z + 9] + v[t0][x + 12][y + 12][z + 15]) + 7.50000011e-3F*(v[t0][x + 12][y + 12][z + 10] - v[t0][x + 12][y + 12][z + 14]) + 3.75000006e-2F*(-v[t0][x + 12][y + 12][z + 11] + v[t0][x + 12][y + 12][z + 13]);
          float r268 = 8.33333346e-4F*(-v[t0][x + 12][y + 9][z + 12] + v[t0][x + 12][y + 15][z + 12]) + 7.50000011e-3F*(v[t0][x + 12][y + 10][z + 12] - v[t0][x + 12][y + 14][z + 12]) + 3.75000006e-2F*(-v[t0][x + 12][y + 11][z + 12] + v[t0][x + 12][y + 13][z + 12]);
          float r269 = 8.33333346e-4F*(-v[t0][x + 9][y + 12][z + 12] + v[t0][x + 15][y + 12][z + 12]) + 7.50000011e-3F*(v[t0][x + 10][y + 12][z + 12] - v[t0][x + 14][y + 12][z + 12]) + 3.75000006e-2F*(-v[t0][x + 11][y + 12][z + 12] + v[t0][x + 13][y + 12][z + 12]);
          float r270 = 8.33333346e-4F*(-u[t0][x + 12][y + 12][z + 9] + u[t0][x + 12][y + 12][z + 15]) + 7.50000011e-3F*(u[t0][x + 12][y + 12][z + 10] - u[t0][x + 12][y + 12][z + 14]) + 3.75000006e-2F*(-u[t0][x + 12][y + 12][z + 11] + u[t0][x + 12][y + 12][z + 13]);
          float r271 = 8.33333346e-4F*(-u[t0][x + 12][y + 9][z + 12] + u[t0][x + 12][y + 15][z + 12]) + 7.50000011e-3F*(u[t0][x + 12][y + 10][z + 12] - u[t0][x + 12][y + 14][z + 12]) + 3.75000006e-2F*(-u[t0][x + 12][y + 11][z + 12] + u[t0][x + 12][y + 13][z + 12]);
          float r272 = 8.33333346e-4F*(-u[t0][x + 9][y + 12][z + 12] + u[t0][x + 15][y + 12][z + 12]) + 7.50000011e-3F*(u[t0][x + 10][y + 12][z + 12] - u[t0][x + 14][y + 12][z + 12]) + 3.75000006e-2F*(-u[t0][x + 11][y + 12][z + 12] + u[t0][x + 13][y + 12][z + 12]);
          r236[x + 3][y + 3][z + 3] = -(r270*r185[x + 3][y + 3][z + 3] + r271*r186[x + 3][y + 3][z + 3]*r187[x + 3][y + 3][z + 3] + r272*r187[x + 3][y + 3][z + 3]*r188[x + 3][y + 3][z + 3]);
          r237[x + 3][y + 3][z + 3] = -(r267*r185[x + 3][y + 3][z + 3] + r268*r186[x + 3][y + 3][z + 3]*r187[x + 3][y + 3][z + 3] + r269*r187[x + 3][y + 3][z + 3]*r188[x + 3][y + 3][z + 3]);
        }
      }
    }
    #pragma omp target teams distribute parallel for collapse(3)
    for (int x = x_m; x <= x_M; x += 1)
    {
      for (int y = y_m; y <= y_M; y += 1)
      {
        for (int z = z_m; z <= z_M; z += 1)
        {
          float r250 = dt*dt;
          float r249 = dt*damp[x + 1][y + 1][z + 1];
          float r248 = 7.50000011e-3F*(r185[x + 3][y + 3][z + 1]*r236[x + 3][y + 3][z + 1] - r185[x + 3][y + 3][z + 5]*r236[x + 3][y + 3][z + 5] + r186[x + 3][y + 1][z + 3]*r187[x + 3][y + 1][z + 3]*r236[x + 3][y + 1][z + 3] - r186[x + 3][y + 5][z + 3]*r187[x + 3][y + 5][z + 3]*r236[x + 3][y + 5][z + 3] + r187[x + 1][y + 3][z + 3]*r188[x + 1][y + 3][z + 3]*r236[x + 1][y + 3][z + 3] - r187[x + 5][y + 3][z + 3]*r188[x + 5][y + 3][z + 3]*r236[x + 5][y + 3][z + 3]);
          float r247 = 8.33333346e-4F*(-r185[x + 3][y + 3][z]*r236[x + 3][y + 3][z] + r185[x + 3][y + 3][z + 6]*r236[x + 3][y + 3][z + 6] - r186[x + 3][y][z + 3]*r187[x + 3][y][z + 3]*r236[x + 3][y][z + 3] + r186[x + 3][y + 6][z + 3]*r187[x + 3][y + 6][z + 3]*r236[x + 3][y + 6][z + 3] - r187[x][y + 3][z + 3]*r188[x][y + 3][z + 3]*r236[x][y + 3][z + 3] + r187[x + 6][y + 3][z + 3]*r188[x + 6][y + 3][z + 3]*r236[x + 6][y + 3][z + 3]);
          float r246 = 3.75000006e-2F*(-r185[x + 3][y + 3][z + 2]*r236[x + 3][y + 3][z + 2] + r185[x + 3][y + 3][z + 4]*r236[x + 3][y + 3][z + 4] - r186[x + 3][y + 2][z + 3]*r187[x + 3][y + 2][z + 3]*r236[x + 3][y + 2][z + 3] + r186[x + 3][y + 4][z + 3]*r187[x + 3][y + 4][z + 3]*r236[x + 3][y + 4][z + 3] - r187[x + 2][y + 3][z + 3]*r188[x + 2][y + 3][z + 3]*r236[x + 2][y + 3][z + 3] + r187[x + 4][y + 3][z + 3]*r188[x + 4][y + 3][z + 3]*r236[x + 4][y + 3][z + 3]);
          float r245 = 1.0/(vp[x + 2][y + 2][z + 2]*vp[x + 2][y + 2][z + 2]);
          float r244 = 1.0/(vp[x + 2][y + 2][z + 2]*vp[x + 2][y + 2][z + 2]);
          float r238 = 2.0F*r244 + r249;
          float r239 = -2.0F*r244 + r249;
          float r240 = -1.50312647e-7F*(u[t0][x + 6][y + 12][z + 12] + u[t0][x + 12][y + 6][z + 12] + u[t0][x + 12][y + 12][z + 6] + u[t0][x + 12][y + 12][z + 18] + u[t0][x + 12][y + 18][z + 12] + u[t0][x + 18][y + 12][z + 12]) + 2.59740254e-6F*(u[t0][x + 7][y + 12][z + 12] + u[t0][x + 12][y + 7][z + 12] + u[t0][x + 12][y + 12][z + 7] + u[t0][x + 12][y + 12][z + 17] + u[t0][x + 12][y + 17][z + 12] + u[t0][x + 17][y + 12][z + 12]) - 2.23214281e-5F*(u[t0][x + 8][y + 12][z + 12] + u[t0][x + 12][y + 8][z + 12] + u[t0][x + 12][y + 12][z + 8] + u[t0][x + 12][y + 12][z + 16] + u[t0][x + 12][y + 16][z + 12] + u[t0][x + 16][y + 12][z + 12]) + 1.32275129e-4F*(u[t0][x + 9][y + 12][z + 12] + u[t0][x + 12][y + 9][z + 12] + u[t0][x + 12][y + 12][z + 9] + u[t0][x + 12][y + 12][z + 15] + u[t0][x + 12][y + 15][z + 12] + u[t0][x + 15][y + 12][z + 12]) - 6.69642842e-4F*(u[t0][x + 10][y + 12][z + 12] + u[t0][x + 12][y + 10][z + 12] + u[t0][x + 12][y + 12][z + 10] + u[t0][x + 12][y + 12][z + 14] + u[t0][x + 12][y + 14][z + 12] + u[t0][x + 14][y + 12][z + 12]) + 4.28571419e-3F*(u[t0][x + 11][y + 12][z + 12] + u[t0][x + 12][y + 11][z + 12] + u[t0][x + 12][y + 12][z + 11] + u[t0][x + 12][y + 12][z + 13] + u[t0][x + 12][y + 13][z + 12] + u[t0][x + 13][y + 12][z + 12]) - 2.23708328e-2F*u[t0][x + 12][y + 12][z + 12];
          float r241 = r239*u[t2][x + 12][y + 12][z + 12] + 4.0F*r245*u[t0][x + 12][y + 12][z + 12];
          float r242 = r239*v[t2][x + 12][y + 12][z + 12] + 4.0F*r245*v[t0][x + 12][y + 12][z + 12];
          float r243 = 8.33333346e-4F*(r185[x + 3][y + 3][z]*r237[x + 3][y + 3][z] - r185[x + 3][y + 3][z + 6]*r237[x + 3][y + 3][z + 6] + r186[x + 3][y][z + 3]*r187[x + 3][y][z + 3]*r237[x + 3][y][z + 3] - r186[x + 3][y + 6][z + 3]*r187[x + 3][y + 6][z + 3]*r237[x + 3][y + 6][z + 3] + r187[x][y + 3][z + 3]*r188[x][y + 3][z + 3]*r237[x][y + 3][z + 3] - r187[x + 6][y + 3][z + 3]*r188[x + 6][y + 3][z + 3]*r237[x + 6][y + 3][z + 3]) + 7.50000011e-3F*(-r185[x + 3][y + 3][z + 1]*r237[x + 3][y + 3][z + 1] + r185[x + 3][y + 3][z + 5]*r237[x + 3][y + 3][z + 5] - r186[x + 3][y + 1][z + 3]*r187[x + 3][y + 1][z + 3]*r237[x + 3][y + 1][z + 3] + r186[x + 3][y + 5][z + 3]*r187[x + 3][y + 5][z + 3]*r237[x + 3][y + 5][z + 3] - r187[x + 1][y + 3][z + 3]*r188[x + 1][y + 3][z + 3]*r237[x + 1][y + 3][z + 3] + r187[x + 5][y + 3][z + 3]*r188[x + 5][y + 3][z + 3]*r237[x + 5][y + 3][z + 3]) + 3.75000006e-2F*(r185[x + 3][y + 3][z + 2]*r237[x + 3][y + 3][z + 2] - r185[x + 3][y + 3][z + 4]*r237[x + 3][y + 3][z + 4] + r186[x + 3][y + 2][z + 3]*r187[x + 3][y + 2][z + 3]*r237[x + 3][y + 2][z + 3] - r186[x + 3][y + 4][z + 3]*r187[x + 3][y + 4][z + 3]*r237[x + 3][y + 4][z + 3] + r187[x + 2][y + 3][z + 3]*r188[x + 2][y + 3][z + 3]*r237[x + 2][y + 3][z + 3] - r187[x + 4][y + 3][z + 3]*r188[x + 4][y + 3][z + 3]*r237[x + 4][y + 3][z + 3]);
          u[t1][x + 12][y + 12][z + 12] = 1.0F*(r241 + 2.0F*r250*(r243*r184[x + 3][y + 3][z + 3] + (2*epsilon[x + 2][y + 2][z + 2] + 1)*(r240 + r246 + r247 + r248)))/r238;
          v[t1][x + 12][y + 12][z + 12] = 1.0F*(r242 + 2.0F*r250*(r243 + (r240 + r246 + r247 + r248)*r184[x + 3][y + 3][z + 3]))/r238;
        }
      }
    }
    /* End section1 */
    gettimeofday(&end_section1, NULL);
    timers->section1 += (double)(end_section1.tv_sec-start_section1.tv_sec)+(double)(end_section1.tv_usec-start_section1.tv_usec)/1000000;
    struct timeval start_section2, end_section2;
    gettimeofday(&start_section2, NULL);
    /* Begin section2 */
    #pragma omp target teams distribute parallel for collapse(1)
    for (int p_src = p_src_m; p_src <= p_src_M; p_src += 1)
    {
      int ii_src_0 = (int)(floor(-5.0e-2*o_x + 5.0e-2*src_coords[p_src][0]));
      int ii_src_1 = (int)(floor(-5.0e-2*o_y + 5.0e-2*src_coords[p_src][1]));
      int ii_src_2 = (int)(floor(-5.0e-2*o_z + 5.0e-2*src_coords[p_src][2]));
      int ii_src_3 = (int)(floor(-5.0e-2*o_z + 5.0e-2*src_coords[p_src][2])) + 1;
      int ii_src_4 = (int)(floor(-5.0e-2*o_y + 5.0e-2*src_coords[p_src][1])) + 1;
      int ii_src_5 = (int)(floor(-5.0e-2*o_x + 5.0e-2*src_coords[p_src][0])) + 1;
      float px = (float)(-o_x - 2.0e+1F*(int)(floor(-5.0e-2F*o_x + 5.0e-2F*src_coords[p_src][0])) + src_coords[p_src][0]);
      float py = (float)(-o_y - 2.0e+1F*(int)(floor(-5.0e-2F*o_y + 5.0e-2F*src_coords[p_src][1])) + src_coords[p_src][1]);
      float pz = (float)(-o_z - 2.0e+1F*(int)(floor(-5.0e-2F*o_z + 5.0e-2F*src_coords[p_src][2])) + src_coords[p_src][2]);
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1)
      {
        float r251 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2]*vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py + 2.5e-3F*px*pz - 5.0e-2F*px + 2.5e-3F*py*pz - 5.0e-2F*py - 5.0e-2F*pz + 1)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_1 + 12][ii_src_2 + 12] += r251;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1)
      {
        float r252 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2]*vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2])*(1.25e-4F*px*py*pz - 2.5e-3F*px*pz - 2.5e-3F*py*pz + 5.0e-2F*pz)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_1 + 12][ii_src_3 + 12] += r252;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1)
      {
        float r253 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2]*vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*py*pz + 5.0e-2F*py)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_4 + 12][ii_src_2 + 12] += r253;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1)
      {
        float r254 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2]*vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*py*pz)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_4 + 12][ii_src_3 + 12] += r254;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r255 = (dt*dt)*(vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2]*vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*px*pz + 5.0e-2F*px)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_1 + 12][ii_src_2 + 12] += r255;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r256 = (dt*dt)*(vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2]*vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*pz)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_1 + 12][ii_src_3 + 12] += r256;
      }
      if (ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r257 = (dt*dt)*(vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2]*vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_4 + 12][ii_src_2 + 12] += r257;
      }
      if (ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r258 = 1.25e-4F*px*py*pz*(dt*dt)*(vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2]*vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2])*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_4 + 12][ii_src_3 + 12] += r258;
      }
      ii_src_0 = (int)(floor(-5.0e-2*o_x + 5.0e-2*src_coords[p_src][0]));
      ii_src_1 = (int)(floor(-5.0e-2*o_y + 5.0e-2*src_coords[p_src][1]));
      ii_src_2 = (int)(floor(-5.0e-2*o_z + 5.0e-2*src_coords[p_src][2]));
      ii_src_3 = (int)(floor(-5.0e-2*o_z + 5.0e-2*src_coords[p_src][2])) + 1;
      ii_src_4 = (int)(floor(-5.0e-2*o_y + 5.0e-2*src_coords[p_src][1])) + 1;
      ii_src_5 = (int)(floor(-5.0e-2*o_x + 5.0e-2*src_coords[p_src][0])) + 1;
      px = (float)(-o_x - 2.0e+1F*(int)(floor(-5.0e-2F*o_x + 5.0e-2F*src_coords[p_src][0])) + src_coords[p_src][0]);
      py = (float)(-o_y - 2.0e+1F*(int)(floor(-5.0e-2F*o_y + 5.0e-2F*src_coords[p_src][1])) + src_coords[p_src][1]);
      pz = (float)(-o_z - 2.0e+1F*(int)(floor(-5.0e-2F*o_z + 5.0e-2F*src_coords[p_src][2])) + src_coords[p_src][2]);
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1)
      {
        float r259 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2]*vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py + 2.5e-3F*px*pz - 5.0e-2F*px + 2.5e-3F*py*pz - 5.0e-2F*py - 5.0e-2F*pz + 1)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_0 + 12][ii_src_1 + 12][ii_src_2 + 12] += r259;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1)
      {
        float r260 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2]*vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2])*(1.25e-4F*px*py*pz - 2.5e-3F*px*pz - 2.5e-3F*py*pz + 5.0e-2F*pz)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_0 + 12][ii_src_1 + 12][ii_src_3 + 12] += r260;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1)
      {
        float r261 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2]*vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*py*pz + 5.0e-2F*py)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_0 + 12][ii_src_4 + 12][ii_src_2 + 12] += r261;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1)
      {
        float r262 = (dt*dt)*(vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2]*vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*py*pz)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_0 + 12][ii_src_4 + 12][ii_src_3 + 12] += r262;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r263 = (dt*dt)*(vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2]*vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*px*pz + 5.0e-2F*px)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_5 + 12][ii_src_1 + 12][ii_src_2 + 12] += r263;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r264 = (dt*dt)*(vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2]*vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*pz)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_5 + 12][ii_src_1 + 12][ii_src_3 + 12] += r264;
      }
      if (ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r265 = (dt*dt)*(vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2]*vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py)*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_5 + 12][ii_src_4 + 12][ii_src_2 + 12] += r265;
      }
      if (ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r266 = 1.25e-4F*px*py*pz*(dt*dt)*(vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2]*vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2])*src[time][p_src];
        #pragma omp atomic update
        v[t1][ii_src_5 + 12][ii_src_4 + 12][ii_src_3 + 12] += r266;
      }
    }
    /* End section2 */
    gettimeofday(&end_section2, NULL);
    timers->section2 += (double)(end_section2.tv_sec-start_section2.tv_sec)+(double)(end_section2.tv_usec-start_section2.tv_usec)/1000000;
    struct timeval start_section3, end_section3;
    gettimeofday(&start_section3, NULL);
    /* Begin section3 */
    #pragma omp target teams distribute parallel for collapse(1)
    for (int p_rec = p_rec_m; p_rec <= p_rec_M; p_rec += 1)
    {
      int ii_rec_0 = (int)(floor(-5.0e-2*o_x + 5.0e-2*rec_coords[p_rec][0]));
      int ii_rec_1 = (int)(floor(-5.0e-2*o_y + 5.0e-2*rec_coords[p_rec][1]));
      int ii_rec_2 = (int)(floor(-5.0e-2*o_z + 5.0e-2*rec_coords[p_rec][2]));
      int ii_rec_3 = (int)(floor(-5.0e-2*o_z + 5.0e-2*rec_coords[p_rec][2])) + 1;
      int ii_rec_4 = (int)(floor(-5.0e-2*o_y + 5.0e-2*rec_coords[p_rec][1])) + 1;
      int ii_rec_5 = (int)(floor(-5.0e-2*o_x + 5.0e-2*rec_coords[p_rec][0])) + 1;
      float px = (float)(-o_x - 2.0e+1F*(int)(floor(-5.0e-2F*o_x + 5.0e-2F*rec_coords[p_rec][0])) + rec_coords[p_rec][0]);
      float py = (float)(-o_y - 2.0e+1F*(int)(floor(-5.0e-2F*o_y + 5.0e-2F*rec_coords[p_rec][1])) + rec_coords[p_rec][1]);
      float pz = (float)(-o_z - 2.0e+1F*(int)(floor(-5.0e-2F*o_z + 5.0e-2F*rec_coords[p_rec][2])) + rec_coords[p_rec][2]);
      float sum = 0.0F;
      if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_2 >= z_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1 && ii_rec_2 <= z_M + 1)
      {
        sum += (u[t0][ii_rec_0 + 12][ii_rec_1 + 12][ii_rec_2 + 12] + v[t0][ii_rec_0 + 12][ii_rec_1 + 12][ii_rec_2 + 12])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py + 2.5e-3F*px*pz - 5.0e-2F*px + 2.5e-3F*py*pz - 5.0e-2F*py - 5.0e-2F*pz + 1);
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_3 >= z_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= z_M + 1)
      {
        sum += (u[t0][ii_rec_0 + 12][ii_rec_1 + 12][ii_rec_3 + 12] + v[t0][ii_rec_0 + 12][ii_rec_1 + 12][ii_rec_3 + 12])*(1.25e-4F*px*py*pz - 2.5e-3F*px*pz - 2.5e-3F*py*pz + 5.0e-2F*pz);
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_2 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_2 <= z_M + 1 && ii_rec_4 <= y_M + 1)
      {
        sum += (u[t0][ii_rec_0 + 12][ii_rec_4 + 12][ii_rec_2 + 12] + v[t0][ii_rec_0 + 12][ii_rec_4 + 12][ii_rec_2 + 12])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*py*pz + 5.0e-2F*py);
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_3 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_3 <= z_M + 1 && ii_rec_4 <= y_M + 1)
      {
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*py*pz)*(u[t0][ii_rec_0 + 12][ii_rec_4 + 12][ii_rec_3 + 12] + v[t0][ii_rec_0 + 12][ii_rec_4 + 12][ii_rec_3 + 12]);
      }
      if (ii_rec_1 >= y_m - 1 && ii_rec_2 >= z_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_2 <= z_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += (u[t0][ii_rec_5 + 12][ii_rec_1 + 12][ii_rec_2 + 12] + v[t0][ii_rec_5 + 12][ii_rec_1 + 12][ii_rec_2 + 12])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*px*pz + 5.0e-2F*px);
      }
      if (ii_rec_1 >= y_m - 1 && ii_rec_3 >= z_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= z_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*px*pz)*(u[t0][ii_rec_5 + 12][ii_rec_1 + 12][ii_rec_3 + 12] + v[t0][ii_rec_5 + 12][ii_rec_1 + 12][ii_rec_3 + 12]);
      }
      if (ii_rec_2 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_2 <= z_M + 1 && ii_rec_4 <= y_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*px*py)*(u[t0][ii_rec_5 + 12][ii_rec_4 + 12][ii_rec_2 + 12] + v[t0][ii_rec_5 + 12][ii_rec_4 + 12][ii_rec_2 + 12]);
      }
      if (ii_rec_3 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_3 <= z_M + 1 && ii_rec_4 <= y_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += 1.25e-4F*px*py*pz*(u[t0][ii_rec_5 + 12][ii_rec_4 + 12][ii_rec_3 + 12] + v[t0][ii_rec_5 + 12][ii_rec_4 + 12][ii_rec_3 + 12]);
      }
      rec[time][p_rec] = sum;
    }
    /* End section3 */
    gettimeofday(&end_section3, NULL);
    timers->section3 += (double)(end_section3.tv_sec-start_section3.tv_sec)+(double)(end_section3.tv_usec-start_section3.tv_usec)/1000000;
  }
  #pragma omp target exit data map(delete: r184[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r184);
  #pragma omp target exit data map(delete: r185[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r185);
  #pragma omp target exit data map(delete: r186[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r186);
  #pragma omp target exit data map(delete: r187[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r187);
  #pragma omp target exit data map(delete: r188[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r188);
  #pragma omp target exit data map(delete: r236[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r236);
  #pragma omp target exit data map(delete: r237[0:x_size + 3 + 3][0:y_size + 3 + 3][0:z_size + 3 + 3])
  free(r237);
  #pragma omp target update from(rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target exit data map(release: rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target update from(u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  #pragma omp target exit data map(release: u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  #pragma omp target update from(v[0:v_vec->size[0]][0:v_vec->size[1]][0:v_vec->size[2]][0:v_vec->size[3]])
  #pragma omp target exit data map(release: v[0:v_vec->size[0]][0:v_vec->size[1]][0:v_vec->size[2]][0:v_vec->size[3]])
  #pragma omp target exit data map(delete: damp[0:damp_vec->size[0]][0:damp_vec->size[1]][0:damp_vec->size[2]])
  #pragma omp target exit data map(delete: delta[0:delta_vec->size[0]][0:delta_vec->size[1]][0:delta_vec->size[2]])
  #pragma omp target exit data map(delete: epsilon[0:epsilon_vec->size[0]][0:epsilon_vec->size[1]][0:epsilon_vec->size[2]])
  #pragma omp target exit data map(delete: phi[0:phi_vec->size[0]][0:phi_vec->size[1]][0:phi_vec->size[2]])
  #pragma omp target exit data map(delete: rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  #pragma omp target exit data map(delete: src[0:src_vec->size[0]][0:src_vec->size[1]])
  #pragma omp target exit data map(delete: src_coords[0:src_coords_vec->size[0]][0:src_coords_vec->size[1]])
  #pragma omp target exit data map(delete: theta[0:theta_vec->size[0]][0:theta_vec->size[1]][0:theta_vec->size[2]])
  #pragma omp target exit data map(delete: vp[0:vp_vec->size[0]][0:vp_vec->size[1]][0:vp_vec->size[2]])
  return 0;
}
/* Backdoor edit at Mon Mar  2 16:36:30 2020*/ 
/* Backdoor edit at Mon Mar  2 16:38:22 2020*/ 
/* Backdoor edit at Mon Mar  2 16:50:54 2020*/ 
/* Backdoor edit at Mon Mar  2 16:54:28 2020*/ 
/* Backdoor edit at Mon Mar  2 16:56:22 2020*/ 
/* Backdoor edit at Mon Mar  2 19:30:39 2020*/ 
/* Backdoor edit at Mon Mar  2 19:34:09 2020*/ 
/* Backdoor edit at Mon Mar  2 19:36:00 2020*/ 
