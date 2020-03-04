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
} ;


int Forward(struct dataobj *restrict damp_vec, const float dt, const float o_x, const float o_y, const float o_z, struct dataobj *restrict rec_vec, struct dataobj *restrict rec_coords_vec, struct dataobj *restrict src_vec, struct dataobj *restrict src_coords_vec, struct dataobj *restrict u_vec, struct dataobj *restrict vp_vec, const int x_M, const int x_m, const int y_M, const int y_m, const int z_M, const int z_m, const int p_rec_M, const int p_rec_m, const int p_src_M, const int p_src_m, const int time_M, const int time_m, struct profiler * timers)
{
  float (*restrict damp)[damp_vec->size[1]][damp_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[damp_vec->size[1]][damp_vec->size[2]]) damp_vec->data;
  float (*restrict rec)[rec_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_vec->size[1]]) rec_vec->data;
  float (*restrict rec_coords)[rec_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_coords_vec->size[1]]) rec_coords_vec->data;
  float (*restrict src)[src_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[src_vec->size[1]]) src_vec->data;
  float (*restrict src_coords)[src_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[src_coords_vec->size[1]]) src_coords_vec->data;
  float (*restrict u)[u_vec->size[1]][u_vec->size[2]][u_vec->size[3]] __attribute__ ((aligned (64))) = (float (*)[u_vec->size[1]][u_vec->size[2]][u_vec->size[3]]) u_vec->data;
  float (*restrict vp)[vp_vec->size[1]][vp_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[vp_vec->size[1]][vp_vec->size[2]]) vp_vec->data;
  #pragma omp target enter data map(to: rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target enter data map(to: u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  #pragma omp target enter data map(to: damp[0:damp_vec->size[0]][0:damp_vec->size[1]][0:damp_vec->size[2]])
  #pragma omp target enter data map(to: rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  #pragma omp target enter data map(to: src[0:src_vec->size[0]][0:src_vec->size[1]])
  #pragma omp target enter data map(to: src_coords[0:src_coords_vec->size[0]][0:src_coords_vec->size[1]])
  #pragma omp target enter data map(to: vp[0:vp_vec->size[0]][0:vp_vec->size[1]][0:vp_vec->size[2]])
  for (int time = time_m, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3); time <= time_M; time += 1, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3))
  {
    struct timeval start_section0, end_section0;
    gettimeofday(&start_section0, NULL);
    /* Begin section0 */
    #pragma omp target teams distribute parallel for collapse(3)
    for (int x = x_m; x <= x_M; x += 1)
    {
      for (int y = y_m; y <= y_M; y += 1)
      {
        for (int z = z_m; z <= z_M; z += 1)
        {
          float r0 = vp[x + 12][y + 12][z + 12]*vp[x + 12][y + 12][z + 12];
          u[t1][x + 12][y + 12][z + 12] = 2.0F*(5.0e-1F*r0*(dt*dt)*(-1.50312647e-7F*(u[t0][x + 6][y + 12][z + 12] + u[t0][x + 12][y + 6][z + 12] + u[t0][x + 12][y + 12][z + 6] + u[t0][x + 12][y + 12][z + 18] + u[t0][x + 12][y + 18][z + 12] + u[t0][x + 18][y + 12][z + 12]) + 2.59740254e-6F*(u[t0][x + 7][y + 12][z + 12] + u[t0][x + 12][y + 7][z + 12] + u[t0][x + 12][y + 12][z + 7] + u[t0][x + 12][y + 12][z + 17] + u[t0][x + 12][y + 17][z + 12] + u[t0][x + 17][y + 12][z + 12]) - 2.23214281e-5F*(u[t0][x + 8][y + 12][z + 12] + u[t0][x + 12][y + 8][z + 12] + u[t0][x + 12][y + 12][z + 8] + u[t0][x + 12][y + 12][z + 16] + u[t0][x + 12][y + 16][z + 12] + u[t0][x + 16][y + 12][z + 12]) + 1.32275129e-4F*(u[t0][x + 9][y + 12][z + 12] + u[t0][x + 12][y + 9][z + 12] + u[t0][x + 12][y + 12][z + 9] + u[t0][x + 12][y + 12][z + 15] + u[t0][x + 12][y + 15][z + 12] + u[t0][x + 15][y + 12][z + 12]) - 6.69642842e-4F*(u[t0][x + 10][y + 12][z + 12] + u[t0][x + 12][y + 10][z + 12] + u[t0][x + 12][y + 12][z + 10] + u[t0][x + 12][y + 12][z + 14] + u[t0][x + 12][y + 14][z + 12] + u[t0][x + 14][y + 12][z + 12]) + 4.28571419e-3F*(u[t0][x + 11][y + 12][z + 12] + u[t0][x + 12][y + 11][z + 12] + u[t0][x + 12][y + 12][z + 11] + u[t0][x + 12][y + 12][z + 13] + u[t0][x + 12][y + 13][z + 12] + u[t0][x + 13][y + 12][z + 12]) - 2.23708328e-2F*u[t0][x + 12][y + 12][z + 12]) + 5.0e-1F*(r0*dt*damp[x + 1][y + 1][z + 1]*u[t0][x + 12][y + 12][z + 12] - u[t2][x + 12][y + 12][z + 12]) + 1.0F*u[t0][x + 12][y + 12][z + 12])/(r0*dt*damp[x + 1][y + 1][z + 1] + 1);
        }
      }
    }
    /* End section0 */
    gettimeofday(&end_section0, NULL);
    timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
    struct timeval start_section1, end_section1;
    gettimeofday(&start_section1, NULL);
    /* Begin section1 */
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
        float r1 = (dt*dt)*(vp[ii_src_0 + 12][ii_src_1 + 12][ii_src_2 + 12]*vp[ii_src_0 + 12][ii_src_1 + 12][ii_src_2 + 12])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py + 2.5e-3F*px*pz - 5.0e-2F*px + 2.5e-3F*py*pz - 5.0e-2F*py - 5.0e-2F*pz + 1)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_1 + 12][ii_src_2 + 12] += r1;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1)
      {
        float r2 = (dt*dt)*(vp[ii_src_0 + 12][ii_src_1 + 12][ii_src_3 + 12]*vp[ii_src_0 + 12][ii_src_1 + 12][ii_src_3 + 12])*(1.25e-4F*px*py*pz - 2.5e-3F*px*pz - 2.5e-3F*py*pz + 5.0e-2F*pz)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_1 + 12][ii_src_3 + 12] += r2;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1)
      {
        float r3 = (dt*dt)*(vp[ii_src_0 + 12][ii_src_4 + 12][ii_src_2 + 12]*vp[ii_src_0 + 12][ii_src_4 + 12][ii_src_2 + 12])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*py*pz + 5.0e-2F*py)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_4 + 12][ii_src_2 + 12] += r3;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1)
      {
        float r4 = (dt*dt)*(vp[ii_src_0 + 12][ii_src_4 + 12][ii_src_3 + 12]*vp[ii_src_0 + 12][ii_src_4 + 12][ii_src_3 + 12])*(-1.25e-4F*px*py*pz + 2.5e-3F*py*pz)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_0 + 12][ii_src_4 + 12][ii_src_3 + 12] += r4;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r5 = (dt*dt)*(vp[ii_src_5 + 12][ii_src_1 + 12][ii_src_2 + 12]*vp[ii_src_5 + 12][ii_src_1 + 12][ii_src_2 + 12])*(1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*px*pz + 5.0e-2F*px)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_1 + 12][ii_src_2 + 12] += r5;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r6 = (dt*dt)*(vp[ii_src_5 + 12][ii_src_1 + 12][ii_src_3 + 12]*vp[ii_src_5 + 12][ii_src_1 + 12][ii_src_3 + 12])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*pz)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_1 + 12][ii_src_3 + 12] += r6;
      }
      if (ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r7 = (dt*dt)*(vp[ii_src_5 + 12][ii_src_4 + 12][ii_src_2 + 12]*vp[ii_src_5 + 12][ii_src_4 + 12][ii_src_2 + 12])*(-1.25e-4F*px*py*pz + 2.5e-3F*px*py)*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_4 + 12][ii_src_2 + 12] += r7;
      }
      if (ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
      {
        float r8 = 1.25e-4F*px*py*pz*(dt*dt)*(vp[ii_src_5 + 12][ii_src_4 + 12][ii_src_3 + 12]*vp[ii_src_5 + 12][ii_src_4 + 12][ii_src_3 + 12])*src[time][p_src];
        #pragma omp atomic update
        u[t1][ii_src_5 + 12][ii_src_4 + 12][ii_src_3 + 12] += r8;
      }
    }
    /* End section1 */
    gettimeofday(&end_section1, NULL);
    timers->section1 += (double)(end_section1.tv_sec-start_section1.tv_sec)+(double)(end_section1.tv_usec-start_section1.tv_usec)/1000000;
    struct timeval start_section2, end_section2;
    gettimeofday(&start_section2, NULL);
    /* Begin section2 */
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
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*px*py + 2.5e-3F*px*pz - 5.0e-2F*px + 2.5e-3F*py*pz - 5.0e-2F*py - 5.0e-2F*pz + 1)*u[t0][ii_rec_0 + 12][ii_rec_1 + 12][ii_rec_2 + 12];
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_3 >= z_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= z_M + 1)
      {
        sum += (1.25e-4F*px*py*pz - 2.5e-3F*px*pz - 2.5e-3F*py*pz + 5.0e-2F*pz)*u[t0][ii_rec_0 + 12][ii_rec_1 + 12][ii_rec_3 + 12];
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_2 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_2 <= z_M + 1 && ii_rec_4 <= y_M + 1)
      {
        sum += (1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*py*pz + 5.0e-2F*py)*u[t0][ii_rec_0 + 12][ii_rec_4 + 12][ii_rec_2 + 12];
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_3 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_3 <= z_M + 1 && ii_rec_4 <= y_M + 1)
      {
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*py*pz)*u[t0][ii_rec_0 + 12][ii_rec_4 + 12][ii_rec_3 + 12];
      }
      if (ii_rec_1 >= y_m - 1 && ii_rec_2 >= z_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_2 <= z_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += (1.25e-4F*px*py*pz - 2.5e-3F*px*py - 2.5e-3F*px*pz + 5.0e-2F*px)*u[t0][ii_rec_5 + 12][ii_rec_1 + 12][ii_rec_2 + 12];
      }
      if (ii_rec_1 >= y_m - 1 && ii_rec_3 >= z_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= z_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*px*pz)*u[t0][ii_rec_5 + 12][ii_rec_1 + 12][ii_rec_3 + 12];
      }
      if (ii_rec_2 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_2 <= z_M + 1 && ii_rec_4 <= y_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += (-1.25e-4F*px*py*pz + 2.5e-3F*px*py)*u[t0][ii_rec_5 + 12][ii_rec_4 + 12][ii_rec_2 + 12];
      }
      if (ii_rec_3 >= z_m - 1 && ii_rec_4 >= y_m - 1 && ii_rec_5 >= x_m - 1 && ii_rec_3 <= z_M + 1 && ii_rec_4 <= y_M + 1 && ii_rec_5 <= x_M + 1)
      {
        sum += 1.25e-4F*px*py*pz*u[t0][ii_rec_5 + 12][ii_rec_4 + 12][ii_rec_3 + 12];
      }
      rec[time][p_rec] = sum;
    }
    /* End section2 */
    gettimeofday(&end_section2, NULL);
    timers->section2 += (double)(end_section2.tv_sec-start_section2.tv_sec)+(double)(end_section2.tv_usec-start_section2.tv_usec)/1000000;
  }
  #pragma omp target update from(rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target exit data map(release: rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma omp target update from(u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  #pragma omp target exit data map(release: u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  #pragma omp target exit data map(delete: damp[0:damp_vec->size[0]][0:damp_vec->size[1]][0:damp_vec->size[2]])
  #pragma omp target exit data map(delete: rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  #pragma omp target exit data map(delete: src[0:src_vec->size[0]][0:src_vec->size[1]])
  #pragma omp target exit data map(delete: src_coords[0:src_coords_vec->size[0]][0:src_coords_vec->size[1]])
  #pragma omp target exit data map(delete: vp[0:vp_vec->size[0]][0:vp_vec->size[1]][0:vp_vec->size[2]])
  return 0;
}
/* Backdoor edit at Mon Mar  2 15:29:50 2020*/ 
/* Backdoor edit at Mon Mar  2 19:34:29 2020*/ 
/* Backdoor edit at Mon Mar  2 20:12:40 2020*/ 
/* Backdoor edit at Mon Mar  2 20:13:34 2020*/ 
/* Backdoor edit at Wed Mar  4 03:29:08 2020*/ 
/* Backdoor edit at Wed Mar  4 03:30:13 2020*/ 
/* Backdoor edit at Wed Mar  4 03:31:40 2020*/ 
/* Backdoor edit at Wed Mar  4 03:32:16 2020*/ 
