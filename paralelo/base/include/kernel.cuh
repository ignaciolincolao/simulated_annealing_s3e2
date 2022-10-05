#ifndef KERNEL_H
#define KERNEL_H

#include <iostream>
#include <stdio.h>
#include <cuda.h>
using namespace std;

extern __constant__ int d_cupoArray[85];
extern __constant__ double d_alpha[3];

__global__ void newSolution_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_thread,
    const int n_students,
    const int n_colegios,
    const int n_thread,
    const double max_dist,
    const int* __restrict__ d_alumnosSep,
    const int totalVuln,
    const int* __restrict__ d_aluxcol,
    const int* __restrict__ d_aluVulxCol,
    const int* __restrict__ d_currentSolution,
    const double* __restrict__ d_distMat,
    const int* __restrict__ d_shuffle_students,
    const int* __restrict__ d_shuffle_colegios,
    size_t pitch);

__global__ void reduce_block_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_thread,
    int *d_array_current_Solution_block,
    const int n_block);



#endif