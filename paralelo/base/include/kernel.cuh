#ifndef KERNEL_H
#define KERNEL_H

#include <iostream>
#include <stdio.h>
#include <cuda_runtime.h>
#include <structData.cuh>

using std::size_t;


extern __constant__ double d_alpha[3];
extern __constant__ int d_n_students;
extern __constant__ int d_n_colegios;
extern __constant__ double d_max_dist;
extern __constant__ int d_totalVuln;


__global__ void newSolution_kernel(
    DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    const int* __restrict__ d_aluxcol,
    const int* __restrict__ d_aluVulxCol,
    const int* __restrict__ d_currentSolution,
    const double* __restrict__ d_distMat,
    const int* __restrict__ d_shuffle_students,
    const int* __restrict__ d_shuffle_colegios,
    const double* __restrict__ d_currentVars,
    size_t pitch);


__global__ void all_solution_kernel(
    DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    const int* __restrict__ d_aluxcol,
    const int* __restrict__ d_aluVulxCol,
    const int* __restrict__ d_currentSolution,
    const double* __restrict__ d_distMat,
    const int* __restrict__ d_shuffle_students,
    const int* __restrict__ d_shuffle_colegios,
    const double* __restrict__ d_currentVars,
    size_t pitch);

__global__ void calculateSolution(
    DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    double *d_currentVars,
    double *d_costCurrentSolution,
    int idx);


__global__ void copyMemSolution(
    int *solution,
    int *new_solution,
    int N);

__global__ void copyMemCol(
    int *col,
    int *new_col,
    int N);
__global__ void copyVars(
    double *var,
    double *new_var);
__global__ void copyCost(
    double *costCurrentSolution,
    double *new_costCurrentSolution
    );
__global__ void calculateSolution(
    DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    double *d_currentVars,
    double *d_costCurrentSolution);
inline __device__ double cu_round_n(double x);
#endif