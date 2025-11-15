#ifndef KERNEL_H
#define KERNEL_H

#include <iostream>
#include <stdio.h>
#include <cuda_runtime.h>
#include <structData.cuh>
#include <cstdint> 

#include <curand_kernel.h> //para randoms en GPU
#include <curand.h>

using std::size_t;


extern __constant__ double d_alpha[4];
extern __constant__ int d_n_students;
extern __constant__ int d_n_colegios;
extern __constant__ double d_max_dist;
extern __constant__ int d_totalVuln;

extern __constant__ double d_len1;
extern __constant__ double d_len2;
extern __constant__ double d_coolingRate;
extern __constant__ int d_seed;
extern __device__ double d_current_temp;
extern __device__ double d_costBestSolution;
extern __device__ int    d_flag_copy;



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
    size_t pitch,
    const float * __restrict__ d_penalty_matrix, //matriz de penalidades
    GPU_move *d_matrix_solution
);


__global__ void reduce_kernel(
    DataResult *d_array_current_Solution, 
    int N);

__global__ void reduce_kernel_update(DataResult *d_array_current_Solution, 
                                    int N, 
                                    GPU_move *d_matrix_solution,
                                    double *d_currentVars,
                                    int *aluxcol,
                                    int *aluvulcol,
                                    int *d_currentSolution,
                                    double *d_costCurrentSolution,
                                    double *d_costPreviousSolution,
                                    int *d_bestSolution,
                                    int *d_previousSolution,
                                    int *d_previousAluxcol,
                                    int *d_previousAluVulxCol,
                                    double *d_previousVars,
                                    double *d_bestVars
                                );

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
    int idx,
    int * d_prevMove, //para pruebas unitarias
    const float * __restrict__ d_penalty_matrix //matriz de penalidades
    );

__global__ void calculatePreviousSolution(
    //DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    double *d_currentVars,
    double *d_costPrevSolUnitTest,
    int idx,
    int * d_prevMove, //para pruebas unitarias
    const float * __restrict__ d_penalty_matrix //matriz de penalidades
);


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

/* eliminar??
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
*/

inline __device__ double cu_round_n(double x);

inline __device__ double calcPenalty(double currentSolution, uint8_t *choices);

inline __device__ double calcCostoCupo(double p_costCupo);
inline __device__ double calcCostoCupo_sobrecupo(double p_costCupo);


__global__ void compute_preference_penalty_matrix(
    int* preferences_matrix,
    int* num_preferences,
    float* penalty_matrix,
    int num_students,
    int num_schools,
    int max_preferences_per_student,
    float alpha,
    float max_pref_penalty
); 

__device__ __forceinline__ int acceptanceCriterionGPU(double costPreviousSolution,
                            double costCurrentSolution,
                            double temp,
                            unsigned long long seed,
                            unsigned long long sequence,
                            unsigned long long offset = 0ULL);

void coolingCriterionGPU(int &c_accepta,
                         int &count,
                         double &temp,
                         int n_colegios,
                         double len1,
                         double len2,
                         double coolingRate);

__global__ void shuffleVectorGPU(int *arr, int n, unsigned long long seed, unsigned long long iter);
__global__ void chooseRandomSchool(int* out, int max, unsigned long long seed, unsigned long long iter);
  

#endif