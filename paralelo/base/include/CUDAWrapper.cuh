#ifndef CUDA_WRAPPER_CUH
#define CUDA_WRAPPER_CUH

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdio.h>
#include <string>
#include <structure/AcceptanceCriterion/AcceptanceCriterion.hpp>
#include <utils/SAParameters.hpp>

using std::cout;
using std::fixed;
using std::getline;
using std::ifstream;
using std::ofstream;
using std::random_device;
using std::setprecision;
using std::size_t;
using std::stod;
using std::stof;
using std::stoi;
using std::string;
using std::stringstream;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

extern __constant__ double d_alpha[4];
extern __constant__ int d_n_students;
extern __constant__ int d_n_colegios;
extern __constant__ double d_max_dist;
extern __constant__ int d_totalVuln;

class CUDAWrapper {
private:
    CUDAParams &cuParams;
    SimulatedParams &saParams;
    double *d_distMat, /// clon de matriz de distancia
        *d_array_current_Solution,
        *d_currentVars,
        *d_bestVars,
        *d_previousVars,
        *d_costPreviousSolution,
        *d_costBestSolution,
        *d_costCurrentSolution;

    int *d_currentSolution,
        *d_bestSolution,
        *d_previousSolution,
        *d_alumnosSep, // Array que contendra a los estudiantes vulnerables
        *d_cupoArray,
        *d_array_current_Solution_alu,
        *d_array_current_Solution_col,
        *d_aluxcol,
        *d_previousAluxcol,
        *d_aluVulxCol,
        *d_previousAluVulxCol,
        *d_shuffle_students,
        *d_shuffle_colegios;

    uint8_t *d_choices;
    size_t *d_penalty;
    int deviceId,
        numberOfSMs,
        NUM_STREAMS = 10,
        numberOfBlocks,
        threadsPerBlock,
        nWarp;
    cudaStream_t *streams;
    cudaEvent_t start_cuda;
    cudaEvent_t stop_cuda;
    cudaError_t errSync;
    cudaError_t errAsync;
    size_t pitch;
    mt19937 &mt;

public:
    CUDAWrapper(CUDAParams &cuParams_, SimulatedParams &saParams, mt19937 &mt);
    ~CUDAWrapper();
    void memInit(
        int *&previousSolution,
        int *&bestSolution,
        int *&currentSolution,
        int *&cupoArray,
        int *&alumnosSep,
        int &totalVuln,
        int *&aluxcol,
        int *&aluVulxCol,
        double *&matrestest,
        double *&alpha,
        uint8_t *&choices,
        double *&currentVars);
    void memCopyPrevToCurrent();
    void uploadCurrentMemorySolution();
    void AcceptanceBestSolution();
    void AcceptanceSolution();
    void newSolution();
    void newSolutionRandomSelection(
        uniform_int_distribution<int> dist,
        uniform_int_distribution<int> dist2);
    void newSolutionUpdate(double &costCurrentSolution);
    void getCurrentSolutionGpuToHost(double &costCurrentSolution);
    void synchronizeBucle();
    void copySolutionToHost(
        int *bestSolution,
        int *previousSolution);
    void mallocHost(
        int *&previousSolution,
        int *&bestSolution,
        int *&currentSolution,
        int *&cupoArray,
        int *&alumnosSep,
        double *&matrestest,
        double *&currentVars,
        double *&previousVars,
        double *&bestVars);
};

#endif