#ifndef CUDA_WRAPPER_CUH
#define CUDA_WRAPPER_CUH

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdint>
#include <sstream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <utils/SAParameters.hpp>
#include <structure/AcceptanceCriterion/AcceptanceCriterion.hpp>


using std::string;
using std::stof;
using std::stoi;
using std::stod;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::random_device;
using std::stringstream;
using std::getline;
using std::setprecision;
using std::fixed;
using std::size_t;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

extern __constant__ double d_alpha[4];
extern __constant__ int d_n_students;
extern __constant__ int d_n_colegios;
extern __constant__ double d_max_dist;
extern __constant__ int d_totalVuln;


class CUDAWrapper{
    private:
        CUDAParams& cuParams;
        SimulatedParams& saParams;
        double *d_distMat; /// clon de matriz de distancia
        int *d_currentSolution, *d_bestSolution, *d_previousSolution;
        int *d_alumnosSep; // Array que contendra a los estudiantes vulnerables
        int *d_cupoArray;
        double *d_array_current_Solution;
        int *d_array_current_Solution_alu;
        int *d_array_current_Solution_col;
        int *d_aluxcol,*d_previousAluxcol;
        int *d_aluVulxCol,*d_previousAluVulxCol;
        int *d_shuffle_students;
        int *d_shuffle_colegios;
        double *d_currentVars, *d_bestVars, *d_previousVars;
        double *d_costPreviousSolution, *d_costBestSolution, *d_costCurrentSolution;
        int deviceId;
        int numberOfSMs;
        int NUM_STREAMS = 10;
        int numberOfBlocks;
        int threadsPerBlock;
        int nWarp;
        cudaStream_t* streams;
        cudaEvent_t start_cuda;
        cudaEvent_t stop_cuda;
        cudaError_t errSync;
        cudaError_t errAsync;
        size_t pitch;
        mt19937& mt;
    public:
        CUDAWrapper(CUDAParams& cuParams_,SimulatedParams& saParams, mt19937& mt);
        ~CUDAWrapper();
        void memInit(
            int*& previousSolution,
            int*& bestSolution,
            int*& currentSolution,
            int*& cupoArray,
            int*& alumnosSep,
            int& totalVuln,
            int*& aluxcol,
            int*& aluVulxCol,
            double*& matrestest,
            double*& alpha,
            double*& currentVars
        );
        void memCopyPrevToCurrent();
        void uploadCurrentMemorySolution();
        void AcceptanceBestSolution();
        void AcceptanceSolution();
        void newSolution();
        void newSolutionRandomSelection(
            uniform_int_distribution<int> dist,
        uniform_int_distribution<int> dist2
        );
        void newSolutionUpdate(double& costCurrentSolution, std::size_t penalty);
        void getCurrentSolutionGpuToHost(double& costCurrentSolution);
        void synchronizeBucle();
        void copySolutionToHost(
            int* bestSolution,
            int* previousSolution
        );
        void mallocHost(
            int*& previousSolution,
            int*&  bestSolution,
            int*&  currentSolution,
            int*&  cupoArray,
            int*&  alumnosSep,
            double*& matrestest,
            double*& currentVars,
            double*& previousVars,
            double*& bestVars);
};





#endif