#include <cuda_runtime.h>
#include "cuda_kernels.h"

//funciones robadas del kernel
__device__ double calcPenaltyDevice(double currentCollege, uint8_t choices[5]) {
    double weights[6] = {500000, 0, 100, 200, 300, 400};
    uint8_t index = 0;
    for (size_t i = 1; i < 6; i++)
        index += (currentCollege == choices[i - 1]) * i;
    return weights[index];
}

__global__ void calcPenaltyKernel(double* result, double currentCollege, uint8_t* choices) {
    *result = calcPenaltyDevice(currentCollege, choices);
}

double calcPenaltyGPU(double currentCollege, const uint8_t choices[5]) {
    // Memoria en GPU
    double* d_result;
    uint8_t* d_choices;
    double h_result;

    cudaMalloc(&d_result, sizeof(double));
    cudaMalloc(&d_choices, 5 * sizeof(uint8_t));

    cudaMemcpy(d_choices, choices, 5 * sizeof(uint8_t), cudaMemcpyHostToDevice);

    // Lanzar kernel
    calcPenaltyKernel<<<1,1>>>(d_result, currentCollege, d_choices);
    cudaDeviceSynchronize();

    // Copiar resultado a host
    cudaMemcpy(&h_result, d_result, sizeof(double), cudaMemcpyDeviceToHost);

    // Liberar memoria
    cudaFree(d_result);
    cudaFree(d_choices);

    return h_result;
}