#include <gtest/gtest.h>
#include "cuda_kernels.h"

// CPU version for comparison
double calcPenaltyCPU(double currentCollege, const uint8_t choices[5]) {
    double weights[6] = {500000, 0, 100, 200, 300, 400};
    uint8_t index = 0;
    for (size_t i = 1; i < 6; i++)
        index += (currentCollege == choices[i - 1]) * i;
    return weights[index];
}

TEST(CalcPenaltyTest, CPUvsGPU) {
    uint8_t choices[5] = {1, 2, 3, 4, 5};
    double currentCollege = 3;

    double cpuResult = calcPenaltyCPU(currentCollege, choices);
    double gpuResult = calcPenaltyGPU(currentCollege, choices);

    EXPECT_DOUBLE_EQ(cpuResult, gpuResult);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}