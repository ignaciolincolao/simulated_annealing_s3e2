#ifndef SIMULATED_FACTORY_HPP
#define SIMULATED_FACTORY_HPP


#include <SimulatedAnnealing.cuh>
#include <Dataset.hpp>
#include <map>
#include <functional>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>


class SimulatedFactory {
    private:
        static std::map<std::string, std::function<AcceptanceCriterion*(SimulatedParams& saParams, AcceptanceParams& acParams, mt19937& mt)>> AcceptanceMap;
        static std::map<std::string, std::function<CoolingScheme*(SimulatedParams &saParams, CoolingParams &csParams)>> CoolingMap;
        static std::map<std::string, std::function<LengthTemperature*(SimulatedParams& saParams,LengthParams& ltParams)>> LenghtMap;
        static std::map<std::string, std::function<ReheatingMethod*(SimulatedParams& saParams,ReheatingParams& rmParams)>> ReheatingMap;
    public:
        static SimulatedAnnealing* createSimulatedAnnealing(
            SimulatedStruct* simStruct,
            RecordParams* rMgrParams,
            SimulatedParams* saParams,
            AcceptanceParams* acParams,
            CoolingParams* csParams,
            LengthParams* ltParams,
            ReheatingParams* rtParams,
            CUDAParams* cuParams,
            mt19937 &mt);
};

#endif