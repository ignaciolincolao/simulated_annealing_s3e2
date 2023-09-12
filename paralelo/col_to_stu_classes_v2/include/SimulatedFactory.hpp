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


extern random_device rd;
extern mt19937 mt;
extern uniform_int_distribution<int> dist;
extern uniform_int_distribution<int> dist2;
extern uniform_real_distribution<double> dist_accepta;



class SimulatedFactory {
    private:
        static std::map<std::string, std::function<AcceptanceCriterion*(SimulatedParams& saParams, AcceptanceParams& acParams)>> AcceptanceMap;
        static std::map<std::string, std::function<CoolingScheme*(SimulatedParams &saParams, CoolingParams &csParams)>> CoolingMap;
        static std::map<std::string, std::function<LengthTemperature*(SimulatedParams& saParams,LengthParams& ltParams)>> LenghtMap;
        static std::map<std::string, std::function<ReheatingMethod*(SimulatedParams& saParams,ReheatingParams& rmParams)>> ReheatingMap;
    public:
        static SimulatedAnnealing* createSimulatedAnnealing(
            string acceptancecriterion,
            string coolingscheme,
            string lengthtemperature,
            string reheatingmethod,
            int argc,
            char *argv[]);
};

#endif