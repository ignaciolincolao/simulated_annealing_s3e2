#ifndef SIMULATED_FACTORY_HPP
#define SIMULATED_FACTORY_HPP


#include <SimulatedAnnealing.cuh>
#include <map>
#include <functional>

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