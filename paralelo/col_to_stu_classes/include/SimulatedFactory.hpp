#ifndef SIMULATED_FACTORY_HPP
#define SIMULATED_FACTORY_HPP


#include <SimulatedAnnealing.cuh>
#include <map>
#include <functional>

class SimulatedFactory {
    private:
        static std::map<std::string, std::function<AcceptanceCriterion*()>> AcceptanceMap;
        static std::map<std::string, std::function<CoolingScheme*()>> CoolingMap;
        static std::map<std::string, std::function<LengthTemperature*()>> LenghtMap;
        static std::map<std::string, std::function<ReheatingMethod*()>> ReheatingMap;
    public:
        static SimulatedAnnealing* createSimulatedAnnealing(
            string acceptancecriterion,
            string coolingscheme,
            string lengthtemeperature,
            string reheatingmethod);
};

#endif