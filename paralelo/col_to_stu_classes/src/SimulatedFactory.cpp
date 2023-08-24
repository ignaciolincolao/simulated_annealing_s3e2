#include<SimulatedFactory.hpp>

std::map<std::string, std::function<AcceptanceCriterion*()>> SimulatedFactory::AcceptanceMap{
    {"AC1", []() { return new AC1(saParams, acParams) }},
    {"AC3", []() { return new AC3(saParams, acParams) }},
    {"AC6", []() { return new AC6(saParams, acParams) }}
};

std::map<std::string, std::function<CoolingScheme*()>> SimulatedFactory::CoolingMap{
    {"CS2", []() { return new CS2(saParams, csParams); }},
};

std::map<std::string, std::function<LengthTemperature*()>> SimulatedFactory::LenghtMap{
    {"TL7", []() { return new TL7(saParams, ltParams); }},
    {"TL8", []() { return new TL8(saParams, ltParams); }},
    {"TL9", []() { return new TL9(saParams, ltParams); }},
    {"TL11", []() { return new TL11(saParams, ltParams); }}
};

std::map<std::string, std::function<ReheatingMethod*()>> SimulatedFactory::ReheatingMap{
    {"TR11", []() { return new TR11(saParams, rtParams); }},
    {"TR12", []() { return new TR12(saParams, rtParams); }},
    {"TR13", []() { return new TR13(saParams, rtParams); }},
    {"TR14", []() { return new TR14(saParams, rtParams); }}
};



SimulatedAnnealing* SimulatedFactory::createSimulatedAnnealing(

            string acceptancecriterion,
            string coolingscheme,
            string lengthtemperature,
            string reheatingmethod){
                switch (acceptancecriterion.at)
                {
                case "AC1":
                    /* code */
                    break;
                
                default:
                    break;
                }

            }