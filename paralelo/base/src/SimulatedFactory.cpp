#include<SimulatedFactory.hpp>



std::map<std::string, std::function<AcceptanceCriterion*(SimulatedParams& saParams, AcceptanceParams& acParams, std::mt19937& mt)>> SimulatedFactory::AcceptanceMap{
    {"AC1", [](SimulatedParams& saParams, AcceptanceParams& acParams, std::mt19937& mt) { return new AC1(saParams, acParams, mt); }},
    {"AC3", [](SimulatedParams& saParams, AcceptanceParams& acParams, std::mt19937& mt) { return new AC3(saParams, acParams, mt); }},
    {"AC6", [](SimulatedParams& saParams, AcceptanceParams& acParams, std::mt19937& mt) { return new AC6(saParams, acParams, mt); }}
};

std::map<std::string, std::function<CoolingScheme*(SimulatedParams &saParams, CoolingParams &csParams)>> SimulatedFactory::CoolingMap{
    {"CS2", [](SimulatedParams& saParams, CoolingParams& csParams) { return new CS2(saParams, csParams); }},
};

std::map<std::string, std::function<LengthTemperature*(SimulatedParams& saParams,LengthParams& ltParams)>> SimulatedFactory::LenghtMap{
    {"TL7", [](SimulatedParams& saParams,LengthParams& ltParams) { return new TL7(saParams, ltParams); }},
    {"TL8", [](SimulatedParams& saParams,LengthParams& ltParams) { return new TL8(saParams, ltParams); }},
    {"TL9", [](SimulatedParams& saParams,LengthParams& ltParams) { return new TL9(saParams, ltParams); }},
    {"TL11",[](SimulatedParams& saParams,LengthParams& ltParams) { return new TL11(saParams, ltParams); }}
};

std::map<std::string, std::function<ReheatingMethod*(SimulatedParams& saParams,ReheatingParams& rmParams)>> SimulatedFactory::ReheatingMap{
    {"TR0", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR0(saParams, rmParams); }},
    {"TR11", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR11(saParams, rmParams); }},
    {"TR12", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR12(saParams, rmParams); }},
    {"TR13", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR13(saParams, rmParams); }},
    {"TR14", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR14(saParams, rmParams); }},
};


SimulatedAnnealing* SimulatedFactory::createSimulatedAnnealing(
            SimulatedStruct* simStruct,
            RecordParams* rMgrParams,
            SimulatedParams* saParams,
            AcceptanceParams* acParams,
            CoolingParams* csParams,
            LengthParams* ltParams,
            ReheatingParams* rtParams,
            CUDAParams* cuParams,
            mt19937 &mt){

            AcceptanceCriterion *aC; 
            CoolingScheme *cS;
            LengthTemperature *lT;
            ReheatingMethod *rM;
            if (AcceptanceMap.find(simStruct->acceptancecriterion) != AcceptanceMap.end()) {
                aC  = AcceptanceMap[simStruct->acceptancecriterion](*saParams, *acParams, mt);
            }
            
            if (CoolingMap.find(simStruct->coolingscheme) != CoolingMap.end()) {
                cS = CoolingMap[simStruct->coolingscheme](*saParams, *csParams);
            }
            
            if (LenghtMap.find(simStruct->lengthtemperature) != LenghtMap.end()) {
                lT = LenghtMap[simStruct->lengthtemperature](*saParams, *ltParams);
            }
            
            if (ReheatingMap.find(simStruct->reheatingmethod) != ReheatingMap.end()) {
                rM = ReheatingMap[simStruct->reheatingmethod](*saParams, *rtParams);
            }
            RecordManager *rMgr = new RecordManager(*saParams, *rMgrParams);
            Dataset *dS = new Dataset("colegios_test.txt", "alumnos_test.txt");
            SimulatedAnnealing *simulatedAnneling = new SimulatedAnnealing(aC, cS, lT, rM, dS, rMgr, saParams, cuParams, mt);
            return simulatedAnneling;

}