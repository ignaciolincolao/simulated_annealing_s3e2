#include<SimulatedFactory.hpp>

std::map<std::string, std::function<AcceptanceCriterion*(SimulatedParams& saParams, AcceptanceParams& acParams)>> SimulatedFactory::AcceptanceMap{
    {"AC1", [](SimulatedParams& saParams, AcceptanceParams& acParams) { return new AC1(saParams, acParams); }},
    {"AC3", [](SimulatedParams& saParams, AcceptanceParams& acParams) { return new AC3(saParams, acParams); }},
    {"AC6", [](SimulatedParams& saParams, AcceptanceParams& acParams) { return new AC6(saParams, acParams); }}
};

std::map<std::string, std::function<CoolingScheme*(SimulatedParams &saParams, CoolingParams &csParams)>> SimulatedFactory::CoolingMap{
    {"CS2", [](SimulatedParams &saParams, CoolingParams &csParams) { return new CS2(saParams, csParams); }},
};

std::map<std::string, std::function<LengthTemperature*(SimulatedParams& saParams,LengthParams& ltParams)>> SimulatedFactory::LenghtMap{
    {"TL7", [](SimulatedParams& saParams,LengthParams& ltParams) { return new TL7(saParams, ltParams); }},
    {"TL8", [](SimulatedParams& saParams,LengthParams& ltParams) { return new TL8(saParams, ltParams); }},
    {"TL9", [](SimulatedParams& saParams,LengthParams& ltParams) { return new TL9(saParams, ltParams); }},
    {"TL11",[](SimulatedParams& saParams,LengthParams& ltParams) { return new TL11(saParams, ltParams); }}
};

std::map<std::string, std::function<ReheatingMethod*(SimulatedParams& saParams,ReheatingParams& rmParams)>> SimulatedFactory::ReheatingMap{
    {"TR11", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR11(saParams, rmParams); }},
    {"TR12", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR12(saParams, rmParams); }},
    {"TR13", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR13(saParams, rmParams); }},
    {"TR14", [](SimulatedParams& saParams,ReheatingParams& rmParams) { return new TR14(saParams, rmParams); }}
};



SimulatedAnnealing* SimulatedFactory::createSimulatedAnnealing(

            string acceptancecriterion,
            string coolingscheme,
            string lengthtemperature,
            string reheatingmethod,
            int argc,
            char *argv[]){
            string name_exp = "base";
            string ruta_save = "../save/"; // Ruta para guardar los archivos
                                        // Valores del alpha con orden Distancia, Segregación, Costo Cupo
            std::string prefijo_save;
            time_t hora_actual;
            struct tm *time_info;
            time(&hora_actual);
            time_info = localtime(&hora_actual);
            char timestr[20];
            strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);
            prefijo_save = string(timestr);
            RecordParams rMgrParams = {
                .prefijo_save = prefijo_save,
                .ruta_save = "../save/",
                .name_exp = name_exp};

            SimulatedParams saParams;

            saParams = {
                .seed = 123456,
                .n_students = 0,
                .n_colegios = 0,
                .count_rechaso = 0,
                .count = 0,
                .c_cooling_temperature = 0,
                .c_accepta = 0,
                .temp = 1.0,
                .min_temp = 0.0000009,
                .alpha1 = 15.0,
                .alpha2 = 30.0,
                .alpha3 = 25.0,
                .max_dist = 0.0,
                .min_dist = 0.0,
                .init_dist = 0.0,
                .costPrevious = 0.0,
                .costCurrent = 0.0};

            AcceptanceParams acParams = {
                .Th = 1.1};
            CoolingParams csParams = {
                .coolingRate = 0.96};
            LengthParams ltParams = {
                .len1 = 1,
                .len2 = 2,
                .len3 = 1.0,
                .len4 = 0.99};
            ReheatingParams rtParams = {
                .e_const = 0.01,
                .max_temp = 0.0,
                .k_reheating = 30,
                .n_reheating = 1,
                .k_reheating_init = 0};

            CUDAParams cuParams = {
                .n_block = 48,
                .n_thread = 32,
                .selectThread = 0,
                .selectBlock = 0};

            AcceptanceCriterion *aC; 
            CoolingScheme *cS;
            LengthTemperature *lT;
            ReheatingMethod *rM;
            if (AcceptanceMap.find(acceptancecriterion) != AcceptanceMap.end()) {
                aC  = AcceptanceMap[acceptancecriterion](saParams, acParams);
                // haz algo con "acceptance"
            }
            
            if (CoolingMap.find(coolingscheme) != CoolingMap.end()) {
                cS = CoolingMap[coolingscheme](saParams, csParams);
                // haz algo con "cooling"
            }
            
            if (LenghtMap.find(lengthtemperature) != LenghtMap.end()) {
                lT = LenghtMap[lengthtemperature](saParams, ltParams);
                // haz algo con "lengthTemp"
            }
            
            if (ReheatingMap.find(reheatingmethod) != ReheatingMap.end()) {
                rM = ReheatingMap[reheatingmethod](saParams, rtParams);
                // haz algo con "reheating"
            }
            RecordManager *rMgr = new RecordManager(saParams, rMgrParams);
            Dataset *dS = new Dataset("colegios_utm.txt", "alumnos_utm.txt");
            SimulatedAnnealing *simulatedAnneling = new SimulatedAnnealing(aC, cS, lT, rM, dS, rMgr, saParams, cuParams);
            return simulatedAnneling;

}