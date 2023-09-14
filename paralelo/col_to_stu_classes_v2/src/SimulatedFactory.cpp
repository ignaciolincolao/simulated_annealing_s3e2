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

            string acceptancecriterion,
            string coolingscheme,
            string lengthtemperature,
            string reheatingmethod,
            int argc,
            char *argv[],
            mt19937 &mt){
            string name_exp = "base";
            string ruta_save = "../save/"; 
            std::string prefijo_save;
            time_t hora_actual;
            struct tm *time_info;
            time(&hora_actual);
            time_info = localtime(&hora_actual);
            char timestr[20];
            strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);
            prefijo_save = string(timestr);

            RecordParams* rMgrParams = new RecordParams{
                .prefijo_save = prefijo_save,
                .ruta_save = "../save/",
                .name_exp = name_exp};

            
            SimulatedParams* saParams = new SimulatedParams{
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

            AcceptanceParams* acParams = new AcceptanceParams{
                .Th = 1.1};
            CoolingParams* csParams = new CoolingParams{
                .coolingRate = 0.96};
            LengthParams* ltParams = new LengthParams{
                .len1 = 1,
                .len2 = 2,
                .len3 = 1.0,
                .len4 = 0.99};
            ReheatingParams* rtParams = new ReheatingParams{
                .e_const = 0.01,
                .max_temp = 0.0,
                .k_reheating = 30,
                .n_reheating = 1,
                .k_reheating_init = 0};

            CUDAParams* cuParams = new CUDAParams{
                .n_block = 48,
                .n_thread = 32,
                .selectThread = 0,
                .selectBlock = 0};

            if (argc > 1)
            {
                // Config init
                saParams->temp = stod(argv[1]);     // Temperatura inicial
                saParams->min_temp = stod(argv[2]); // Minima temperatura que puede llegar
                saParams->seed = stoi(argv[3]);     // Semilla inicial
                saParams->alpha1 = stod(argv[4]);   //  de distancia
                saParams->alpha2 = stod(argv[5]);   // Alpha de segregación
                saParams->alpha3 = stod(argv[6]);   // Alpha de costocupo
                // Cooling
                csParams->coolingRate = stod(argv[7]); // Tasa de enfriamiento
                // Reheating
                rtParams->k_reheating = stod(argv[8]);
                rtParams->e_const = stod(argv[9]);
                rtParams->n_reheating = stoi(argv[10]);
                // Temperature Length
                ltParams->len1 = stof(argv[11]);
                ltParams->len2 = stof(argv[12]);
                ltParams->len3 = stod(argv[13]);
                ltParams->len4 = stod(argv[14]);
                // Acceptance Criterion
                acParams->Th = stod(argv[15]);
                // Exploration criterion
                cuParams->n_block = stoi(argv[16]);  // Numero de blockes = numeros de alumnos aleatorios
                cuParams->n_thread = stoi(argv[17]); // Numero de threads por bloque = numeros de
                                                    // escuelas aleatorios
                // Ubicacion de archivos
                ruta_save = argv[18];
                prefijo_save = argv[19];
                name_exp = argv[20];
                cout << "entro" << endl;
            }
            rtParams->max_temp = std::numeric_limits<double>::max();

            AcceptanceCriterion *aC; 
            CoolingScheme *cS;
            LengthTemperature *lT;
            ReheatingMethod *rM;
            if (AcceptanceMap.find(acceptancecriterion) != AcceptanceMap.end()) {
                aC  = AcceptanceMap[acceptancecriterion](*saParams, *acParams, mt);
                // haz algo con "acceptance"
            }
            
            if (CoolingMap.find(coolingscheme) != CoolingMap.end()) {
                cS = CoolingMap[coolingscheme](*saParams, *csParams);
                // haz algo con "cooling"
            }
            
            if (LenghtMap.find(lengthtemperature) != LenghtMap.end()) {
                lT = LenghtMap[lengthtemperature](*saParams, *ltParams);
                // haz algo con "lengthTemp"
            }
            
            if (ReheatingMap.find(reheatingmethod) != ReheatingMap.end()) {
                rM = ReheatingMap[reheatingmethod](*saParams, *rtParams);
                // haz algo con "reheating"
            }
            RecordManager *rMgr = new RecordManager(*saParams, *rMgrParams);
            Dataset *dS = new Dataset("colegios_utm.txt", "alumnos_utm.txt");
            SimulatedAnnealing *simulatedAnneling = new SimulatedAnnealing(aC, cS, lT, rM, dS, rMgr, saParams, cuParams, mt);
            cout << simulatedAnneling->saParams.seed << endl;
            return simulatedAnneling;

}