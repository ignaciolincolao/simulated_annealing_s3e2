#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <Dataset.hpp>
#include <SimulatedFactory.hpp>


int seed = time(NULL);
string name_exp = "base";
string ruta_save = "../save/"; // Ruta para guardar los archivos
                               // Valores del alpha con orden Distancia, Segregación, Costo Cupo
random_device rd;
mt19937 mt(rd());
uniform_int_distribution<int> dist(0, 0);
uniform_int_distribution<int> dist2(0, 0);
uniform_real_distribution<double> dist_accepta(0.0, 1.0);
char timestr[20];

int main(int argc, char *argv[])
{

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

    if (argc > 1)
    {
        // Config init
        saParams.temp = stod(argv[1]);     // Temperatura inicial
        saParams.min_temp = stod(argv[2]); // Minima temperatura que puede llegar
        saParams.seed = stoi(argv[3]);     // Semilla inicial
        saParams.alpha1 = stod(argv[4]);   //  de distancia
        saParams.alpha2 = stod(argv[5]);   // Alpha de segregación
        saParams.alpha3 = stod(argv[6]);   // Alpha de costocupo
        // Cooling
        csParams.coolingRate = stod(argv[7]); // Tasa de enfriamiento
        // Reheating
        rtParams.k_reheating = stod(argv[8]);
        rtParams.e_const = stod(argv[9]);
        rtParams.n_reheating = stoi(argv[10]);
        // Temperature Length
        ltParams.len1 = stof(argv[11]);
        ltParams.len2 = stof(argv[12]);
        ltParams.len3 = stod(argv[13]);
        ltParams.len4 = stod(argv[14]);
        // Acceptance Criterion
        acParams.Th = stod(argv[15]);
        // Exploration criterion
        cuParams.n_block = stoi(argv[16]);  // Numero de blockes = numeros de alumnos aleatorios
        cuParams.n_thread = stoi(argv[17]); // Numero de threads por bloque = numeros de
                                            // escuelas aleatorios
        // Ubicacion de archivos
        ruta_save = argv[18];
        prefijo_save = argv[19];
        name_exp = argv[20];
    }
    rtParams.max_temp = std::numeric_limits<double>::max();
    double alpha[3] = {saParams.alpha1, saParams.alpha2, saParams.alpha3};
    alpha[0] = saParams.alpha1;
    alpha[1] = saParams.alpha2;
    alpha[2] = saParams.alpha3;
    mt.seed(time(NULL));

    AcceptanceCriterion *aC = new AC1(saParams, acParams);
    CoolingScheme *cS = new CS2(saParams, csParams);
    LengthTemperature *lT = new TL7(saParams, ltParams);
    ReheatingMethod *rM = new TR0(saParams, rtParams);
    RecordManager *rMgr = new RecordManager(saParams, rMgrParams);
    Dataset *dS = new Dataset("colegios_utm.txt", "alumnos_utm.txt");
    SimulatedAnnealing *simulatedAnneling = new SimulatedAnnealing(aC, cS, lT, rM, dS, rMgr, saParams, cuParams);

    simulatedAnneling->runGPU();
    return (EXIT_SUCCESS);
}
