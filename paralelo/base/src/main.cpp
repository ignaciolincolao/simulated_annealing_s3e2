#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <SimulatedFactory.hpp>





int main(int argc, char *argv[])
{

    random_device rd;
    mt19937 mt(rd());

    // Hora Actual
    time_t hora_actual;
    struct tm *time_info;
    time(&hora_actual);
    time_info = localtime(&hora_actual);
    char timestr[20];
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T_%H-%M", time_info);

    // Configuración del algoritmo
    RecordParams* rMgrParams = new RecordParams{
                .prefijo_save = string(timestr),
                .ruta_save = "../../save/",
                //.ruta_save = "./save/",
                .name_exp = "base",
                .activated_files = {true,true,true,true,true}};

    double alp1 = 0.2; //distancia
    double alp2 = 0.1; //segregacion
    double alp3 = 0.3; //costocupo
    double alp4 = 0.4; //penalty parents
    SimulatedParams* saParams = new SimulatedParams{
        //.seed = 1574067955,
        .seed = 1574067956,
        .n_students = 0,
        .n_colegios = 0,
        .count_rechaso = 0,
        .count = 0,
        .c_cooling_temperature = 0,
        .c_accepta = 0,
        .p = 0.00010,
        .k = 0.01,
        .pMax = 0.3,
        .pInit = 0.01,
        .temp = 32768.0,
        .min_temp = 0.00000009,
        .alpha1 = alp1,
        .alpha2 = alp2,
        .alpha3 = alp3,
        .alpha4 = alp4,
        .max_dist = 0.0,
        .min_dist = 0.0,
        .init_dist = 0.0,
        .costPrevious = 0.0,
        .costCurrent = 0.0,
        .alpha = {alp1, alp2, alp3, alp4},
        //.max_choices = 14 //cantidad maxima de elecciones de padres a considerar
        .max_choices = 14,
        .penalty_curve = 0.5f,      //parametro que indica que tan curva es la maquina
        .penalty_max_pref = 0.5f,   //maxima penalidad por la ultima preferencia
        .costCupoFactorAlpha =0.8f //valor minimo que toma el factor alpha del costcupo (calculo original)
    };

    AcceptanceParams* acParams = new AcceptanceParams{
        .Th = 1.1};
    CoolingParams* csParams = new CoolingParams{
        .coolingRate = 0.98};
    LengthParams* ltParams = new LengthParams{
        .len1 = 5,
        .len2 = 5,
        .len3 = 1.0,
        .len4 = 0.999};
    ReheatingParams* rtParams = new ReheatingParams{
        .e_const = 0.01,
        .max_temp = std::numeric_limits<double>::max(),
        .k_reheating = 30,
        .n_reheating = 1,
        .k_reheating_init = 0};

    CUDAParams* cuParams = new CUDAParams{
        .n_block = 94,
        .n_thread = 32,
        .selectThread = 0,
        .selectBlock = 0};


    SimulatedStruct* simStruct = new SimulatedStruct{
        .acceptancecriterion = "AC1",
        .coolingscheme = "CS2",
        .lengthtemperature ="TL7",
        .reheatingmethod = "TR0"
    };

    SimulatedAnnealing *simulatedAnneling = SimulatedFactory::createSimulatedAnnealing(
            simStruct,
            rMgrParams,
            saParams,
            acParams,
            csParams,
            ltParams,
            rtParams,
            cuParams,
            mt);

    simulatedAnneling->runGPU();
    //simulatedAnneling->ValidateGPU();
    delete simulatedAnneling;
    delete simStruct;
    delete rMgrParams;
    delete saParams;
    delete acParams;
    delete csParams;
    delete ltParams;
    delete rtParams;
    delete cuParams;
    

    
    
    return (EXIT_SUCCESS);
}
