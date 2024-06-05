#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <SimulatedFactory.hpp>



int seed = 0;
int n_block_min = 1;
int n_block_max = 32;
int n_block_factor = 32;
int n_block_idX = 3;//int n_block_idX = 4;
int n_thread_min = 1;
int n_thread_max = 32;
int n_thread_factor = 32;
int n_thread_idX;
double temp_min = 1000;//100;
double temp_max = 100000;//10000;
int temp_idX = 4;
double coolingRate_min = 0.9;
double coolingRate_max = 0.999;
int coolingRate_idX = 0;//int coolingRate_idX = 1;
float len1_min = 1.f;
float len1_max = 10.f;
int len1_idX = 1;//int len1_idX = 2;
float len2_min = 1.f;
float len2_max = 10.f;
int len2_idX = 2;//int len2_idX = 3;
int n_block;
int n_thread;
int it;
double temp;
double coolingRate;
float len1;
float len2;




void algorithm_sample(int indice, string timestr, string pathSave){
    random_device rd;
    mt19937 mt(rd());
    std::ofstream fileData(pathSave, std::ios::app);
    if (!fileData.is_open()) {
        std::cerr << "Error al abrir el archivo para escritura." << std::endl;
        exit(1);
    }

    RecordParams* rMgrParams = new RecordParams{
                .prefijo_save = string(timestr),
                .ruta_save = "../save/",
                .name_exp = "base"};

    double alp1 =  15.0;
    double alp2 = 30.0;
    double alp3 = 25.0;
    seed = mt();
    SimulatedParams* saParams = new SimulatedParams{
        .seed = seed,
        .n_students = 0,
        .n_colegios = 0,
        .count_rechaso = 0,
        .count = 0,
        .c_cooling_temperature = 0,
        .c_accepta = 0,
        .temp = 1.0,
        .min_temp = 0.0000009,
        .alpha1 = alp1,
        .alpha2 = alp2,
        .alpha3 = alp3,
        .max_dist = 0.0,
        .min_dist = 0.0,
        .init_dist = 0.0,
        .costPrevious = 0.0,
        .costCurrent = 0.0,
        .alpha = {alp1, alp2, alp3}};

    AcceptanceParams* acParams = new AcceptanceParams{
        .Th = 1.1};
    CoolingParams* csParams = new CoolingParams{
        .coolingRate = 0.97};
    LengthParams* ltParams = new LengthParams{
        .len1 = 1,
        .len2 = 2,
        .len3 = 1.0,
        .len4 = 0.999};
    ReheatingParams* rtParams = new ReheatingParams{
        .e_const = 0.01,
        .max_temp = std::numeric_limits<double>::max(),
        .k_reheating = 30,
        .n_reheating = 1,
        .k_reheating_init = 0};

    CUDAParams* cuParams = new CUDAParams{
        .n_block = 32,
        .n_thread = 32,
        .selectThread = 0,
        .selectBlock = 0};


    SimulatedStruct* simStruct = new SimulatedStruct{
        .acceptancecriterion = "AC1",
        .coolingscheme = "CS2",
        .lengthtemperature ="TL7",
        .reheatingmethod = "TR0"
    };


    int block_threads[36][2] = {
                                {1,32}, // 32
                                {2,32}, // 64
                                {1,64}, 
                                {8,32}, // 256
                                {4,64},
                                {2,128},
                                {1,256},
                                {16,32}, // 512
                                {8,64},
                                {4,128},
                                {2,256},
                                {1,512},
                                {32,32}, // 1024
                                {16,64},
                                {8,128},
                                {4,256},
                                {2,512},
                                {1,1024},
                                {64,32},  // 2048
                                {32,64},
                                {16,128},
                                {8,256},
                                {4,512},
                                {2,1024},
                                {128,32}, // 4096
                                {64,64},
                                {32,128},
                                {16,256},
                                {8,512},
                                {4,1024},
                                {256,32}, // 8192
                                {128,64},
                                {64,128},
                                {32,256},
                                {16,512},
                                {8,1024}
                            };
    n_block = block_threads[indice][0];
    n_thread = block_threads[indice][1];
    temp= 1602.26;
    coolingRate= 0.98;
    len1= 2.78021;
    len2= 9.89461;
    cuParams->n_block = n_block;
    cuParams->n_thread = n_thread;
    saParams->temp = temp;
    csParams->coolingRate = coolingRate;
    ltParams->len1 = len1;
    ltParams->len2 = len2;
    cout << "temp= "<< saParams->temp 
        <<" | coolingRate= " << csParams->coolingRate 
        << " | len1= " << ltParams->len1
        << " | len2= " << ltParams->len2
        << " | n_block= " << cuParams->n_block
        << " | n_thread= " << cuParams->n_thread
        << endl;

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

    double val = simulatedAnneling->runGPU();
    it = simulatedAnneling->saParams.count;

    fileData << seed << ","
            << val << ","
            << it << ","
            << n_block << ","
            << n_thread << ","
            << temp << ","
            << coolingRate << ","
            << len1 << ","
            << len2 << endl;
    delete simulatedAnneling;
    delete simStruct;
    delete rMgrParams;
    delete saParams;
    delete acParams;
    delete csParams;
    delete ltParams;
    delete rtParams;
    delete cuParams;
}

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
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);
    const std::string file_name = "../save/"+string(timestr)+"seed_iteration.csv";
    std::ofstream fileData(file_name, std::ios::app);
    fileData << "seed,"
        <<  "z,"
        <<  "it,"
        <<  "n_block,"
        <<  "n_thread,"
        <<  "temp,"
        <<  "coolingRate,"
        <<  "len1,"
        <<  "len2"
        << std::endl;
    fileData.close();

    
    for (int x=0; x < 36; x++){
        for(int y=0;y < 30; y++){
        algorithm_sample(x, timestr, file_name);
        }
    }


    
    
    return (EXIT_SUCCESS);
}
