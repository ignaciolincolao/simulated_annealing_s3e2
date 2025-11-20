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
int count = 0;

void algorithm_sample(const double config[4], int seed,  string timestr, string pathSave,int argc, char *argv[]){
    random_device rd;
    mt19937 mt(rd());

    RecordParams* rMgrParams = new RecordParams{
                .prefijo_save = string(timestr),
                .ruta_save = "../save/",
                //.ruta_save = "../../save/",
                .name_exp = "base",
                .activated_files = {false,true,false,false,false}};

    //seed = mt();
    SimulatedParams* saParams = new SimulatedParams{
        .seed = seed,
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
        .alpha1 = config[0],
        .alpha2 = config[1],
        .alpha3 = config[2],
        .alpha4 = config[3],
        .max_dist = 0.0,
        .min_dist = 0.0,
        .init_dist = 0.0,
        .costPrevious = 0.0,
        .costCurrent = 0.0,
        .alpha = {config[0], config[1], config[2], config[3]},
        .max_choices = 12,
        .penalty_curve = 0.5f,      //parametro que indica que tan curva es la maquina
        .penalty_max_pref = 0.5f   //maxima penalidad por la ultima preferencia
    };

    AcceptanceParams* acParams = new AcceptanceParams{
        .Th = 1.1};
    CoolingParams* csParams = new CoolingParams{
        .coolingRate = 0.98};
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
        .n_block = 31,
        .n_thread = 32,
        .selectThread = 0,
        .selectBlock = 0};


    SimulatedStruct* simStruct = new SimulatedStruct{
        .acceptancecriterion = "AC1",
        .coolingscheme = "CS2",
        .lengthtemperature ="TL7",
        .reheatingmethod = "TR0"
    };

    if(argc < 2){
        temp= 1602.26;
        coolingRate= 0.98;
        len1= 2.78021;
        len2= 9.89461;
    }
    else{
        temp= std::stod(argv[1]);
        coolingRate= std::stod(argv[2]);
        len1= std::stod(argv[3]);
        len2= std::stod(argv[4]);

    }
    saParams->temp = temp;
    csParams->coolingRate = coolingRate;
    ltParams->len1 = len1;
    ltParams->len2 = len2;
    cout<< " | n_iters= " << count
        << " | a_dist= " << config[0]
        << " | a_seg= "  << config[1]
        << " | a_costcup= " << config[2]
        << " | a_penalty= " << config[3]
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
    count++;


    std::ifstream inFile(pathSave);
    bool isEmpty = inFile.peek() == std::ifstream::traits_type::eof();
    inFile.close();
    std::ofstream fileData(pathSave, std::ios::app);
    if (!fileData.is_open()) {
        std::cerr << "Error al abrir el archivo para escritura." << std::endl;
        exit(1);
    }
    if (fileData.is_open()) {
        if (isEmpty) {
            fileData << "seed,"
                <<  "z,"
                <<  "it,"
                <<  "n_block,"
                <<  "n_thread,"
                <<  "temp,"
                <<  "coolingRate,"
                <<  "len1,"
                <<  "len2";
                for (int i=0; i < simulatedAnneling->recordManager->vector_percentage.size(); i++){
                    fileData << "," << simulatedAnneling->recordManager->vector_percentage.at(i);
                }
                fileData << endl;
        }
        fileData << seed << ","
                    << val << ","
                    << it << ","
                    << n_block << ","
                    << n_thread << ","
                    << temp << ","
                    << coolingRate << ","
                    << len1 << ","
                    << len2;
        for (int i=0; i < simulatedAnneling->recordManager->vector_percentage.size(); i++){
            fileData << "," << simulatedAnneling->recordManager->vector_it_percentage.at(i);
        }
        fileData << endl;
    }
    else{
        std::cerr << "No se pudo abrir el archivo " << pathSave << std::endl;
        }





    
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
    //const std::string file_name = "../../save/"+string(timestr)+"seed_iteration.csv";
    int init = 0; 




    

    double delta_in = 0.10;  // paso deseado (puedes cambiarlo)

    // Elegimos L redondeando 0.6/delta_in y definimos delta_eff := 0.6/L
    // para garantizar suma exacta 1: 0.4 + L*delta_eff = 1.0
    double L_real = 0.6 / delta_in;
    long L = lround(L_real);               // redondeo al entero más cercano
    double delta_eff = 0.6 / (double)L;    // paso ajustado exacto
    double eps = 1e-12;
    // Rejilla en el simplex truncado: p_i = 0.1 + k_i*delta_eff, sum k_i = L
    // Iteramos con 3 bucles y cerramos con k4 = L - k1 - k2 - k3.

    int seed = 1000000;
    for (long k1 = 0; k1 <= L; ++k1) {
        for (long k2 = 0; k2 <= L - k1; ++k2) {
            for (long k3 = 0; k3 <= L - k1 - k2; ++k3) {
                long k4 = L - k1 - k2 - k3;

                double p1 = 0.1 + k1 * delta_eff;
                double p2 = 0.1 + k2 * delta_eff;
                double p3 = 0.1 + k3 * delta_eff;
                double p4 = 0.1 + k4 * delta_eff;

                // Verificación numérica de suma 1 (tolerancia eps)
                for (int i=0; i<50; i++){
                    if (count < init){
                        continue;
                    }else{
                        double config[4] = {p1, p2, p3, p4};
                        algorithm_sample(config, seed, timestr, file_name,argc,argv);
                        seed ++;
                    }
                }
            }
        }
    }



    
    return (EXIT_SUCCESS);
}
