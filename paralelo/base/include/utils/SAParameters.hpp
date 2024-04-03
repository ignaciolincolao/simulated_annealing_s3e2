#ifndef SA_PARAMETERS_HPP
#define SA_PARAMETERS_HPP

#include <string>

class SimulatedAnnealing;

struct SimulatedStruct {
    std::string acceptancecriterion;
    std::string coolingscheme;
    std::string explorationCriterion;
    std::string lengthtemperature;
    std::string reheatingmethod;
};

struct SimulatedParams {
    int seed;
    int n_students;
    int n_colegios;
    int count_rechaso;
    int count;
    int count_trials;
    int c_cooling_temperature;
    int c_accepta;
    int max_changes_school;
    int max_changes_students;
    double p;
    double k;
    double pMax;
    double pInit;
    double temp; // Temperatura Inicial
    double temp_init;
    double min_temp;// 0.00000009; // Minima temperatura que puede llegar
    double alpha1; // Alpha de distancia
    double alpha2; // Alpha de segregación
    double alpha3; // Alpha de costocupo
    double max_dist;
    double min_dist;
    double init_dist;
    double costPrevious;
    double costCurrent;
    int* shuffle_student;
    int* shuffle_colegios;
    double alpha[3]; // Valores del alpha con orden Distancia, Segregación, Costo Cupo

};

struct CUDAParams {   
    int n_block;
    int n_thread;
    int selectThread;
    int selectBlock;
};


#endif