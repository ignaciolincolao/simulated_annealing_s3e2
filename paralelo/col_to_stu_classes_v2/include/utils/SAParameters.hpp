#ifndef SA_PARAMETERS_HPP
#define SA_PARAMETERS_HPP

class SimulatedAnnealing;

struct SimulatedParams {
    int seed;
    int n_students;
    int n_colegios;
    int count_rechaso;
    int count;
    int count_trials;
    int c_cooling_temperature;
    int c_accepta=0;
    int max_changes_school = 0;
    int max_changes_students = 0;
    double temp; // Temperatura Inicial
    double temp_init;
    double min_temp;// 0.00000009; // Minima temperatura que puede llegar
    double alpha1; // Alpha de distancia
    double alpha2; // Alpha de segregación
    double alpha3; // Alpha de costocupo
    double max_dist;
    double min_dist;
    double init_dist;
    double costPrevious=0.0;
    double costCurrent=0.0;
    int* shuffle_student;
    int* shuffle_colegios;
    double alpha[3] = {alpha1, alpha2, alpha3}; // Valores del alpha con orden Distancia, Segregación, Costo Cupo



};

struct CUDAParams {   
    int n_block;
    int n_thread;
    int selectThread;
    int selectBlock;
};


#endif