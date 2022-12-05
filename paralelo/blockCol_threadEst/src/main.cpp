
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <random>
#include <cstdint>
#include <cstring>
#include <sas.cuh>
#include <tuple>
#include <ctime>
///////////////////////////////////////////////////
///////////////////////////////////////////////////
using namespace std;


int n_students, n_colegios;


///////////////////////////////////////////////////
/// Funciones generales
///////////////////////////////////////////////////

///////////////////////////////////////////////////
/// Parametros de configuración Default
///////////////////////////////////////////////////



double alpha1 = 15; // Alpha de distancia valor 1 < alpha1
double alpha2 = 30; // Alpha de segregación valor 1 < alpha2
double alpha3 = 25; // Alpha de costocupo valor 1 < alpha3
double coolingRate = 0.98; // Tasa de enfriamiento valores entre 0 < coolingRate < 1
double temp = 100000.0; // Temperatura inicial
double min_temp = 0.00000009; // Minima temperatura que puede llegar
double max_temp = 0;
double k_reheating = 1;
int n_reheating = 1; // Variable ligada a cuanto debe esperar para iniciar recalentamiento
int max_reheating = 100;
int seed = 12315;
float len1 =1;// 0.00000009; // Minima temperatura que puede llegar
float len2 =2;
double len3 = 1.0;
double len4 = 0.99;
double e_const=0.01;
double Th = 1.1;
string name_exp= "base";
string ruta_save = "../save/"; // Ruta para guardar los archivos
double alpha[3]={alpha1,alpha2,alpha3}; // Valores del alpha con orden Distancia, Segregación, Costo Cupo
random_device rd;
mt19937 mt(rd());
uniform_int_distribution<int> dist(0,0);
uniform_int_distribution<int> dist2(0,0);
uniform_real_distribution<double> dist_accepta(0.0, 1.0);
double max_dist=0.0;
double min_dist=0.0;
double init_dist=0.0;
char timestr[20];
string prefijo_save;


////////////////////////////////
// VARIABLES GLOBALES PARA CUDA
////////////////////////////////

int selectThread=0,
    selectBlock = 0,
    n_block = 32, // Numero de estudiantes simultaneos
    n_thread = 32; // Numero de escuelas simultaneos



int main(int argc, char *argv[]) {
    time_t hora_actual;
    struct tm * time_info;
    time(&hora_actual);
    time_info = localtime(&hora_actual);
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);
    prefijo_save = string(timestr);
    if (argc>1) {
        // Config init
        temp = stod(argv[1]); // Temperatura inicial
        min_temp = stod(argv[2]); // Minima temperatura que puede llegar
        seed= stoi(argv[3]); // Semilla inicial
        alpha1 = stod(argv[4]); //  de distancia
        alpha2 = stod(argv[5]); // Alpha de segregación
        alpha3 = stod(argv[6]); // Alpha de costocupo
        // Cooling
        coolingRate = stod(argv[7]); // Tasa de enfriamiento
        // Reheating
        k_reheating = stod(argv[8]);
        e_const = stod(argv[9]);
        n_reheating = stoi(argv[10]);
        // Temperature Length
        len1 = stof(argv[11]);
        len2 = stof(argv[12]);
        len3 = stod(argv[13]);
        len4 = stod(argv[14]);
        // Acceptance Criterion
        Th = stod(argv[15]);
        // Exploration criterion
        n_block = stoi(argv[16]); // Numero de blockes = numeros de alumnos aleatorios
        n_thread = stoi(argv[17]); // Numero de threads por bloque = numeros de escuelas aleatorios
        // Ubicacion de archivos
        ruta_save = argv[18];
        prefijo_save = argv[19];
        name_exp = argv[20];
    }
    max_temp= pow(10,300);

    
    double bestSolution;
    int count;
    for(int x=0; x<10; x++){
        prefijo_save = to_string(x)+"_32_32";
        sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,32,32,-1,2);
    }
    for(int x=0; x<10; x++){
        prefijo_save = to_string(x)+"_256_256";
        sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,256,256,-1,2);
    }
    for(int x=0; x<10; x++){
        prefijo_save = to_string(x)+"_512_512";
        sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,512,512,-1,2);
    }
        for(int x=0; x<10; x++){
        prefijo_save = to_string(x)+"_1024_1024";
        sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,1024,1024,-1,2);
    }
    //tie(bestSolution, count)  = sasFunc(len1,len2,coolingRate,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,temp,min_temp,n_block,n_thread,seed,2);
    //tie(bestSolution, count)  = 
    //tie(bestSolution, count)  = sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,64,64,seed,1);
    //tie(bestSolution, count)  = sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,256,256,seed,1);
    //tie(bestSolution, count)  = sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,512,512,seed,1);
    //tie(bestSolution, count)  = sasFunc(1,2,0.98,alpha1,alpha2,alpha3,k_reheating,n_reheating,max_reheating,100000,min_temp,1024,32,seed,1);
    return (EXIT_SUCCESS);

}
