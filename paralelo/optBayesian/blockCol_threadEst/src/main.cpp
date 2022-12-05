
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
#include "bayesopt.hpp"    
#include "parameters.hpp"            // For the C++ API
#include <boost/numeric/ublas/assignment.hpp> // <<= op assigment
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
    n_thread = 1024; // Numero de escuelas simultaneos


class paralelOptimization: public bayesopt::ContinuousModel
{
 public:
  paralelOptimization(size_t dim,bayesopt::Parameters param):
    ContinuousModel(dim,param) { }

  double evaluateSample( const vectord &Xi ) 
  {
    double bestSolution;
    double result;
    int count;
    float a1 = 0.6, a2 = 0.4;
    // Transforma Valores
    float var1 = round(Xi(0)*84)+1; // Len1    [1,85]
    float var2  = round(Xi(1)*84)+1; // Len2   [1,85]
    double var3 = min(0.099*Xi(2)+0.900,1.0); // CoolingRate [0.9,0.999]
    double var4 = Xi(3); // k_reheating (0,1] 
    int var5 = round(Xi(4)*999)+1; // n_reheating [1, 1000] 
    int var6 = round(Xi(5)*31)*32+32; // n_thread [] multiplos de 32
    int  var7 = round(Xi(6)*31)*32+32; // n_block [] multiplos de 32
    int var8 = round(Xi(7)*100); // max numero de recalentamiento
    double var9 = round(Xi(7)*100000000)+0.0;
    double var10 = round(Xi(7)*0.99999990)+0.00000009;
    cout << Xi(0) << " " << Xi(1) << " " << Xi(2) << " " << Xi(3) << " " << Xi(4) << endl;
    cout << var1 << " " << var2 << " " << var3 << " " << var4 << " " << var5 << endl;
    cout << Xi(5) << " " << Xi(6) << " " << Xi(7) << " " << Xi(8) << " " << Xi(9) << endl;
    cout << var6 << " " << var7 << " " << var8 << " " << var9 << " " << var10 << endl;
    // Ejecuta
    tie(bestSolution, count)  = sasFunc(var1,var2,var3,alpha1,alpha2,alpha3,var4,var5,var8,var9,var10,var7,var6,-1,1);
    result = (bestSolution*a1)+(a2*((double)count/30000000));
    cout << "Solución Actual: " << bestSolution << " | Iteraciones: "
         << count << " | resultado: " << result <<  " | a1: " << a1 
         <<   " | a2: " << a2 << " (count/30000000): " << ((double)count/30000000) 
         << " (bestSolution*a1): " <<  (bestSolution*a1) << " (a2*(count/30000000)): " <<  (a2*((double)count/30000000)) << endl;
    cout << "---------------------------------------------------------------------------------------------------------" << endl;
    return result;
  };
  bool checkReachability( const boost::numeric::ublas::vector<double> &query )
  { 
    { return true; };
  };
};


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
    int dim = 10;
    bayesopt::Parameters params = initialize_parameters_to_default();

    params.sc_type = SC_MAP;

    params.l_type = L_MCMC; // L_MCMC mayor precición pero mayor tiempo, L_EMPIRICAL mas rapido pero con menor precición
    params.noise = 0.001; 
    params.n_iterations = 60;    // Number of iterations
    params.random_seed = 0; // Si el valor es positivo se usa como semilla para el generador de numeros aleatorios, si es negativo se usa como semilla el tiempo.
    params.n_init_samples = 5; //
    params.n_iter_relearn = 1; 
    params.verbose_level = 4;
    params.log_filename = "../logs/bayesopt"+prefijo_save;// Falta agregar el current time
    
    params.load_save_flag = 2;
    params.save_filename = "../logs/bayesopt"+prefijo_save+".dat";
    paralelOptimization optimizer(dim,params);
    boost::numeric::ublas::vector<double> bestPoint(dim);
    optimizer.optimize(bestPoint);
    std::cout << "Final result: " << bestPoint << std::endl;
    double bestSolution;
    int count;

    return (EXIT_SUCCESS);

}
