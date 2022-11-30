
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




///////////////////////////////////////////////////
/// Funciones generales
///////////////////////////////////////////////////



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
    float var1 = round(Xi(0)*9)+1; // Len1    [1,10]
    float var2  = round(Xi(1)*84)+1; // Len2   [1,85]
    double var3 = 0.099*Xi(2)+0.900; // CoolingRate [0.9,0.999]
    double var4 = Xi(3)+0.001; // k_reheating (0,1] 
    int var5 = round(Xi(4)*99)+1; // n_reheating [1, 100] 
    int var6 = round(Xi(5)*16)*32; // n_thread [] multiplos de 32
    int  var7 = round(Xi(6)*16)*32; // n_block [] multiplos de 32
    cout << Xi(0) << " " << Xi(1) << " " << Xi(2) << " " << Xi(3) << " " << Xi(4) << " " << Xi(5) << " " << Xi(6) << endl;
    cout << var1 << " " << var2 << " " << var3 << " " << var4 << " " << var5 << " " << var6 << " " << var7 << endl;
    // Ejecuta
    tie(bestSolution, count) = sasFunc(var1, var2, var3, var4, var5, var6, var7);
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
    char timestr[20];
    string prefijo_save;
    time_t hora_actual;
    struct tm * time_info;
    time(&hora_actual);
    time_info = localtime(&hora_actual);
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);
    prefijo_save = string(timestr);
    /*
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
    alpha[0]=alpha1;
    alpha[1]=alpha2;
    alpha[2]=alpha3;
    mt.seed(seed);
    */


    int dim = 7;
    bayesopt::Parameters params = initialize_parameters_to_default();
    // Configuración de parametros
    //params.kernel.name = "kSum(kSEISO,kConst)";
    //params.kernel.hp_mean <<= 1.0, 1.0;
    //params.kernel.hp_std <<= 1.0, 1.0;

    //params.mean.name = "mConst";
    //params.mean.coef_mean <<= 1.0;
    //params.mean.coef_std <<= 1.0;
    

    //params.surr_name = "sStudentTProcessJef";
    

    params.sc_type = SC_MAP;

    params.l_type = L_MCMC; // L_MCMC mayor precición pero mayor tiempo, L_EMPIRICAL mas rapido pero con menor precición
    params.noise = 0.001; 
    params.n_iterations = 100;    // Number of iterations
    params.random_seed = 0; // Si el valor es positivo se usa como semilla para el generador de numeros aleatorios, si es negativo se usa como semilla el tiempo.
    params.n_init_samples = 10; //
    params.n_iter_relearn = 1; 
    params.verbose_level = 4;
    params.log_filename = "../logs/bayesopt"+prefijo_save;// Falta agregar el current time




    paralelOptimization optimizer(dim,params);
    //Define bounds and prepare result.
    boost::numeric::ublas::vector<double> bestPoint(dim);
    boost::numeric::ublas::vector<double> lowerBound(dim);
    boost::numeric::ublas::vector<double> upperBound(dim);
    //Set the bounds. This is optional. Default is [0,1]
    //Only required because we are doing continuous optimization
    //optimizer.setBoundingBox(lowerBounds,upperBounds);
    //Collect the result in bestPoint
    optimizer.optimize(bestPoint);
    std::cout << "Final result: " << bestPoint << std::endl;


    

    return (EXIT_SUCCESS);

}
