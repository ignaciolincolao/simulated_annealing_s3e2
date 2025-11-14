#ifndef CONFIG_H
#define CONFIG_H


#include <structure/AcceptanceCriterion/AcceptanceCriterion.hpp>
#include <structure/CoolingScheme/CoolingScheme.hpp>
#include <structure/LengthTemperature/LengthTemperature.hpp>
#include <structure/ReheatingMethod/ReheatingMethod.hpp>

using std::mt19937;
using std::random_device;
using std::stof;
using std::string;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

///////////////////////////////////////////////////
/// Parametros de configuración Default
///////////////////////////////////////////////////



// double temp = 10000000000000.0; // Temperatura inicial
// double min_temp = 0.00000009;   // Minima temperatura que puede llegar
double min_temp = 0.0000009; // Minima temperatura que puede llegar
double max_temp = 0;
double k_reheating = 30;
int n_reheating = 1; // Variable ligada a cuanto debe esperar para iniciar recalentamiento
int seed = 12315;

string ruta_save = "../save/"; // Ruta para guardar los archivos
string name_exp = "base";
double alpha[3] = {alpha1, alpha2, alpha3}; // Valores del alpha con orden Distancia, Segregación, Costo Cupo
double max_dist = 0.0;
double min_dist = 0.0;
double init_dist = 0.0;
char timestr[20];
string prefijo_save;

uniform_real_distribution<double> dist_accepta(0.0, 1.0);

#endif