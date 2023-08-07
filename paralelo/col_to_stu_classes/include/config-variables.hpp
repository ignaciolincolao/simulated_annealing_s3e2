#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>

using std::mt19937;
using std::random_device;
using std::stof;
using std::string;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

///////////////////////////////////////////////////
/// Parametros de configuración Default
///////////////////////////////////////////////////

int n_students, n_colegios;

double alpha1 = 15; // Alpha de distancia valor 1 < alpha1
double alpha2 = 30; // Alpha de segregación valor 1 < alpha2
double alpha3 = 25; // Alpha de costocupo valor 1 < alpha3

float len1 = 10; // 0.00000009; // Minima temperatura que puede llegar
float len2 = 85;
double len3 = 1.0;
double len4 = 0.99;

double e_const = 0.01;
double Th = 1.1;
double coolingRate = 0.01; // Tasa de enfriamiento valores entre 0 < coolingRate < 1

// double temp = 10000000000000.0; // Temperatura inicial
double temp = 1000000.0; // Temperatura inicial
// double min_temp = 0.00000009;   // Minima temperatura que puede llegar
double min_temp = 1; // Minima temperatura que puede llegar
double max_temp = 0;
double k_reheating = 30;
int n_reheating = 1; // Variable ligada a cuanto debe esperar para iniciar recalentamiento
int seed = 12315;

string ruta_save = "../save/"; // Ruta para guardar los archivos
string name_exp = "base";
double alpha[3] = {alpha1, alpha2, alpha3}; // Valores del alpha con orden Distancia, Segregación, Costo Cupo
random_device rd;
mt19937 mt(rd());
uniform_int_distribution<int> dist(0, 0);
uniform_int_distribution<int> dist2(0, 0);
double max_dist = 0.0;
double min_dist = 0.0;
double init_dist = 0.0;
char timestr[20];
string prefijo_save;

uniform_real_distribution<double> dist_accepta(0.0, 1.0);

////////////////////////////////
// VARIABLES GLOBALES PARA CUDA
////////////////////////////////

int selectThread = 0,
    selectBlock = 0,
    n_block = 1024,  // Numero de estudiantes simultaneos
    n_thread = 1024; // Numero de escuelas simultaneos

#endif