#ifndef MAIN_S3E2_H
#define MAIN_S3E2_H



#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdint>
#include <stdio.h>
#include <sstream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <string>

using namespace std;

///////////////////////////////////////////////////
/// Estructura de datos de los colegios.
///////////////////////////////////////////////////
///////////////////////////////////////////////////
/// Estructura de datos de los colegios.
///////////////////////////////////////////////////
struct Info_colegio{
    long double latitude = 0.0;
    long double longitude = 0.0;
    int num_alu = 0;
    int rbd = 0;
    int prioritario = 0;
};
///////////////////////////////////////////////////
/// Estructura de alumnos
///////////////////////////////////////////////////
struct Info_alu{
    int rbd = 0;
    int sep = 0;
    long double latitude = 0.0;
    long double longitude = 0.0;
};

extern int n_students, n_colegios;
extern float len1, len2;
extern long double len3, len4;
extern long double e_const, Th;
extern long double coolingRate; // Tasa de enfriamiento
extern long double temp; // Temperatura Inicial
extern long double min_temp;// 0.00000009; // Minima temperatura que puede llegar
extern long double alpha1; // Alpha de distancia
extern long double alpha2; // Alpha de segregación
extern long double alpha3; // Alpha de costocupo
extern long double max_temp;
extern long double k_reheating;
extern int n_reheating;
extern int seed;


extern string ruta_save; // Ruta para guardar los archivos
extern long double alpha[3]; // Valores del alpha con orden Distancia, Segregación, Costo Cupo
extern std::mt19937 mt;
extern std::uniform_int_distribution<int> dist;
extern std::uniform_int_distribution<int> dist2;
extern long double max_dist;
extern long double min_dist;
extern long double init_dist;

extern char timestr[20];
extern int n_block;
extern int n_thread;
extern string prefijo_save;
extern string name_exp;

///////////////////////////////////////////////////
/// Funciones generales
///////////////////////////////////////////////////
long double sasFunc();
long double calCosto(int currentSolution[], long double **distMat, const long double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]);
long double meanDist(const int currentSolution[], long double  **distMat);
long double sumDist(const int currentSolution[], long double  **distMat);
long double S(const int currentSolution[],const int alumnosSep[], int totalVuln);
long double sumS(const int currentSolution[],const int alumnosSep[], int totalVuln);
long double costCupo(int currentSolution[],int cupoArray[]);
long double sumCostCupo(int currentSolution[],int cupoArray[]);
void newSolution(int currentSolution[],const int previousSolution[]);
long double newSolution_v2(int n_students,int n_colegios,int totalVuln,int aluxcol[],int aluVulxCol[],int cupoArray[],long double **distMat, int currentSolution[],const long double ptr_alpha[]);
void assignSchoolToArray(int previousSolution[], int bestSolution[], int currentSolution[], Info_colegio *ptr_colegios, Info_alu *ptr_students, int cupoArray[]);
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, long double **distMat);
void shuffle(int[],int,std::uniform_int_distribution<int>);
void getDataSchool(std::vector<Info_colegio> &colegios);
void getDataStudents(std::vector<Info_alu> &students, int &totalVuln);
long double getMaxDistance(long double **distMat);
void normalizedAlpha(long double alpha[3]);
void initializeArray(int *aluxcol, int *previousAluxCol, int *bestAluxCol, int *aluVulxCol, int *previousAluVulxCol, int *bestAluVulxCol, int *alumnosSep, std::vector<Info_alu> &students,std::vector<Info_colegio> &colegios);
long double round_n(long double x, int n);

#endif


