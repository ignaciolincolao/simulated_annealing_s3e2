#ifndef MAIN_S3E2_H
#define MAIN_S3E2_H


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


using std::stoi;
using std::stod;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::string;
using std::cout;
using std::random_device;
using std::stringstream;
using std::getline;
using std::setprecision;
using std::fixed;


///////////////////////////////////////////////////
/// Estructura de datos de los colegios.
///////////////////////////////////////////////////
struct Info_colegio {
    double latitude = 0.0;
    double longitude = 0.0;
    int num_alu = 0;
    int rbd = 0;
    int prioritario = 0;
};
///////////////////////////////////////////////////
/// Estructura de alumnos
///////////////////////////////////////////////////
struct Info_alu {
    int rbd = 0;
    int sep = 0;
    double latitude = 0.0;
    double longitude = 0.0;
};

extern int n_students, n_colegios;

extern double alpha1; // Alpha de distancia
extern double alpha2; // Alpha de segregación
extern double alpha3; // Alpha de costocupo

extern float len1, len2;
extern double len3, len4;

extern double e_const, Th;
extern double coolingRate; // Tasa de enfriamiento

extern double temp; // Temperatura Inicial
extern double min_temp;// 0.00000009; // Minima temperatura que puede llegar
extern double max_temp;
extern double k_reheating;
extern int n_reheating;
extern int seed;

extern string ruta_save; // Ruta para guardar los archivos
extern string name_exp;
extern double alpha[3]; // Valores del alpha con orden Distancia, Segregación, Costo Cupo
extern std::mt19937 mt;
extern std::uniform_int_distribution<int> dist;
extern std::uniform_int_distribution<int> dist2;
extern double max_dist;
extern double min_dist;
extern double init_dist;

extern char timestr[20];
extern int n_block;
extern int n_thread;
extern int selectThread;
extern int selectBlock;
extern string prefijo_save;
extern string name_exp;

///////////////////////////////////////////////////
/// Funciones generales
///////////////////////////////////////////////////
double sasFunc();
double calCosto(int currentSolution[], double **distMat, const double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]);
double meanDist(const int currentSolution[], double  **distMat);
double sumDist(const int currentSolution[], double  **distMat);
double S(const int currentSolution[],const int alumnosSep[], int totalVuln);
double sumS(const int currentSolution[],const int alumnosSep[], int totalVuln);
double costCupo(int currentSolution[],int cupoArray[]);
double sumCostCupo(int currentSolution[],int cupoArray[]);
void newSolution(int currentSolution[],const int previousSolution[]);
double newSolution_v2(int n_students,int n_colegios,int totalVuln,int aluxcol[],int aluVulxCol[],int cupoArray[],double **distMat, int currentSolution[],const double ptr_alpha[]);
void assignSchoolToArray(int previousSolution[], int bestSolution[], int currentSolution[], Info_colegio *ptr_colegios, Info_alu *ptr_students, int cupoArray[]);
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, double **distMat);
void shuffle(int[],int,std::uniform_int_distribution<int>);
void getDataSchool(std::vector<Info_colegio> &colegios);
void getDataStudents(std::vector<Info_alu> &students, int &totalVuln);
double getMaxDistance(double **distMat);
void normalizedAlpha(double alpha[3]);
void initializeArray(int *aluxcol, int *previousAluxCol, int *bestAluxCol, int *aluVulxCol, int *previousAluVulxCol, int *bestAluVulxCol, int *alumnosSep, std::vector<Info_alu> &students,std::vector<Info_colegio> &colegios);

#endif


