#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H


#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdint>
#include <future>
#include <stdio.h>
#include <sstream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <string>
#include <structure/AcceptanceCriterion/AcceptanceCriterion.hpp>
#include <structure/CoolingScheme/CoolingScheme.hpp>
#include <structure/ExplorationCriterion/ExplorationCriterion.hpp>
#include <structure/LengthTemperature/LengthTemperature.hpp>
#include <structure/ReheatingMethod/ReheatingMethod.hpp>
#include <RecordManager.hpp>
#include <Dataset.hpp>

using std::string;
using std::stof;
using std::stoi;
using std::stod;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::random_device;
using std::stringstream;
using std::getline;
using std::setprecision;
using std::fixed;
using std::mt19937;
using std::random_device;
using std::uniform_int_distribution;
using std::uniform_real_distribution;




class SimulatedAnnealing {

private:
    AcceptanceCriterion* acceptanceCriterion;
    CoolingScheme* coolingScheme;
    LengthTemperature* lengthTemperature;
    ReheatingMethod* reheatingMethod;
    Dataset* dataSet;
    RecordManager* recordManager;

    int* previousSolution;
    int* bestSolution;
    int* currentSolution;
    int* cupoArray;
    int* alumnosSep;
    int* aluxcol;
    int* aluVulxCol;
    int* previousAluxCol;
    int* previousAluVulxCol;
    int* bestAluxCol;
    int* bestAluVulxCol;
    double* alpha;
    double* currentVars;
    double* previousVars;
    double* bestVars;
    double* ptr_alpha;
    double *matrestest;
    double **distMat;
    double *probSelection;
    mt19937& mt;
    

public:
    AcceptanceParams& acParams;
    CoolingParams& csParams;
    LengthParams& ltParams;
    ReheatingParams& rmParams;
    RecordParams& rmgrParams;
    SimulatedParams& saParams;
    CUDAParams& cuParams;
    int totalVuln;
    double costCurrentSolution;
    double costBestSolution;
    double costPreviousSolution;
    
    uniform_int_distribution<int> dist;
    uniform_int_distribution<int> dist2;
    uniform_real_distribution<double> dist_accepta;
    std::vector<std::future<void>> futures;

    SimulatedAnnealing(AcceptanceCriterion* AC,
        CoolingScheme* CS,
        LengthTemperature* LT,
        ReheatingMethod* RM,
        Dataset* DS,
        RecordManager* RMgr,
        SimulatedParams* saParams_,
        CUDAParams* cuParams_,
        mt19937& mt
    );
    ~SimulatedAnnealing();
    SimulatedParams& getSaParams() { return saParams; };
    double runGPU();
    template <typename T>
    void inicializationValues(T* wrapper);
    double calCosto(int *currentSolution, double **distMat, const double *ptr_alpha, int *alumnosSep, int totalVuln, int *cupoArray);
    double meanDist(const int *currentSolution, double  **distMat);
    double sumDist(const int *currentSolution, double  **distMat);
    double S(const int *currentSolution,const int *alumnosSep, int totalVuln);
    double sumS(const int *currentSolution,const int *alumnosSep, int totalVuln);
    double costCupo(int *currentSolution,int *cupoArray);
    double sumCostCupo(int *currentSolution,int *cupoArray);
    void newSolution(int *currentSolution,const int *previousSolution);
    void assignSchoolToArray(int *previousSolution, int *bestSolution, int *currentSolution, Info_colegio *ptr_colegios, Info_alu *ptr_students, int *cupoArray);
    void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, double **distMat);
    void shuffle(int *values, const int max_change, uniform_int_distribution<int> distri);
    double getMaxDistance(double **distMat);
    void normalizedAlpha(double *alpha);
    void initializeArray(int *aluxcol, int *previousAluxCol, int *bestAluxCol, int *aluVulxCol, int *previousAluVulxCol, int *bestAluVulxCol, int *alumnosSep, vector<Info_alu> &students,vector<Info_colegio> &colegios);
    double round_n(double x);
    int acceptanceCriterionApply();
    int selecSolution();
    void UpdateProb(int it);
};

#endif