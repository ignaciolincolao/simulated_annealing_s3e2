#ifndef EXPLORATION_CRITERION_H
#define EXPLORATION_CRITERION_H

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <random>
#include <cstdint>
#include <cstring>
#include <iomanip> 
using namespace std;

extern int n_block;
extern int n_thread;
extern mt19937 mt;
extern uniform_int_distribution<int> dist;
extern uniform_int_distribution<int> dist2;
extern double newSolution_v2(int n_students,int n_colegios,int totalVuln,int aluxcol[],int aluVulxCol[],int cupoArray[],double **distMat, int currentSolution[], const double ptr_alpha[]);
extern void shuffle(int values[], const int max_change, std::uniform_int_distribution<int> distri);
extern double calCosto(int currentSolution[], double **distMat, const double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]);
extern const int DECIMAL;
extern double round_n(double x);

double solutionNE1(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep);

double solutionNE1_v2(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep,
        double *currentVars,
        double max_dist);

double solutionNE3(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep);

double solutionNE4(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep);



#endif