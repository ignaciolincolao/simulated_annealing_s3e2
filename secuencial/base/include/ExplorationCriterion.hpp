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
extern long double newSolution_v2(int n_students,int n_colegios,int totalVuln,int aluxcol[],int aluVulxCol[],int cupoArray[],long double **distMat, int currentSolution[], const long double ptr_alpha[]);
extern void shuffle(int values[], const int max_change, std::uniform_int_distribution<int> distri);
extern long double calCosto(int currentSolution[], long double **distMat, const long double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]);
extern const int DECIMAL;
extern long double round_n(long double x);

long double solutionNE1(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        long double **distMat, 
        int *currentSolution,
        long double costCurrentSolution,
        const long double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep);

long double solutionNE1_v2(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        long double **distMat, 
        int *currentSolution,
        long double costCurrentSolution,
        const long double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep,
        long double *currentVars,
        long double max_dist);

long double solutionNE3(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        long double **distMat, 
        int *currentSolution,
        long double costCurrentSolution,
        const long double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep);

long double solutionNE4(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        long double **distMat, 
        int *currentSolution,
        long double costCurrentSolution,
        const long double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep);



#endif