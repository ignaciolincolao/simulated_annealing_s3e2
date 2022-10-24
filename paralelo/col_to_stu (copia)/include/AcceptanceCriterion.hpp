#ifndef ACCEPTANCE_CRITERION_H
#define ACCEPTANCE_CRITERION_H

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <random>

using namespace std;

extern std::mt19937 mt;
extern double temp;
int metropolisAC1(double costPrevious, double costCurrent);
int metropolisAC3(double costPrevious, double costCurrent,double Th);
int dCriteriaAC6(double costPrevious, double costCurrent);
double p(double costPrevious,double costCurrent);

#endif