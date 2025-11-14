#ifndef ACCEPTANCE_CRITERION_H
#define ACCEPTANCE_CRITERION_H

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <random>

using namespace std;

extern std::mt19937 mt;
extern long double temp;
int metropolisAC1(long double costPrevious, long double costCurrent);
int metropolisAC3(long double costPrevious, long double costCurrent,long double Th);
int dCriteriaAC6(long double costPrevious, long double costCurrent);
long double p(long double costPrevious,long double costCurrent);

#endif