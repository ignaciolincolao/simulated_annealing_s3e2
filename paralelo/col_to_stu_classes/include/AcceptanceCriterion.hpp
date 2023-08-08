#ifndef ACCEPTANCE_CRITERION_H
#define ACCEPTANCE_CRITERION_H

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <random>

using std::mt19937;
using std::uniform_real_distribution;

extern mt19937 mt;
extern double temp;

/*
int metropolisAC1(double costPrevious, double costCurrent);
int metropolisAC3(double costPrevious, double costCurrent,double Th);
int dCriteriaAC6(double costPrevious, double costCurrent);
double p(double costPrevious,double costCurrent);
*/

class AcceptanceCriterion
{
private:
    double *costPrevious_, *costCurrent_;

public:
    AcceptanceCriterion(double *costPrevious, double *costCurrent);
    int metropolisAC1();
    int metropolisAC3(double Th);
    int dCriteriaAC6();
    double p();
};

#endif
