#ifndef REHEATING_METHODS_H
#define REHEATING_METHODS_H

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::max;

class ReHeating
{
private:
    double *temp_, *k_reheating_;
    int *n_reheating_;

public:
    ReHeating(double *temp, double *k_reheating, int *n_reheating);
    void TR11(int count_rechaso);

    void TR12(int count);

    void TR13(int c_cooling_temperature);

    void TR14(int k_reheating_init, int count_rechaso, double e_const);
};

#endif