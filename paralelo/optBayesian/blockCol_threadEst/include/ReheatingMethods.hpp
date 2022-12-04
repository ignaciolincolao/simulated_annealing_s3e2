#ifndef REHEATING_METHODS_H
#define REHEATING_METHODS_H

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
using namespace std;


void reheatingTR11(double &temp,
    int &count_reheating,
    double k_reheating,
    int n_reheating,
    int count_rechaso,
    int max_reheating);
void reheatingTR12(double &temp, 
    double k_reheating,
    int n_reheating,
    int count);
void reheatingTR13(double &temp,
    double k_reheating,
    int n_reheating,
    int c_cooling_temperature);
void reheatingTR14(double &temp,
    double k_reheating,
    int k_reheating_init,
    int n_reheating,
    int count_rechaso,
    double e_const);
#endif