#ifndef REHEATING_METHODS_H
#define REHEATING_METHODS_H

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
using namespace std;


void reheatingTR11(long double &temp, 
    long double k_reheating,
    int n_reheating,
    int count_rechaso);
void reheatingTR12(long double &temp, 
    long double k_reheating,
    int n_reheating,
    int count);
void reheatingTR13(long double &temp,
    long double k_reheating,
    int n_reheating,
    int c_cooling_temperature);
void reheatingTR14(long double &temp,
    long double k_reheating,
    int k_reheating_init,
    int n_reheating,
    int count_rechaso,
    long double e_const);
#endif