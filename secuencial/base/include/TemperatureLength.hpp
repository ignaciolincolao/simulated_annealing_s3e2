#ifndef TEMPERATURE_LENGTH_H
#define TEMPERATURE_LENGTH_H
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <math.h> 

using namespace std;

bool temperatureTL7(long double &temp,
    int &c_cooling_temperature,
    int &c_accepta,
    float len1,
    float len2,
    int n_colegios,
    long double coolingRate,
    int count);

bool temperatureTL8(long double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    float &len1,
    float len2,
    long double coolingRate);

bool temperatureTL9(long double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    long double &len3,
    long double &len4,
    long double coolingRate);

bool temperatureTL11(long double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    long double &len3,
    long double &len4,
    long double coolingRate);

#endif