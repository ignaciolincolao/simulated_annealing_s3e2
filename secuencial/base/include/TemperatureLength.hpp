#ifndef TEMPERATURE_LENGTH_H
#define TEMPERATURE_LENGTH_H
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <math.h> 

using namespace std;

bool temperatureTL7(double &temp,
    int &c_cooling_temperature,
    int &c_accepta,
    float len1,
    float len2,
    int n_colegios,
    double coolingRate,
    int count);

bool temperatureTL8(double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    float &len1,
    float len2,
    double coolingRate);

bool temperatureTL9(double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    double &len3,
    double &len4,
    double coolingRate);

bool temperatureTL11(double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    double &len3,
    double &len4,
    double coolingRate);

#endif