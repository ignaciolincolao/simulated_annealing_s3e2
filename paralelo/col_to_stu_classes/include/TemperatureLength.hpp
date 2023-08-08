#ifndef TEMPERATURE_LENGTH_H
#define TEMPERATURE_LENGTH_H
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <math.h>

class TemperatureLength
{
public:
    bool TL7(int &c_cooling_temperature,
             int &c_accepta,
             float len1,
             float len2,
             int n_colegios,
             int count);

    bool TL8(int &c_cooling_temperature,
             int &count_trials,
             float &len1,
             float len2);

    bool TL9(int &c_cooling_temperature,
             int &count_trials,
             double &len3,
             double &len4);

    bool TL11(int &c_cooling_temperature,
              int &count_trials,
              double &len3,
              double &len4);
};

#endif