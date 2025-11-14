#ifndef LENGTH_TEMPERATURE_HPP
#define LENGTH_TEMPERATURE_HPP

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <utils/SAParameters.hpp>

struct LengthParams{
    float len1;
    float len2;
    double len3;
    double len4;
    float len1_init;
    float len2_init;
    double len3_init;
    double len4_init;
};

class LengthTemperature
{
protected:
    SimulatedParams& saParams;
    LengthParams& ltParams;
public:

    LengthTemperature(SimulatedParams& saParams,LengthParams& ltParams)
    : saParams(saParams), ltParams(ltParams){};
    LengthParams& getLtParams() { return ltParams; };
    virtual bool apply() = 0;
    virtual ~LengthTemperature() = default;
};

class TL7: public LengthTemperature{
    public:
        TL7(SimulatedParams& saParams,LengthParams& ltParams) 
        :LengthTemperature(saParams,ltParams){};
        bool apply() override;
};
class TL8: public LengthTemperature{
    public:
        TL8(SimulatedParams& saParams,LengthParams& ltParams) 
        :LengthTemperature(saParams,ltParams){};
        bool apply() override;
};
class TL9: public LengthTemperature{
    public:
        TL9(SimulatedParams& saParams,LengthParams& ltParams) 
        :LengthTemperature(saParams,ltParams){};
        bool apply() override;
};
class TL11: public LengthTemperature{
    public:
        TL11(SimulatedParams& saParams,LengthParams& ltParams) 
        :LengthTemperature(saParams,ltParams){};
        bool apply() override;
};


#endif