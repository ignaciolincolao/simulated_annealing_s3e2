#ifndef COOLING_SCHEME_H
#define COOLING_SCHEME_H

class CoolingScheme
{
private:
    double *temp_, coolingRate_;

public:
    CoolingScheme(double *temp, double coolingRate);
    void CS2();
    double getTemp();
    double getCoolingRate();
};

#endif