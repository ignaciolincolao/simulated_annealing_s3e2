#include <CoolingScheme.hpp>

CoolingScheme::CoolingScheme(double *temp, double coolingRate) : temp_{temp}, coolingRate_{coolingRate} {}

void CoolingScheme::CS2()
{
    *temp_ = (*temp_ * (coolingRate_));
}

double CoolingScheme::getTemp() { return *temp_; }
double CoolingScheme::getCoolingRate() { return coolingRate_; }