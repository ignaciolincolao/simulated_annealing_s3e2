#include <structure/CoolingScheme/CoolingScheme.hpp>

void CS2::apply()
{
    saParams.temp = (saParams.temp * (csParams.coolingRate));
}