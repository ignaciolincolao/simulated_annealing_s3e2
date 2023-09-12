#include <structure/LengthTemperature/LengthTemperature.hpp>

/**
 * len1 =  length initial
 * len2 =  length const int
 * len3 =  length initial double
 * len4 =  length const double (less than 1)
 **/


bool TL7::apply()
{
    if (saParams.c_accepta >= ltParams.len1 * saParams.n_colegios)
    {
        saParams.c_accepta = 0;
        saParams.c_cooling_temperature++;
        return true;
    }
    if (saParams.count % (int(ltParams.len2 * saParams.n_colegios)) == 0)
    {
        saParams.c_cooling_temperature++;
        return true;
    }
    return false;
}

/**
 * Arithmetic
 */
bool TL8::apply()
{
    if (saParams.count_trials >= ltParams.len1)
    {
        ltParams.len1 = ltParams.len1 + ltParams.len2;
        saParams.count_trials = 0;
        saParams.c_cooling_temperature++;
        return true;
    }
    return false;
}


bool TL9::apply()
{
    if (saParams.count_trials >= ltParams.len3)
    {
        ltParams.len3 = ltParams.len3 / ltParams.len4;
        saParams.count_trials = 0;
        saParams.c_cooling_temperature++;
        return true;
    }
    return false;
}

/**
 * Exponential
 */
bool TL11::apply()
{
    if (saParams.count_trials >= ltParams.len3)
    {
        ltParams.len3 = pow(ltParams.len3, 1.0 / ltParams.len4);
        saParams.count_trials = 0;
        saParams.c_cooling_temperature++;
        return true;
    }
    return false;
}
