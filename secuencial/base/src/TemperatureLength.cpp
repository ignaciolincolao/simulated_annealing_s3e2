#include <TemperatureLength.hpp>

/**
 * len1 =  length initial
 * len2 =  length const int
 * len3 =  length initial long double
 * len4 =  length const long double (less than 1)
**/
bool temperatureTL7(long double &temp,
    int &c_cooling_temperature,
    int &c_accepta,
    float len1,
    float len2,
    int n_colegios,
    long double coolingRate,
    int count)
{
    if(c_accepta>=len1*n_colegios){
        c_accepta=0;
        c_cooling_temperature++;
        return true;
    }
    if(count%(int(len2*n_colegios))==0){
        c_cooling_temperature++;
        return true;
    }
    return false;
}


/**
 * Arithmetic 
*/
bool temperatureTL8(long double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    float &len1,
    float len2,
    long double coolingRate)
{
    if(count_trials>=len1){
        len1 = len1+len2;
        count_trials=0;
        c_cooling_temperature++;
        return true;
    }  
    return false; 
}
/**
 * Geometric 
*/
bool temperatureTL9(long double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    long double &len3,
    long double &len4,
    long double coolingRate)
{
    if(count_trials>=len3){
        len3 = len3/len4;
        count_trials=0;
        c_cooling_temperature++;
        return true;
    }   
    return false;
}

/**
 * Exponential 
*/
bool temperatureTL11(long double &temp,
    int &c_cooling_temperature,
    int &count_trials,
    long double &len3,
    long double &len4,
    long double coolingRate)
{
    if(count_trials>=len3){
        len3 = pow(len3,1.0/len4);
        count_trials=0;
        c_cooling_temperature++;
        return true;
    }   
    return false;
}

