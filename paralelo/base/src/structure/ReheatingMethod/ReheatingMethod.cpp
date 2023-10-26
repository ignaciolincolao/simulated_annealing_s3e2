#include <structure/ReheatingMethod/ReheatingMethod.hpp>

void TR0::apply()
{
    
}


/**
 * TR11 Realiza un recalentamiento despues de no aceptar cualquier nueva solución en los ultimos k movimientos dado
 **/

void TR11::apply()
{
    if (saParams.count_rechaso >= rmParams.n_reheating)
    {
        saParams.temp = saParams.temp / rmParams.k_reheating;
    }
}

/**
 * TR12 Realiza un recalentamiento despues de una cierta cantidad de movimientos
 **/

void TR12::apply()
{
    if (saParams.count % rmParams.n_reheating == 0)
    {
        saParams.temp = saParams.temp / rmParams.k_reheating;
    }
}

/**
 * TR13 Realiza un recalentamiento despues de realizar n cantidad de enfrimiento
 **/

void TR13::apply()
{
    if (saParams.c_cooling_temperature % rmParams.n_reheating == 0)
    {
        saParams.temp = saParams.temp / rmParams.k_reheating;
    }
}

/**
 * TR14 Aplica a una de los de recalentamiento convinando con T/K donde k se reduce por una constante E usando max(E,K-E)
 **/

void TR14::apply()
{
    if (saParams.count_rechaso >= rmParams.n_reheating)
    {
        saParams.temp = saParams.temp / rmParams.k_reheating;
        rmParams.k_reheating = max(rmParams.e_const, rmParams.k_reheating - rmParams.e_const);
    }
    if (saParams.count_rechaso == 0)
    {
        rmParams.k_reheating = rmParams.k_reheating_init;
    }
}
