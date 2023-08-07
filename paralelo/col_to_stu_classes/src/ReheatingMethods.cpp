#include <ReheatingMethods.hpp>

ReHeating::ReHeating(double *temp, double *k_reheating, int *n_reheating)
    : temp_{temp}, k_reheating_{k_reheating}, n_reheating_{n_reheating} {}

/**
 * TR11 Realiza un recalentamiento despues de no aceptar cualquier nueva solución en los ultimos k movimientos dado
 **/

void ReHeating::TR11(int count_rechaso)
{
    if (count_rechaso >= *n_reheating_)
    {
        *temp_ = *temp_ / *k_reheating_;
    }
}

/**
 * TR12 Realiza un recalentamiento despues de una cierta cantidad de movimientos
 **/

void ReHeating::TR12(int count)
{
    if (count % *n_reheating_ == 0)
    {
        *temp_ = *temp_ / *k_reheating_;
    }
}

/**
 * TR13 Realiza un recalentamiento despues de realizar n cantidad de enfrimiento
 **/

void ReHeating::TR13(int c_cooling_temperature)
{
    if (c_cooling_temperature % *n_reheating_ == 0)
    {
        *temp_ = *temp_ / *k_reheating_;
    }
}

/**
 * TR14 Aplica a una de los de recalentamiento convinando con T/K donde k se reduce por una constante E usando max(E,K-E)
 **/

void ReHeating::TR14(int k_reheating_init, int count_rechaso, double e_const)
{
    if (count_rechaso >= *n_reheating_)
    {
        *temp_ = *temp_ / *k_reheating_;
        *k_reheating_ = max(e_const, *k_reheating_ - e_const);
    }
    if (count_rechaso == 0)
    {
        *k_reheating_ = k_reheating_init;
    }
}
