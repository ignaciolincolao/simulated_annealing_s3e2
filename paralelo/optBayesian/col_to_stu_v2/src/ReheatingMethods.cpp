#include <ReheatingMethods.hpp>




/**
 * TR11 Realiza un recalentamiento despues de no aceptar cualquier nueva solución en los ultimos k movimientos dado
**/

void reheatingTR11(double &temp, 
    double k_reheating,
    int n_reheating,
    int count_rechaso)
{
    if(count_rechaso >= n_reheating){
        temp= temp/k_reheating;
        cout << "reheating" << temp << endl;
    }
}


/**
 * TR12 Realiza un recalentamiento despues de una cierta cantidad de movimientos
**/

void reheatingTR12(double &temp, 
    double k_reheating,
    int n_reheating,
    int count)
{
    if(count % n_reheating == 0){
        temp= temp/k_reheating;
    }
}

/**
 * TR13 Realiza un recalentamiento despues de realizar n cantidad de enfrimiento
**/

void reheatingTR13(double &temp,
    double k_reheating,
    int n_reheating,
    int c_cooling_temperature)
{
    if(c_cooling_temperature % n_reheating == 0){
        temp= temp/k_reheating;
    }
}

/**
 * TR14 Aplica a una de los de recalentamiento convinando con T/K donde k se reduce por una constante E usando max(E,K-E)
**/

void reheatingTR14(double &temp,
    double &k_reheating,
    int k_reheating_init,
    int n_reheating,
    int count_rechaso,
    double e_const)
{
    if(count_rechaso >= n_reheating){
        temp= temp/k_reheating;
        k_reheating = max(e_const, k_reheating-e_const);
    }
    if(count_rechaso==0){
        k_reheating = k_reheating_init;
    }
    
}
