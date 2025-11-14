#include <AcceptanceCriterion.hpp>



///////////////////////////////////////////////////
/// Metropolis-based
/// Función de aceptación en base a mayor temperatura mayor probabilidad que acepte a una solución peor
/// en caso de menor temperatura menor probabibilidad que acepte una solución peor.
///////////////////////////////////////////////////

int metropolisAC1(long double costPrevious, long double costCurrent){
    uniform_real_distribution<long double> dist_accepta(0.0, 1.0);
    if(costCurrent < costPrevious){
        return 1;
    }
    else{
        long double valor=p(costPrevious,costCurrent);
        long double nrandom=dist_accepta(mt);
        if(nrandom<valor){
            return 1;
        }
        else{
            return 0;
        }
    }
}
long double p(long double costPrevious,long double costCurrent){
    long double po;
    po = exp(-(costCurrent-costPrevious)/((long double)temp));
    return po;
}


/*
*  
*
*/
int metropolisAC3(long double costPrevious, long double costCurrent,long double Th){
    uniform_real_distribution<long double> dist_accepta(0.0, 1.0);
    //long double Th=1.2;
    if(costCurrent < costPrevious){
        return 1;
    }
    else{
        if(costPrevious <=costCurrent*Th){
            long double valor=p(costPrevious,costCurrent);
            long double nrandom=dist_accepta(mt);
            if(nrandom<valor){
                return 1;
            }
            else{
                return 0;
            }
        }
        else{
            return 0;
        }
    }
}

///////////////////////////////////////////////////
/// Deterministic criteria
///////////////////////////////////////////////////

int dCriteriaAC6(long double costPrevious, long double costCurrent){
    if(costCurrent - costPrevious <= temp){
        return 1;
    }
    else{
        return 0;
    }
}
