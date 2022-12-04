#include <AcceptanceCriterion.hpp>



///////////////////////////////////////////////////
/// Metropolis-based
/// Función de aceptación en base a mayor temperatura mayor probabilidad que acepte a una solución peor
/// en caso de menor temperatura menor probabibilidad que acepte una solución peor.
///////////////////////////////////////////////////

int metropolisAC1(double costPrevious, double costCurrent){
    uniform_real_distribution<double> dist_accepta(0.0, 1.0);
    if(costCurrent < costPrevious){
        return 1;
    }
    else{
        double valor=p(costPrevious,costCurrent);
        double nrandom=dist_accepta(mt);
        if(nrandom<valor){
            return 1;
        }
        else{
            return 0;
        }
    }
}
double p(double costPrevious,double costCurrent){
    double po;
    po = exp(-(costCurrent-costPrevious)/((double)temp));
    return po;
}


/*
*  
*
*/
int metropolisAC3(double costPrevious, double costCurrent,double Th){
    uniform_real_distribution<double> dist_accepta(0.0, 1.0);
    //double Th=1.2;
    if(costCurrent < costPrevious){
        return 1;
    }
    else{
        if(costPrevious <=costCurrent*Th){
            double valor=p(costPrevious,costCurrent);
            double nrandom=dist_accepta(mt);
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

int dCriteriaAC6(double costPrevious, double costCurrent){
    if(costCurrent - costPrevious <= temp){
        return 1;
    }
    else{
        return 0;
    }
}
