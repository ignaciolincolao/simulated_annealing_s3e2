#include <structure/AcceptanceCriterion/AcceptanceCriterion.hpp>


///////////////////////////////////////////////////
/// Metropolis-based
/// Función de aceptación en base a mayor temperatura mayor probabilidad que acepte a una solución peor
/// en caso de menor temperatura menor probabibilidad que acepte una solución peor.
///////////////////////////////////////////////////

double AcceptanceCriterion::p(double costPrevious,double costCurrent,double temp){
    double po;
    po = exp(-(costCurrent-costPrevious)/((double)temp));
    return po;
}

int AC1::apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta){
    if (costCurrentSolution < costPreviousSolution)
    {
        return 1;
    }
    else
    {
        double valor = p(costPreviousSolution,costCurrentSolution,saParams.temp);
        double nrandom = dist_accepta(mt);
        if (nrandom < valor)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

int AC3::apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta) {
    // double Th=1.2;
    if (costCurrentSolution < costPreviousSolution)
    {
        return 1;
    }
    else
    {
        if (costPreviousSolution <= costCurrentSolution * acParams.Th)
        {
        double valor = p(costPreviousSolution,costCurrentSolution,saParams.temp);
        double nrandom = dist_accepta(mt);
            if (nrandom < valor)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return 0;
        }
    }
}
///////////////////////////////////////////////////
/// Deterministic criteria
///////////////////////////////////////////////////

int AC6::apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta) {
    if (costCurrentSolution - costPreviousSolution <= saParams.temp)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}