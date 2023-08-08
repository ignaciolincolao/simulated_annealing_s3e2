#include <AcceptanceCriterion.hpp>

///////////////////////////////////////////////////
/// Metropolis-based
/// Función de aceptación en base a mayor temperatura mayor probabilidad que acepte a una solución peor
/// en caso de menor temperatura menor probabibilidad que acepte una solución peor.
///////////////////////////////////////////////////

AcceptanceCriterion::AcceptanceCriterion(double *costPrevious, double *costCurrent)
    : costPrevious_{costPrevious}, costCurrent_{costCurrent} {}

int AcceptanceCriterion::metropolisAC1()
{
    uniform_real_distribution<double> dist_accepta(0.0, 1.0);
    if (*costCurrent_ < *costPrevious_)
    {
        return 1;
    }
    else
    {
        double valor = p();
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

double AcceptanceCriterion::p()
{
    return exp(-(*costCurrent_ - *costPrevious_) / ((double)temp));
}

/*
 *
 *
 */
int AcceptanceCriterion::metropolisAC3(double Th)
{
    uniform_real_distribution<double> dist_accepta(0.0, 1.0);
    // double Th=1.2;
    if (*costCurrent_ < *costPrevious_)
    {
        return 1;
    }
    else
    {
        if (*costPrevious_ <= *costCurrent_ * Th)
        {
            double valor = p();
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

int AcceptanceCriterion::dCriteriaAC6()
{
    if (*costCurrent_ - *costPrevious_ <= temp)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
