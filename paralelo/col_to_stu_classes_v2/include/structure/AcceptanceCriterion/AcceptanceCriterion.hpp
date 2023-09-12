#ifndef ACCEPTANCE_CRITERION_HPP
#define ACCEPTANCE_CRITERION_HPP


#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <random>
#include <utils/SAParameters.hpp>

using std::mt19937;
using std::uniform_real_distribution;

extern mt19937 mt;
extern double temp;

struct AcceptanceParams{
    double Th;
};

class AcceptanceCriterion
{
protected:
    SimulatedParams& saParams;
    AcceptanceParams& acParams;

public:

    AcceptanceCriterion(SimulatedParams& saParams, AcceptanceParams& acParams)
    :saParams(saParams), acParams(acParams){};
    AcceptanceParams& getAcParams() { return acParams; };
    virtual ~AcceptanceCriterion() = default;
    virtual int apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta) = 0;
    double p(double costPrevious,double costCurrent, double temp);
};


class AC1: public AcceptanceCriterion{
    public:
        AC1(SimulatedParams& saParams, AcceptanceParams& acParams) 
        : AcceptanceCriterion(saParams, acParams) {}
        int apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta) override;
};

class AC3: public AcceptanceCriterion{
    public:
        AC3(SimulatedParams& saParams, AcceptanceParams& acParams) 
        : AcceptanceCriterion(saParams, acParams) {}
        int apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta) override;
};

class AC6: public AcceptanceCriterion{
    public:
        AC6(SimulatedParams& saParams, AcceptanceParams& acParams) 
        : AcceptanceCriterion(saParams, acParams) {}
        int apply(double &costPreviousSolution,double &costCurrentSolution, uniform_real_distribution<double> dist_accepta) override;
};



#endif