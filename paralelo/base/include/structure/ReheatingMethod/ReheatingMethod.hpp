#ifndef REHEATING_METHODS_HPP
#define REHEATING_METHODS_HPP

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <cstdlib>
#include <utils/SAParameters.hpp>

using std::cout;
using std::endl;
using std::max;

struct ReheatingParams
{
    double e_const;
    double max_temp;
    double k_reheating;
    int n_reheating;
    int k_reheating_init;
};


class ReheatingMethod
{
public:
    SimulatedParams& saParams;
    ReheatingParams& rmParams;
    ReheatingMethod(SimulatedParams& saParams,ReheatingParams& rmParams)
    : saParams(saParams),rmParams(rmParams){};
    ReheatingParams& getRmParams() { return rmParams; };
    virtual void apply() = 0;
    virtual ~ReheatingMethod() = default;
};
class TR0: public ReheatingMethod{
    public:
        TR0(SimulatedParams& saParams,ReheatingParams& rmParams):
        ReheatingMethod(saParams,rmParams){};
        void apply() override;
};

class TR11: public ReheatingMethod{
    public:
        TR11(SimulatedParams& saParams,ReheatingParams& rmParams):
        ReheatingMethod(saParams,rmParams){};
        void apply() override;
};
class TR12: public ReheatingMethod{
    public:
        TR12(SimulatedParams& saParams,ReheatingParams& rmParams):
        ReheatingMethod(saParams,rmParams){};
        void apply() override;
};
class TR13: public ReheatingMethod{
    public:
        TR13(SimulatedParams& saParams,ReheatingParams& rmParams):
        ReheatingMethod(saParams,rmParams){};
        void apply() override;
};
class TR14: public ReheatingMethod{
    public:
        TR14(SimulatedParams& saParams,ReheatingParams& rmParams):
        ReheatingMethod(saParams,rmParams){};
        void apply() override;
};

#endif