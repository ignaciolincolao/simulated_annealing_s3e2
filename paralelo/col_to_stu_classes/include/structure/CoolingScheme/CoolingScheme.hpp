#ifndef COOLING_SCHEME_HPP
#define COOLING_SCHEME_HPP
#include <utils/SAParameters.hpp>

struct CoolingParams{
    double coolingRate;
};

class CoolingScheme
{

protected:
    SimulatedParams& saParams;
    CoolingParams& csParams;
public:
   
    CoolingScheme(SimulatedParams& saParams,CoolingParams& csParams)
    : saParams(saParams), csParams(csParams){};
    CoolingParams& getCsParams() { return csParams; };
    virtual ~CoolingScheme() = default;
    virtual void apply() = 0;
};

class CS2: public CoolingScheme{
    public: 
        CS2(SimulatedParams& saParams,CoolingParams& csParams)
        : CoolingScheme(saParams, csParams){};
        void apply() override;
};

#endif