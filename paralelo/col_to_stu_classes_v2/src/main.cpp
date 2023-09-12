#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <SimulatedFactory.hpp>





int main(int argc, char *argv[])
{


    SimulatedAnnealing *simulatedAnneling = SimulatedFactory::createSimulatedAnnealing(
            "AC1",
            "CS2",
            "TL7",
            "TR11",
            argc,
            argv);
    cout << "test" << endl;
    simulatedAnneling->runGPU();
    return (EXIT_SUCCESS);
}
