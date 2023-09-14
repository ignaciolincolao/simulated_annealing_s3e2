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

    random_device rd;
    mt19937 mt(rd());
    SimulatedAnnealing *simulatedAnneling = SimulatedFactory::createSimulatedAnnealing(
            "AC1",
            "CS2",
            "TL7",
            "TR0",
            argc,
            argv,
            mt);
    simulatedAnneling->runGPU();
    return (EXIT_SUCCESS);
}
