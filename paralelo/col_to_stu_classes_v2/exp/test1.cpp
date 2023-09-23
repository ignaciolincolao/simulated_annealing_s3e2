// please see the explanation in the documentation
// http://www.resibots.eu/limbo

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <SimulatedFactory.hpp>


// you can also include <limbo/limbo.hpp> but it will slow down the compilation
#include <limbo/bayes_opt/boptimizer.hpp>

using namespace limbo;

char timestr[20];
struct Params {
    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer {
        
    };

// depending on which internal optimizer we use, we need to import different parameters
#ifdef USE_NLOPT
    struct opt_nloptnograd : public defaults::opt_nloptnograd {
    };
#elif defined(USE_LIBCMAES)
    struct opt_cmaes : public defaults::opt_cmaes {
    };
#else
    struct opt_gridsearch : public defaults::opt_gridsearch {
    };
#endif

    struct kernel : public defaults::kernel {
        BO_PARAM(double, noise, 0.001);
    };

    struct bayes_opt_bobase : public defaults::bayes_opt_bobase {
        
    };

    struct kernel_maternfivehalves : public defaults::kernel_maternfivehalves {
    };

    struct init_randomsampling : public defaults::init_randomsampling {
        
    };

    struct stop_maxiterations : public defaults::stop_maxiterations {
        
    };

    // we use the default parameters for acqui_ucb
    struct acqui_ucb : public defaults::acqui_ucb {
    };
};

struct Eval {
    // number of input dimension (x.size())
    BO_PARAM(size_t, dim_in, 1);
    // number of dimensions of the result (res.size())
    BO_PARAM(size_t, dim_out, 1);

    // the function to be optimized
    Eigen::VectorXd operator()(const Eigen::VectorXd& x) const
    {
        double y = 0;

            // YOUR CODE HERE

            // return a 1-dimensional vector
        return tools::make_vector(y);
    }
};

int main(int argc, char *argv[])
{
    random_device rd;
    mt19937 mt(rd());

    // Hora Actual
    time_t hora_actual;
    struct tm *time_info;
    time(&hora_actual);
    time_info = localtime(&hora_actual);
    char timestr[20];
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);

    // Configuración del algoritmo
    RecordParams* rMgrParams = new RecordParams{
                .prefijo_save = string(timestr),
                .ruta_save = "../save/",
                .name_exp = "base"};

            
    SimulatedParams* saParams = new SimulatedParams{
        .seed = 123456,
        .n_students = 0,
        .n_colegios = 0,
        .count_rechaso = 0,
        .count = 0,
        .c_cooling_temperature = 0,
        .c_accepta = 0,
        .temp = 1.0,
        .min_temp = 0.0000009,
        .alpha1 = 15.0,
        .alpha2 = 30.0,
        .alpha3 = 25.0,
        .max_dist = 0.0,
        .min_dist = 0.0,
        .init_dist = 0.0,
        .costPrevious = 0.0,
        .costCurrent = 0.0};

    AcceptanceParams* acParams = new AcceptanceParams{
        .Th = 1.1};
    CoolingParams* csParams = new CoolingParams{
        .coolingRate = 0.99};
    LengthParams* ltParams = new LengthParams{
        .len1 = 1,
        .len2 = 2,
        .len3 = 1.0,
        .len4 = 0.999};
    ReheatingParams* rtParams = new ReheatingParams{
        .e_const = 0.01,
        .max_temp = std::numeric_limits<double>::max(),
        .k_reheating = 30,
        .n_reheating = 1,
        .k_reheating_init = 0};

    CUDAParams* cuParams = new CUDAParams{
        .n_block = 32,
        .n_thread = 32,
        .selectThread = 0,
        .selectBlock = 0};


    SimulatedStruct* simStruct = new SimulatedStruct{
        .acceptancecriterion = "AC1",
        .coolingscheme = "CS2",
        .lengthtemperature ="TL7",
        .reheatingmethod = "TR0"
    };

    SimulatedAnnealing *simulatedAnneling = SimulatedFactory::createSimulatedAnnealing(
            simStruct,
            rMgrParams,
            saParams,
            acParams,
            csParams,
            ltParams,
            rtParams,
            cuParams,
            mt);
    simulatedAnneling->runGPU();
    // we use the default acquisition function / model / stat / etc.
    bayes_opt::BOptimizer<Params> boptimizer;
    // run the evaluation
    boptimizer.optimize(Eval());
    // the best sample found
    std::cout << "Best sample: " << boptimizer.best_sample()(0) << " - Best observation: " << boptimizer.best_observation()(0) << std::endl;
    return 0;
}
