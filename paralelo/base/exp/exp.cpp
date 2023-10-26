// please see the explanation in the documentation
// http://www.resibots.eu/limbo
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <SimulatedFactory.hpp>
#include <fstream>

#include <limbo/limbo.hpp>
#include <fstream>
#include <limbo/kernel/exp.hpp>
#include <limbo/kernel/squared_exp_ard.hpp>
#include <limbo/mean/data.hpp>
#include <limbo/model/gp.hpp>
#include <limbo/model/gp/kernel_lf_opt.hpp>
#include <limbo/tools.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/serialize/binary_archive.hpp>
#include <limbo/serialize/text_archive.hpp>

using namespace limbo;

template <typename Params>
double valRange(Params a, Params b,const Eigen::VectorXd& x,int n_id){
    return a + x[n_id] * ( b - a );
}

struct Params {
    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer {
        BO_PARAM(int, hp_period, 10);
    };
    struct bayes_opt_bobase : public defaults::bayes_opt_bobase {
        BO_PARAM(int, stats_enabled, true);
    };
    struct stop_maxiterations {
        BO_PARAM(int, iterations, 50);
    };
    struct acqui_ei {
        BO_PARAM(double, jitter, 0.0);
    };
    struct init_randomsampling {
        BO_PARAM(int, samples, 10);
    };
    struct kernel : public defaults::kernel {
        BO_PARAM(double, noise, 1e-10);
        BO_PARAM(bool, optimize_noise, true);
    };
    struct kernel_squared_exp_ard : public defaults::kernel_squared_exp_ard {
    };
    struct opt_rprop : public defaults::opt_rprop {
    };
    struct opt_parallelrepeater : public defaults::opt_parallelrepeater {
    };
    struct opt_cmaes : public defaults::opt_cmaes {
        BO_PARAM(int, restarts, 1);
        BO_PARAM(int, max_fun_evals, -1);
    };
};

Eigen::VectorXd simAnnealing(const Eigen::VectorXd& x)
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

    double alp1 =  15.0;
    double alp2 = 30.0;
    double alp3 = 25.0;
    SimulatedParams* saParams = new SimulatedParams{
        .seed = static_cast<int>(time(nullptr)),
        .n_students = 0,
        .n_colegios = 0,
        .count_rechaso = 0,
        .count = 0,
        .c_cooling_temperature = 0,
        .c_accepta = 0,
        .temp = 1.0,
        .min_temp = 0.0000009,
        .alpha1 = alp1,
        .alpha2 = alp2,
        .alpha3 = alp3,
        .max_dist = 0.0,
        .min_dist = 0.0,
        .init_dist = 0.0,
        .costPrevious = 0.0,
        .costCurrent = 0.0,
        .alpha = {alp1, alp2, alp3}};

    AcceptanceParams* acParams = new AcceptanceParams{
        .Th = 1.1};
    CoolingParams* csParams = new CoolingParams{
        .coolingRate = 0.97};
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


    /*
    Realiza los cambios de las variables
    */

    cuParams->n_block = 32*int(valRange(1,32,x,4));
    cuParams->n_thread = 32*int(valRange(1,32,x,5));
    saParams->temp = valRange(100,10000,x,0);
    csParams->coolingRate = valRange(0.9,0.999,x,1);//valRange(0.9,0.99,x,1);//
    ltParams->len1 = valRange(1,100,x,2);//valRange(1,10,x,2);//
    ltParams->len2 = valRange(1,100,x,3);//valRange(1,10,x,3);//
    cout << "valores de x: " << "temp= "<< saParams->temp 
        <<" | coolingRate= " << csParams->coolingRate 
        << " | len1= " << ltParams->len1
        << " | len2= " << ltParams->len2
        << " | n_block= " << cuParams->n_block
        << " | n_thread= " << cuParams->n_thread
        << endl;

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

    double val = simulatedAnneling->runGPU();
    Eigen::VectorXd res(2);
    res[0] = val;
    res[1] = simulatedAnneling->saParams.count;
    delete simulatedAnneling;
    delete simStruct;
    delete rMgrParams;
    delete saParams;
    delete acParams;
    delete csParams;
    delete ltParams;
    delete rtParams;
    delete cuParams;

    return res;
}


template <typename Params>
struct LinearWeighting {
    using result_type = double;
    LinearWeighting(const int& weight) : _weight(weight){}

    double operator()(const Eigen::VectorXd& x) const
    {
        return -(1*x[0]+(x[1]/_weight));
        //return -(1*x[0]+0.000001*x[1]);
    }
    protected:
        int _weight;

};

template <typename Params>
struct eval_func {
    BO_PARAM(size_t, dim_in, 6);
    BO_PARAM(size_t, dim_out, 2);

    eval_func() {}

    Eigen::VectorXd operator()(const Eigen::VectorXd& x) const
    {
        Eigen::VectorXd values = simAnnealing(x);
        return values;
    }
};

int main()
{
    using kernel_t = kernel::SquaredExpARD<Params>;

    using mean_t = mean::Data<Params>;

    using gp_opt_t = model::gp::KernelLFOpt<Params>;
    using gp_t = model::GP<Params, kernel_t, mean_t,gp_opt_t>;

    using acqui_t = acqui::EI<Params, gp_t>;
    using acqui_opt_t = opt::Cmaes<Params>;

    using init_t = init::RandomSampling<Params>;

    using stop_t = boost::fusion::vector<stop::MaxIterations<Params>>;

    using stat_t = boost::fusion::vector<stat::ConsoleSummary<Params>, stat::Samples<Params>, stat::Observations<Params>, stat::AggregatedObservations<Params>, stat::GPAcquisitions<Params>, stat::BestAggregatedObservations<Params>, stat::GPKernelHParams<Params>>;

    bayes_opt::BOptimizer<Params, modelfun<gp_t>, acquifun<acqui_t>, acquiopt<acqui_opt_t>, initfun<init_t>, statsfun<stat_t>> boptimizer;
    // Instantiate aggregator

    // Optimization
    /*
    LinearWeighting<Params> aggregator(1000000);
    boptimizer.optimize(eval_func<Params>(), aggregator);
    std::cout << "New target!" << std::endl;
    aggregator = LinearWeighting<Params>(10000000);
    boptimizer.optimize(eval_func<Params>(), aggregator, false);
    */
    // Do not forget to pass `false` as the last parameter in `optimize`,

    

    gp_t gp_ard(6, 2);
    limbo::serialize::TextArchive f1("base_val");
    gp_ard.load(f1);

    //boptimizer.optimize_hyperparams();
    //std::cout << "Optimizando Hyperparametros" << std::endl;

    
    return 1;
}
