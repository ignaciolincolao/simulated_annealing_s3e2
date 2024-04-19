// please see the explanation in the documentation
// http://www.resibots.eu/limbo
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdio.h>
#include <SimulatedFactory.hpp>
#include <fstream>
#include <chrono>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <limbo/limbo.hpp>
#include <fstream>
#include <limits>

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
int seed = 0;
int n_block_min = 1;
int n_block_max = 32;
int n_block_factor = 32;
int n_block_idX = 3;//int n_block_idX = 4;
int n_thread_min = 1;
int n_thread_max = 32;
int n_thread_factor = 32;
int n_thread_idX;
double temp_min = 1000;//100;
double temp_max = 100000;//10000;
int temp_idX = 4;
double coolingRate_min = 0.9;
double coolingRate_max = 0.999;
int coolingRate_idX = 0;//int coolingRate_idX = 1;
float len1_min = 1.f;
float len1_max = 10.f;
int len1_idX = 1;//int len1_idX = 2;
float len2_min = 1.f;
float len2_max = 10.f;
int len2_idX = 2;//int len2_idX = 3;
int n_block;
int n_thread;
int it;
double temp;
double coolingRate;
float len1;
float len2;
/*
    cuParams->n_block = n_block;//32*int(valRange(1,32,x,4));
    cuParams->n_thread = n_thread;//32*int(valRange(1,32,x,5));
    saParams->temp = temp;
    csParams->coolingRate = coolingRate;//valRange(0.9,0.999,x,1);//valRange(0.9,0.99,x,1);//
    ltParams->len1 = len1;//valRange(1,100,x,2);////
    ltParams->len2 = len2;//valRange(1,100,x,3);//valRange(1,10,x,3);//

*/



template <typename Params>
double valRange(Params a, Params b,const Eigen::VectorXd& x,int n_id){
    return a + x[n_id] * ( b - a );
}

struct Params {
    struct bayes_opt_boptimizer : public defaults::bayes_opt_boptimizer {
        BO_PARAM(int, hp_period, 1);
    };
    struct bayes_opt_bobase : public defaults::bayes_opt_bobase {
        BO_PARAM(int, stats_enabled, true);
    };
    struct stop_maxiterations {
        BO_PARAM(int, iterations, 1000);
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
    seed = mt();
    SimulatedParams* saParams = new SimulatedParams{
        .seed = seed,
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

    int block_threads[26][2] = {
                                {1,32},
                                {2,32},
                                {1,64},
                                {8,32},
                                {4,64},
                                {2,128},
                                {1,256},
                                {8,64},
                                {4,128},
                                {2,256},
                                {1,512},
                                {32,32},
                                {8,128},
                                {4,256},
                                {2,512},
                                {1,1024},
                                {64,32},
                                {32,64},
                                {8,256},
                                {4,512},
                                {2,1024},
                                {128,32},
                                {64,64},
                                {32,128},
                                {8,512},
                                {4,1024}
                            };

    /*
    Realiza los cambios de las variables
    */
    //n_block = pow(2,floor(x[n_block_idX]*7-std::numeric_limits<float>::epsilon()) +1);//n_block_factor*int(valRange(n_block_min,n_block_max,x,n_block_idX));
    //n_thread= pow(2,floor(x[n_thread_idX]*6-std::numeric_limits<float>::epsilon()) +5);//n_thread_factor*int(valRange(n_thread_min,n_thread_max,x,n_thread_idX));
    int indice= max(int(floor(x[n_block_idX]*26-std::numeric_limits<float>::epsilon())),0);
    //cout << floor(1.00*13-std::numeric_limits<float>::epsilon()) << endl;
    n_block = block_threads[indice][0];
    n_thread = block_threads[indice][1];
    temp= valRange(temp_min,temp_max,x,temp_idX);//pow(n_block*n_thread,1.5);
    coolingRate= valRange(coolingRate_min,coolingRate_max,x,coolingRate_idX);
    len1= valRange(len1_min,len1_max,x,len1_idX);
    len2= valRange(len2_min,len2_max,x,len2_idX);
    cuParams->n_block = n_block;//32*int(valRange(1,32,x,4));
    cuParams->n_thread = n_thread;//32*int(valRange(1,32,x,5));
    saParams->temp = temp;
    csParams->coolingRate = coolingRate;//valRange(0.9,0.999,x,1);//valRange(0.9,0.99,x,1);//
    ltParams->len1 = len1;//valRange(1,100,x,2);////
    ltParams->len2 = len2;//valRange(1,100,x,3);//valRange(1,10,x,3);//
    cout << "valores de x: " << "temp (No cuenta como x)= "<< saParams->temp 
        <<" | coolingRate= " << csParams->coolingRate 
        << " | len1= " << ltParams->len1
        << " | len2= " << ltParams->len2
        << " | n_block= " << cuParams->n_block
        << " | n_thread= " << cuParams->n_thread
        << endl;

    cout << "valores de x: "  
        <<" coolingRate= " << x[coolingRate_idX]
        << " | len1= " << x[len1_idX]
        << " | len2= " << x[len2_idX]
        << " | n_block_thread= " << x[n_block_idX]
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
    Eigen::VectorXd res(1);
    res[0] = val;
    it = simulatedAnneling->saParams.count;
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
        return -(x[0]);
        //return -(1*x[0]+0.000001*x[1]);
    }
    protected:
        int _weight;

};

template <typename Params>
struct eval_func {
    BO_PARAM(size_t, dim_in, 5);
    BO_PARAM(size_t, dim_out, 1);

    eval_func() {}

    Eigen::VectorXd operator()(const Eigen::VectorXd& x) const
    {
        Eigen::VectorXd values = simAnnealing(x);
        return values;
    }
};


template <typename Params>
struct SeedObservation : public stat::StatBase<Params> {
    template <typename BO, typename AggregatorFunction>
    void operator()(const BO& bo, const AggregatorFunction& afun)
    {
        // [optional] if statistics have been disabled or if there are no observations, we do not do anything
        if (!bo.stats_enabled() || bo.observations().empty())
            return;

        // [optional] we create a file to write / you can use your own file but remember that this method is called at each iteration (you need to create it in the constructor)
        this->_create_log_file(bo, "seed_iteration.dat");

        if (bo.total_iterations() == 0)
            (*this->_log_file) << "id,"
                                << "seed,"
                                <<  "x0,"
                                <<  "x1," 
                                <<  "x2," 
                                <<  "x3," 
                                <<  "x4,"
                                <<  "lz,"
                                <<  "z,"
                                <<  "it,"
                                <<  "n_block,"
                                <<  "n_thread,"
                                <<  "temp,"
                                <<  "coolingRate,"
                                <<  "len1,"
                                <<  "len2"
                                << std::endl;

        // ----- write what we have found ------
        // the file is (*this->_log_file)
        auto samples = bo.samples().back().transpose();
        auto observations = bo.observations().back().transpose();
        auto lZ = afun(bo.observations().back());


        (*this->_log_file) << bo.total_iterations() << "," 
                           << seed << ",";
        for (auto& sample : samples){
            (*this->_log_file) << std::fixed << std::setprecision(15) << sample << ",";
        }
        (*this->_log_file) << lZ << ",";
        for (auto& observation : observations){
            (*this->_log_file) << std::defaultfloat << std::setprecision(6) << observation << ",";
        }
        (*this->_log_file)  << it << ","
                            << n_block << ","
                            << n_thread << ","
                            << temp << ","
                            << coolingRate << ","
                            << len1 << ","
                            << len2;
        (*this->_log_file) << endl;
                           
                           
                           
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

    using stat_t = boost::fusion::vector<stat::ConsoleSummary<Params>, 
                    stat::Samples<Params>, stat::Observations<Params>,
                    stat::AggregatedObservations<Params>,
                    stat::GPAcquisitions<Params>,
                    stat::BestAggregatedObservations<Params>,
                    stat::GPKernelHParams<Params>,
                    SeedObservation<Params>
                    >;

    bayes_opt::BOptimizer<Params, modelfun<gp_t>, acquifun<acqui_t>, acquiopt<acqui_opt_t>, initfun<init_t>, statsfun<stat_t>> boptimizer;
    // Instantiate aggregator

    // Optimization
    
    LinearWeighting<Params> aggregator(10000000);
    boptimizer.optimize(eval_func<Params>(), aggregator);
    //std::cout << "New target!" << std::endl;
    //aggregator = LinearWeighting<Params>(100000000);
    //boptimizer.optimize(eval_func<Params>(), aggregator, false);
    
    // Do not forget to pass `false` as the last parameter in `optimize`,

    
    /*
    gp_t gp_ard(6, 2);

    limbo::serialize::TextArchive f1(directoryName);
    gp_ard.load(f1);
    */
    //boptimizer.optimize_hyperparams();
    //std::cout << "Optimizando Hyperparametros" << std::endl;
    
    return 1;
}
