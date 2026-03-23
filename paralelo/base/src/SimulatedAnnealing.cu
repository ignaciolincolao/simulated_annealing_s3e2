#include <SimulatedAnnealing.cuh>
#include <CUDAWrapper.cuh>

#include <limits>
#include <iomanip> //para imprimir mas decimales

//borrar
#include <fstream>
#include <filesystem>   // C++17 para verificar existencia del archivo

#ifndef SAVE_DATA
#define SAVE_DATA 0
#endif
#define DECIMAL 16

typedef std::numeric_limits<double> dbl;


SimulatedAnnealing::SimulatedAnnealing(AcceptanceCriterion* AC,
    CoolingScheme* CS,
    LengthTemperature* LT,
    ReheatingMethod* RM,
    Dataset* DS,
    RecordManager* RMgr,
    SimulatedParams* saParams,
    CUDAParams* cuParams,
    mt19937& mt)
    : 
    acceptanceCriterion(AC), 
    coolingScheme(CS),
    lengthTemperature(LT),
    reheatingMethod(RM),
    dataSet(DS),
    recordManager(RMgr),
    saParams(*saParams),
    cuParams(*cuParams),
    acParams(AC->getAcParams()),
    csParams(CS->getCsParams()),
    ltParams(LT->getLtParams()),
    rmParams(RM->getRmParams()),
    rmgrParams(RMgr->getRmgrParams()),
    mt(mt), 
    dist(0, 0), 
    dist2(0, 0), 
    dist_accepta(0.0, 1.0)
    {      
        mt.seed(saParams->seed);
        probSelection = (double *)malloc(sizeof(double)*cuParams->n_block*cuParams->n_thread);
        UpdateProb(0);
    }
SimulatedAnnealing::~SimulatedAnnealing(){
    delete acceptanceCriterion;
    delete coolingScheme;
    delete lengthTemperature;
    delete reheatingMethod;
    delete dataSet;
    delete recordManager;

    cudaFreeHost(previousSolution);
    cudaFreeHost(bestSolution);
    cudaFreeHost(currentSolution);
    cudaFreeHost(cupoArray);
    cudaFreeHost(alumnosSep);
    free(aluxcol);
    free(aluVulxCol);
    free(previousAluxCol);
    free(previousAluVulxCol);
    free(bestAluxCol);
    free(bestAluVulxCol);
    cudaFreeHost(currentVars);
    cudaFreeHost(previousVars);
    cudaFreeHost(bestVars);
    cudaFreeHost(matrestest);
    for(int i=0; i < saParams.n_students; i++ ){
        free(distMat[i]);
    }
    free(distMat);
    cudaFreeHost(saParams.shuffle_student);
    cudaFreeHost(saParams.shuffle_colegios);
}   

double SimulatedAnnealing::runGPU(){
    CUDAWrapper* cudaWrapper = new CUDAWrapper(cuParams, saParams, mt);
    // cout << "test" << endl;
    inicializationValues(cudaWrapper);
    cudaWrapper->memInit(previousSolution,
        bestSolution,
        currentSolution,
        cupoArray,
        alumnosSep,
        totalVuln,
        aluxcol,
        aluVulxCol,
        matrestest,
        alpha,
        currentVars);

    cout << "--------------- Primeros datos -------------\n";
    cout << "Primer costo de solución: " << costBestSolution << "\n";
    cout << "Primer distancia: " << meanDist(currentSolution, distMat)/saParams.max_dist << "\n"; //lo normalize
    cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    cout << "Primer CostoCupo: " << costCupo(currentSolution, cupoArray) << "\n";
    cout << "Penalty inicial: " << penaltyParents(currentSolution,h_penalty_matrix)/saParams.n_students << "\n\n";
#if SAVE_DATA
    #if ENABLE_OPEN_RECORD_INFO
    recordManager->openRecordInfo();
    recordManager->SaveInfoInit(costBestSolution,
        meanDist(currentSolution, distMat),
        S(currentSolution, alumnosSep, totalVuln),
        costCupo(currentSolution, cupoArray),
        penaltyParents(currentSolution,h_penalty_matrix));
    recordManager->closeRecordInfo();
    #endif
    #ifdef ENABLE_OPEN_RECORD_GRAPHICS
    recordManager->openRecordGraphics();
    recordManager->SaveGraphicsInit(meanDist(currentSolution, distMat),
    S(currentSolution, alumnosSep, totalVuln),
    costCupo(currentSolution, cupoArray),
    costCurrentSolution,
    penaltyParents(currentSolution,h_penalty_matrix));
    recordManager->closeRecordGraphics();
    #endif

    #ifdef ENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION
    recordManager->openRecordGraphicsBestSolution();
    recordManager->SaveGraphicsBestSolution(currentSolution);
    recordManager->closeRecordGraphicsBestSolution();
    #endif 
    
    #ifdef ENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION_RBD
    recordManager->openRecordGraphicsBestSolutionRBD();
    recordManager->SaveGraphicsFirstSolutionRBD(currentSolution, dataSet->ptr_colegios, dataSet->ptr_students);
    recordManager->closeRecordGraphicsBestSolutionRBD();
    #endif
#endif
    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////
    int id_select;
    auto start = std::chrono::high_resolution_clock::now();
    ///////////////////////////////////////////////////
    /// Comienza a ejecutarse el algoritmo de SA
    ///////////////////////////////////////////////////

    #if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_REGISTER
    float diff_temp = saParams.temp - saParams.min_temp;
    for(int i=0; i < recordManager->vector_percentage.size();i++){
        recordManager->vector_temp_percentage[i] = diff_temp * (recordManager->vector_percentage[i]/100) + saParams.min_temp;
    }
    #endif

    #endif
    saParams.count++;

    #ifdef ENABLE_GPU_RECORD_TIME
    std::ofstream tiempoGPU;
    tiempoGPU.open("../../save/tiempo_GPU.txt", std::ios::app);
    #endif


    while(saParams.temp > saParams.min_temp){
        #ifdef ENABLE_GPU_RECORD_TIME
        auto start_shuffle = std::chrono::high_resolution_clock::now();
        #endif
        ///////////////////////////////////////////////////
        /// Copia Solución Anterior a la actual
        ///////////////////////////////////////////////////
        cudaWrapper->memCopyPrevToCurrent();
        ///////////////////////////////////////////////////
        ///  Selecciona aleatoria mente a los alumnos
        ///////////////////////////////////////////////////
        shuffle(saParams.shuffle_student, saParams.max_changes_students, dist);
        shuffle(saParams.shuffle_colegios, saParams.max_changes_school, dist2);

        #ifdef ENABLE_GPU_RECORD_TIME
        auto end_shuffle = std::chrono::high_resolution_clock::now();
        double time_taken_shuffle = std::chrono::duration_cast<std::chrono::nanoseconds>(end_shuffle - start_shuffle).count();
        tiempoGPU << std::setprecision(10)<<time_taken_shuffle << ",";
        tiempoGPU.flush();
        #endif
        ///////////////////////////////////////////////////
        ///  Envia datos a GPU
        ///////////////////////////////////////////////////
    
        cudaWrapper->uploadCurrentMemorySolution();
        ///////////////////////////////////////////////////
        ///  Ejecuta los kernel
        //////////////////////////////////////////////////
        
        cudaWrapper->newSolution();
        cudaWrapper->find_minimum();

        //cudaWrapper->sortSolutions();
        //UpdateProb(saParams.count);
        id_select= 0;//selecSolution();
    
        ///////////////////////////////////////////////////
        ///  Actualiza la nueva solución en la GPU
        //////////////////////////////////////////////////
        cudaWrapper->newSolutionUpdate(costCurrentSolution, id_select);
        
        #ifdef ENABLE_GPU_RECORD_TIME
        auto start_post = std::chrono::high_resolution_clock::now();
        #endif
        ///////////////////////////////////////////////////
        ///  Verifica Error
        //////////////////////////////////////////////////
        if(costCurrentSolution<0.00 || isnan(costCurrentSolution)){
            cout << "error" << endl;
            cout << saParams.count << endl;
            std::cout << saParams.shuffle_colegios[cuParams.selectThread] << "\n";
            std::cout << saParams.shuffle_student[cuParams.selectBlock] << "\n";
            std::cout << "distancia: " << meanDist(currentSolution,distMat) << "\n";
            std::cout << "Segregación: " << S(currentSolution,alumnosSep, totalVuln) << "\n";
            std::cout << "CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";
            std::cout << "Penalty: " << penaltyParents(currentSolution,h_penalty_matrix) << "\n";
            std::cout << costCurrentSolution;
            exit(1);
        }
        
#if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_MOVE_SOLUTION
        auto move = cudaWrapper->getMovementDeviceToHost(id_select);
        recordManager->vector_historyCostSolution.emplace_back(costCurrentSolution);
        recordManager->vector_historyTemp.emplace_back(saParams.temp);
        recordManager->vector_historystu.emplace_back(std::get<0>(move));
        recordManager->vector_historycol.emplace_back(std::get<1>(move));
    #endif

    #ifdef ENABLE_OPEN_RECORD_REGISTER
        if (saParams.temp < recordManager->vector_temp_percentage[recordManager->threshold_count]){
            recordManager->vector_it_percentage[recordManager->threshold_count] = saParams.count++;
            recordManager->threshold_count+=1;
        }
    #endif
        
#endif
        
        if(costCurrentSolution < costBestSolution){
            cudaWrapper->AcceptanceBestSolution();
            costBestSolution = costCurrentSolution;
            costPreviousSolution = costCurrentSolution;
            saParams.c_accepta++;
            saParams.count_rechaso = 0;
            //cout << costCurrentSolution << " | " << saParams.count <<  endl;

#if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_REGISTER
        #ifdef ENABLE_OPEN_ITERATION_REGISTER      
            cudaWrapper->copySolutionToHost(bestSolution, previousSolution);
            recordManager->vector_costCurrentSolution.emplace_back(costBestSolution);
            recordManager->vector_meanDist.emplace_back(meanDist(bestSolution, distMat));
            recordManager->vector_segregation.emplace_back(S(bestSolution, alumnosSep, totalVuln));
            recordManager->vector_costoCupo.emplace_back(costCupo(bestSolution, cupoArray));
            recordManager->vector_penalty.emplace_back(penaltyParents(bestSolution,h_penalty_matrix));
            recordManager->vector_temp.emplace_back(saParams.temp);
            recordManager->vector_count.emplace_back(saParams.count);
        #endif
    #endif
    #ifdef ENABLE_OPEN_RECORD_MOVE_SOLUTION
        recordManager->vector_historyAcceptSolution.emplace_back(true);
    #endif
#endif

        }
        else {
            if(acceptanceCriterion->apply(costPreviousSolution,costCurrentSolution,dist_accepta ) == 1) {

                cudaWrapper->AcceptanceSolution();
                costPreviousSolution = costCurrentSolution;
#if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_MOVE_SOLUTION
                recordManager->vector_historyAcceptSolution.emplace_back(true);
    #endif
#endif
                saParams.count_rechaso = 0;
                saParams.c_accepta++;
            }
            else {
                saParams.count_rechaso++;

#if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_MOVE_SOLUTION

                recordManager->vector_historyAcceptSolution.emplace_back(false);
    #endif
#endif
                
            }
        }

        if(lengthTemperature->apply()){
            coolingScheme->apply();
        }
        reheatingMethod->apply();
        cudaWrapper->synchronizeBucle();
        saParams.count_trials++;
        saParams.count++;

        #ifdef ENABLE_GPU_RECORD_TIME
        auto end_post = std::chrono::high_resolution_clock::now();
        double time_taken_post = std::chrono::duration_cast<std::chrono::nanoseconds>(end_post - start_post).count();
        tiempoGPU << std::setprecision(10)<< time_taken_post << "\n";
        #endif

    }
    ///////////////////////////////////////////////////
    /// Obtiene el tiempo de ejecución
    ///////////////////////////////////////////////////
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cudaWrapper->copySolutionToHost(bestSolution, previousSolution);

    cout << "--------------- Resultado Final ----------------" << "\n";
    cout << "Numero de Ciclos: " << saParams.count << "\n";
    cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
    cout << "Costo de la mejor solución: " << costBestSolution << "\n";
    cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
    cout << "Tiempo de ejecución de SA: " << time_taken << "\n";
    cout << "distancia: " << meanDist(bestSolution, distMat)/saParams.max_dist << "\n"; //lo normalice
    cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << "CostoCupo: " << costCupo(bestSolution, cupoArray) << "\n";
    cout << "Penalty final: " << penaltyParents(bestSolution,h_penalty_matrix)/saParams.n_students << "\n";

    int* summaryPrefs = summaryPreferences(bestSolution, dataSet->students);
    int* data_costocupo = summaryCostoCupo(bestSolution, dataSet->colegios);
    cout << "--------------- Finalizo con exito ----------------" << "\n";

    int unassigned =summaryPrefs[saParams.max_choices];
    //int unassigned = balanceCostoCupo(bestSolution,dataSet->students, dataSet->colegios);
    
    //llamar a la funcion que lo calcula por el algoritmo original del SAE
    //std::vector<int> solution;
    //asignacionSAE(dataSet->students, dataSet->colegios, solution); //743 sin asignar en alguna pref
    //fin llamada

    

#if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_INFO
    recordManager->openRecordInfo();
    recordManager->SaveInfoFinish(costPreviousSolution,
        costBestSolution,
        costCurrentSolution,
        time_taken,
        meanDist(bestSolution, distMat),
        S(bestSolution, alumnosSep, totalVuln),
        costCupo(bestSolution, cupoArray),
        penaltyParents(bestSolution,h_penalty_matrix));
    recordManager->closeRecordInfo();
    #endif
    #ifdef ENABLE_OPEN_RECORD_GRAPHICS
    recordManager->openRecordGraphics();
    recordManager->SaveGraphicsFinish();
    recordManager->closeRecordGraphics();
    #endif
    #ifdef ENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION
    recordManager->openRecordGraphicsBestSolution();
    recordManager->SaveGraphicsBestSolution(bestSolution);
    recordManager->closeRecordGraphicsBestSolution();
    #endif
    #ifdef ENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION_RBD
    recordManager->openRecordGraphicsBestSolutionRBD();
    recordManager->SaveGraphicsUpdateSolutionRBD(bestSolution, dataSet->ptr_colegios);
    recordManager->closeRecordGraphicsBestSolutionRBD();
    #endif
    #ifdef ENABLE_OPEN_RECORD_REGISTER
    recordManager->openRecordRegister();
    recordManager->SaveInfoRegister(
        time_taken,
        costBestSolution,
        meanDist(bestSolution, distMat),
        S(bestSolution, alumnosSep, totalVuln),
        costCupo(bestSolution, cupoArray),
        penaltyParents(bestSolution,h_penalty_matrix),
        csParams.coolingRate,
        rmParams.k_reheating_init,
        rmParams.e_const,
        rmParams.k_reheating,
        ltParams.len1_init,
        ltParams.len2_init,
        ltParams.len3_init,
        ltParams.len4_init,
        ltParams.len1,
        ltParams.len2,
        ltParams.len3,
        ltParams.len4,
        acParams.Th,
        cuParams.n_block,
        cuParams.n_thread,
        bestSolution,
        summaryPrefs,
        data_costocupo[3],
        data_costocupo[0],
        data_costocupo[2]

    );
    recordManager->closeRecordRegister();
    #endif
    #ifdef ENABLE_OPEN_RECORD_MOVE_SOLUTION
    recordManager->openRecordMoveSolution();
    recordManager->AllMovementFinish();
    recordManager->closeRecordMoveSolution();
    #endif
#endif
    delete cudaWrapper;
    // cout << "finalizo con :" << costBestSolution << endl;
    return (costBestSolution);
}

template <typename T>
void SimulatedAnnealing::inicializationValues(T* wrapper){
    int x = 0, z = 0;
    totalVuln = dataSet->totalVuln;
    saParams.n_colegios = dataSet->n_colegios;
    saParams.n_students = dataSet->n_students;
    //cout << fixed << setprecision(70) << endl;
    //srand(time(NULL));


    ///////////////////////////////////////////////////
    /// Inicializa Variables y arreglos
    ///////////////////////////////////////////////////

    aluxcol= (int *)malloc(sizeof(int)*saParams.n_colegios);
    aluVulxCol = (int *)malloc(sizeof(int)*saParams.n_colegios);
    previousAluxCol = (int *)malloc(sizeof(int)*saParams.n_colegios);
    previousAluVulxCol = (int *)malloc(sizeof(int)*saParams.n_colegios);
    bestAluxCol = (int *)malloc(sizeof(int)*saParams.n_colegios);
    bestAluVulxCol = (int *)malloc(sizeof(int)*saParams.n_colegios);
    alpha = saParams.alpha;
    ptr_alpha = &saParams.alpha[0];
    
    saParams.count = 0;


    distMat=(double **)malloc(sizeof(double)*saParams.n_students);
    for(x=0; x < saParams.n_students; x++) {
        distMat[ x ]=(double *)malloc(sizeof(double)*saParams.n_colegios);
    }

    wrapper->mallocHost(
        previousSolution,
        bestSolution,
        currentSolution,
        cupoArray,
        alumnosSep,
        matrestest,
        currentVars,
        previousVars,
        bestVars);


    ///////////////////////////////////////////////////
    /// Se asignan las escuelas un arreglo que y estudiantes a la escuela
    /// las escuelas tendran como identificación el indice
    /// y currentSolution tiene como indice al estudiante y el valor del indice a la escuela que asignada
    ///////////////////////////////////////////////////
    initializeArray(aluxcol, 
                previousAluxCol, 
                bestAluxCol, 
                aluVulxCol, 
                previousAluVulxCol, 
                bestAluVulxCol, 
                alumnosSep,
                dataSet->students,
                dataSet->colegios);
    assignSchoolToArray(previousSolution, bestSolution, currentSolution, dataSet->ptr_colegios, dataSet->ptr_students, cupoArray);
    calcDist(dataSet->ptr_colegios, dataSet->ptr_students, distMat);
    saParams.max_dist = getMaxDistance(distMat);
    normalizedAlpha(alpha);

    ///////////////////////////////////////////////////
    /// Registro de datos
    ///////////////////////////////////////////////////

    //calculo de la matriz de penalidades
    std::vector<int> num_ele_vec;
    num_ele_vec.reserve(dataSet->students.size());

    for (const auto& alu : dataSet->students) {
        num_ele_vec.push_back(alu.num_ele);
    }

    std::vector<int> choices_matrix;
    choices_matrix.reserve(saParams.n_students * saParams.max_choices);

    for (const auto& alu : dataSet->students) {
        //solo tomo los primeros saParams.max_choices elementos
        for (int j = 0; j < saParams.max_choices; ++j) {
            choices_matrix.push_back(static_cast<int>(alu.choices[j]));
        }
    }

    //std::vector<float> penalty_matrix(saParams.n_students * saParams.n_colegios);
    h_penalty_matrix = new float[saParams.n_students * saParams.n_colegios];

    wrapper->compute_penalty_matrix(choices_matrix.data(),
                            num_ele_vec.data(),
                            h_penalty_matrix);
    
/*   
    //borrable???
    //creo que aqui se calculo de forma normal, pero mas adelante se sobreescribe el valor con uno normalizado
    //entonces calCosto no es la version CPU
    costBestSolution = calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;

    //fin borrable ------
*/
    saParams.count++;

    ///////////////////////////////////////////////////
    /// Genera distribuciones para seleccionar un estudiante y una escuela al azar
    ///////////////////////////////////////////////////
    
    dist = uniform_int_distribution<int>(0, saParams.n_students-1);
    dist2 = uniform_int_distribution<int>(0, saParams.n_colegios-1);

    saParams.count_rechaso=0;
    saParams.c_accepta = 0;
    saParams.c_cooling_temperature = 0;
    rmParams.k_reheating_init = rmParams.k_reheating;
    saParams.temp_init = saParams.temp;
    saParams.count_trials = 0;
    ltParams.len1_init = ltParams.len1;
    ltParams.len2_init = ltParams.len2;
    ltParams.len3_init = ltParams.len3;
    ltParams.len4_init = ltParams.len4;

    ////////////////////////////////////////////////////////////////////////
    // VARIABLES DE PRUEBA
    ////////////////////////////////////////////////////////////////////////

    
    // double costCurrentSolutionV2 = costCurrentSolution;
    
    currentVars[0] = sumDist(currentSolution,distMat);
    currentVars[1] = sumS(currentSolution, alumnosSep, totalVuln);
    currentVars[2] = sumCostCupo(currentSolution,cupoArray);
    currentVars[3] = penaltyParents(currentSolution,h_penalty_matrix);
    previousVars[0] = currentVars[0];
    previousVars[1] = currentVars[1];
    previousVars[2] = currentVars[2];
    previousVars[3] = currentVars[3];
    
    //porque llama a sumdist y luego le hace la misma division
    //o sea hay una funcion adicional que lo unico que hace es no hacer una division
    //y esa misma division se hace despues
    //lo mismo para las demas funciones
    //porque sumS divide en 2, que significa ese 2

    double var1,var2,var3,var4;
    var1 = (currentVars[0]/saParams.n_students);
    var1= (var1/saParams.max_dist);
    //cout << var1 << "\n";
    var2 = (currentVars[1]/2.0);
    //cout << var2 << "\n";
    var3 = (currentVars[2] /saParams.n_colegios);
    //cout << var3 << "\n";
    var4 = currentVars[3] /saParams.n_students;
    costBestSolution = (double)((ptr_alpha[0] * var1) + (ptr_alpha[1] * var2) + (ptr_alpha[2] * var3) + (ptr_alpha[3] * var4));
    //cout << "F2: " << costBestSolution << "\n";

    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;
    auto start_compare = std::chrono::high_resolution_clock::now();
    auto end_compare = std::chrono::high_resolution_clock::now();
    double time_taken_v1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_compare - start_compare).count();
    //double time_taken_v2;
    //double vector_time1 =0;
    // double vector_time2=0;
    saParams.max_changes_students = min(cuParams.n_thread*cuParams.n_block, saParams.n_students);
    saParams.max_changes_school = saParams.n_colegios;//min(cuParams.n_block, saParams.n_colegios);

    double *array_costCurrentSolution = (double *) malloc(sizeof(double) * cuParams.n_block * cuParams.n_thread);
    for (x = 0; x < cuParams.n_block; x++){
        for (z = 0; z < cuParams.n_thread; z++){
            array_costCurrentSolution[cuParams.n_thread * x + z] = 0.0;
        }
    }
    ///////////////////////////////////////////////////
    /// Inicializa las distribuciónes
    ///////////////////////////////////////////////////
    for (int x = 0; x < saParams.n_students; x++) {
        for (int z = 0; z < saParams.n_colegios; z++) {
            matrestest[saParams.n_colegios * x + z] = distMat[x][z];
        }
    }

    dist = std::uniform_int_distribution<int>(0, saParams.n_students-1);
    dist2 = std::uniform_int_distribution<int>(0, saParams.n_colegios-1);
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////


///////////////////////////////////////////////////
/// Calcula el costo
///////////////////////////////////////////////////
double SimulatedAnnealing::calCosto(int *currentSolution, double **distMat, const double *ptr_alpha, int *alumnosSep, int totalVuln, int *cupoArray){
    double var1 = meanDist(currentSolution,distMat)/saParams.max_dist;
    //cout << "distancia: " << var1 << "\n";
    double var2 = S(currentSolution, alumnosSep, totalVuln);
    //cout << "Segregación: " << var2 << "\n";
    double var3 = costCupo(currentSolution,cupoArray);
    //cout << "CostoCupo: " << var3 << "\n";
    double var4 = penaltyParents(currentSolution, h_penalty_matrix)/(saParams.n_students);
    //cout << "Penalty: " << var4 << "\n";

    return (double)((ptr_alpha[0] * var1) + (ptr_alpha[1] * var2) + (ptr_alpha[2] * var3) + (ptr_alpha[3] * var4));
}

///////////////////////////////////////////////////
/// Distancia promedio que recorreran los estudiantes
///////////////////////////////////////////////////
double SimulatedAnnealing::meanDist(const int *currentSolution, double  **distMat){
    double sumDist=0.0;
    for(int i=0;i<saParams.n_students;i++){
        sumDist+=distMat[i][currentSolution[i]]; // distMat[estudiante][escuela]
    }
    //cout << "meanDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << saParams.n_students << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist/saParams.n_students;
}

double SimulatedAnnealing::sumDist(const int *currentSolution, double  **distMat){
    double sumDist=0.0;
    for(int i=0;i<saParams.n_students;i++){
        sumDist+=distMat[i][currentSolution[i]]; // distMat[estudiante][escuela]
    }
    //cout << "sumDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << saParams.n_students << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist;
}


///////////////////////////////////////////////////
/// Calcula segregación por duncan
///////////////////////////////////////////////////

double SimulatedAnnealing::S(const int *currentSolution,const int *alumnosSep, int totalVuln){
    double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<saParams.n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < saParams.n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(saParams.n_students-totalVuln)));
        }
    }
    return totalSesc/2.0;
}

double SimulatedAnnealing::sumS(const int *currentSolution,const int *alumnosSep, int totalVuln){
    double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<saParams.n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < saParams.n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(saParams.n_students-totalVuln)));
        }
    }
    return totalSesc;
}


///////////////////////////////////////////////////
/// Calcula el costo de tener los estudiantes en las escuelas
///////////////////////////////////////////////////

double SimulatedAnnealing::costCupo(int *currentSolution,int *cupoArray){
    double totalcostCupo = 0.0;
    int totalAluCol = 0;
    for(int j=0;j<saParams.n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<saParams.n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        double p_costCupo = double(totalAluCol)/cupoArray[j];
        totalcostCupo += calcCostoCupo(p_costCupo);
        //totalcostCupo+=(double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2));
        //revisar cuantos colegios ha cerrado
    }
    return totalcostCupo/saParams.n_colegios;
}

double SimulatedAnnealing::calcCostoCupo(double p_costCupo) {
    double r = 0.6;
    double s = 6.0;
    int l_izq = p_costCupo <= 0.5 && p_costCupo >= 0.0; 
    int l_der = p_costCupo > 0.5 && p_costCupo <= 1.0; 

    double penalty_sobrecupo = calcCostoCupo_sobrecupo(p_costCupo);

    double costCupoEscuela =  l_izq*(pow(2.0, r)*pow(p_costCupo, r)) + l_der*(pow(2.0, s) * pow(1.0 - p_costCupo, s))+penalty_sobrecupo;
    //cout <<l_izq<< l_der<<"p_costCupo: "<<p_costCupo<<" | result: "<<costCupoEscuela<<" l_der"<< l_der*pow(2.0, s) * pow(1.0 - p_costCupo, s)<<"\n";
    return costCupoEscuela;
}

double SimulatedAnnealing::calcCostoCupo_sobrecupo(double p_costCupo) {
    const double k = 10.0;  // controla la pendiente del crecimiento exponencial
    const double S = 1.0;   // factor de escala (ajústalo al resto de tus costos)

    double over = fmax(0.0, p_costCupo - 1.0);
    double penalty = S * (exp(k * over) - 1.0);

    return penalty;
}



double SimulatedAnnealing::sumCostCupo(int* currentSolution,int *cupoArray){
    double totalcostCupo = 0.0;
    int totalAluCol = 0;
    for(int j=0;j<saParams.n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<saParams.n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        double p_costCupo = double(totalAluCol)/cupoArray[j];
        totalcostCupo += calcCostoCupo(p_costCupo);
        //totalcostCupo+= (double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2));
    }
    return totalcostCupo;
}


///////////////////////////////////////////////////
/// Genera una nueva solución en donde asigna a un estudiante a una escuela
/// aleatoriamente
///////////////////////////////////////////////////

/*
void SimulatedAnnealing::newSolution(int *currentSolution,const int *previousSolution){
    //random_device rd;
    //mt19937 mt(rd());
    uniform_int_distribution<int> dist(0, saParams.n_students);
    random_device rd2;
    mt19937 mt2(rd2());
    uniform_int_distribution<int> dist2(0, saParams.n_colegios);
    int selectStudent=dist(mt);
    int selectSchool = dist2(mt2);
    for(int x=0; x<saParams.n_students; x++){
        if(x == selectStudent) {
            currentSolution[x] = selectSchool;
        }
        else {
            currentSolution[x] = previousSolution[x];
        }
    }

}

*/


///////////////////////////////////////////////////
/// Asigna a las soluciones la escuela actual Solo se utiliza al inicio
///////////////////////////////////////////////////
void SimulatedAnnealing::assignSchoolToArray(int *previousSolution, int *bestSolution, int *currentSolution, Info_colegio *ptr_colegios, Info_alu *ptr_students, int *cupoArray){
    Info_alu *ptr_aux = ptr_students;
    for(int x=0;x < saParams.n_colegios;x++){
        for(int y=0; y < saParams.n_students; y++){
            if(ptr_colegios->rbd == ptr_students->rbd){
                previousSolution[y] = x;
                bestSolution[y] = x;
                currentSolution[y] = x;
            }
            ptr_students++;

        }
        /*
         * cupoArray sera un arreglo que por indice es la escuela y su valor sera el cupo que posee esa escuela
         * se asume que las escuelas pueden tener sobre cupo.
         */

        //cupoArray[x] = ptr_colegios->num_alu+ ((int)((ptr_colegios->num_alu*10)/100));
        cupoArray[x] = ptr_colegios->cupos;
        ptr_students = ptr_aux;
        ptr_colegios++;
    }
}
///////////////////////////////////////////////////
/// Crea una matriz de distancia donde x es el estudiante, y es la escuela
///////////////////////////////////////////////////
void SimulatedAnnealing::calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, double **distMat){
    Info_colegio *ptr_aux = ptr_colegios;
    for(int x=0;x < saParams.n_students ;x++){
        for(int y=0; y < saParams.n_colegios; y++){
            distMat[x][y] = sqrt( pow((ptr_students->latitude - ptr_colegios->latitude),2)+pow((ptr_students->longitude - ptr_colegios->longitude),2))/1000;
            ptr_colegios++;

        }
        ptr_colegios = ptr_aux;
        ptr_students++;
    }
}



void SimulatedAnnealing::shuffle(int *values, const int max_change, uniform_int_distribution<int> distri) {
    int randvalue1,randvalue2,tem_value;
    for (int i = 0; i<max_change; i++) {
        randvalue1 = distri(mt);
        randvalue2 = i;
        tem_value = values[randvalue1];
        values[randvalue1] = values[randvalue2];
        values[randvalue2] = tem_value;
    }
}

////////////////////////////////////////////////
////// Obtiene la maxima distancia que un estudiante podria llegar a recorrer
///////////////////////////////////////////////////
double SimulatedAnnealing::getMaxDistance(double **distMat){
    double max = 0;
    for(int i=0;i<saParams.n_students;i++){
        for(int x=0;x<saParams.n_colegios;x++){
            if(distMat[i][x]>max){
                max = distMat[i][x];
            }
        }
    }
    return max;
}

///////////////////////////////////////////////////
/// Calcula el valor de los alpha
///////////////////////////////////////////////////
void SimulatedAnnealing::normalizedAlpha(double *alpha)
{
    double sumaAlpha = 0.0;
    for(int x=0; x<4; x++){
        sumaAlpha +=alpha[x];
    }
    for(int x=0; x<4; x++){
        alpha[x]= alpha[x]/(double)sumaAlpha;
    }
}


///////////////////////////////////////////////////
/// Asigna Información de las escuelas a best, previus y current soluciones
///////////////////////////////////////////////////
void SimulatedAnnealing::initializeArray(int *aluxcol, int *previousAluxCol, int *bestAluxCol, int *aluVulxCol, int *previousAluVulxCol, int *bestAluVulxCol, int *alumnosSep, vector<Info_alu> &students,vector<Info_colegio> &colegios)
{
    for(int x = 0; x < saParams.n_colegios; x++){
        aluxcol[x] = colegios[x].num_alu;
        previousAluxCol[x] = colegios[x].num_alu;
        bestAluxCol[x] = colegios[x].num_alu;
        aluVulxCol[x] = colegios[x].prioritario;
        previousAluVulxCol[x] = colegios[x].prioritario;
        bestAluVulxCol[x] = colegios[x].prioritario;

    }
    ///////////////////////////////////////////////////
    /// Se crear un arreglo donde el el valor es la posición del estudiante sep
    ///////////////////////////////////////////////////
    for(int x=0; x < saParams.n_students; x++) {
        alumnosSep[x] = students[x].sep;
    }
}


double SimulatedAnnealing::round_n(double x)
{
    double digits = pow(10.0, DECIMAL);
    return trunc(x * digits) / digits;
}


int SimulatedAnnealing::acceptanceCriterionApply() {
    return acceptanceCriterion->apply(costPreviousSolution,costCurrentSolution,dist_accepta);
}


int SimulatedAnnealing::selecSolution(){
    
    size_t size = 16;//cuParams.n_block*cuParams.n_thread;
    for (size_t x = 0; x < size; x++){
        double select = dist_accepta(mt);
        if (select<probSelection[x]) return x* ((cuParams.n_block*cuParams.n_thread)/16);
    }
    return size-1;
   
}

void SimulatedAnnealing::UpdateProb(int it){
        saParams.p = saParams.pMax - (saParams.pMax - saParams.pInit) * exp(-saParams.k*it);

        int size = 16;//cuParams.n_block*cuParams.n_thread;
        probSelection[0] = saParams.p;
        double sum = saParams.p; 

        for (int x = 1; x < size; x++) {
            probSelection[x] = probSelection[x - 1] * (1 - saParams.p);
            sum += probSelection[x];
        }

        // Normalización y cálculo acumulativo en un paso
        double acumulado = 0.0;
        for (int x = 0; x < size; x++) {
            probSelection[x] /= sum; // Normaliza
            acumulado += probSelection[x]; // Acumula
            probSelection[x] = acumulado;
        }

        // Asegurar que el último elemento sea 1
        probSelection[size - 1] = 1;
}

/*
funcion original del oscar
std::size_t SimulatedAnnealing::penaltyParents(int *currentSolution) {
    std::array<std::size_t, 6> weights{0, 100, 200, 300, 400, 500000};

    std::size_t penalty = 0;
    bool find = false;

    for (std::size_t i = 0; i < saParams.n_students; i++) {
        for (std::size_t j = 0; j < 5; j++) {
            if (currentSolution[i] == choices_parents[i * 5 + j]) {
                penalty += weights[j];
                find = true;
                break;
            }
        }

        if (!find)
            penalty += weights[5];

        find = false;
    }

    return penalty;
}
*/

double SimulatedAnnealing::penaltyParents(int *currentSolution, float* h_penalty_matrix) {

    double penalty = 0;

    for (std::size_t i = 0; i < saParams.n_students; i++) {
        penalty += h_penalty_matrix[i*saParams.n_colegios +currentSolution[i]];
    }

    return penalty;
}


void SimulatedAnnealing::ValidateGPU(){
    CUDAWrapper* cudaWrapper = new CUDAWrapper(cuParams, saParams, mt);
    inicializationValues(cudaWrapper);
    cudaWrapper->memInit(previousSolution,
        bestSolution,
        currentSolution,
        cupoArray,
        alumnosSep,
        totalVuln,
        aluxcol,
        aluVulxCol,
        matrestest,
        alpha,
        currentVars);


    std::vector<int> solution;
    int n_ciclos = 100000;
    std::vector<DataResult> sol;
    int id_select= 0;
    DataResult cpu;
    double costTempSol;
    double dif;
    double max_error = -1;
    int ciclo_max_error = -1;
    double previousSolutionUnitTest;
    double epsilon = 1e-15;

    for (int ciclo = 0; ciclo < n_ciclos; ++ciclo) {
        cudaWrapper->memCopyPrevToCurrent();
            
        shuffle(saParams.shuffle_student, saParams.max_changes_students, dist);
        shuffle(saParams.shuffle_colegios, saParams.max_changes_school, dist2);


        cudaWrapper->uploadCurrentMemorySolution();

        //primer kernel (mejor solucion por bloque)
        cudaWrapper->newSolution();
        cudaWrapper->find_minimum(); //encontrar la mejor solucion del bloque
        cudaWrapper->getSolution(sol); //mi wrapper que obtiene los datos


        auto r = sol[id_select];
        //mostrar accion realizada
        cout << "iteracion: "<< ciclo << " accion: -> stu= " << r.stu << " -> col= " << r.col << " costo=" << r.costSolution << "\n";
        
        cudaWrapper->newSolutionUpdate(costCurrentSolution, id_select);
        cudaWrapper->previousSolution(id_select);
        cudaWrapper->getPreviousSolutionUnitTest(previousSolutionUnitTest);

        if((previousSolutionUnitTest-costPreviousSolution) > epsilon){
            cout << "En el ciclo " << ciclo << " el calculo del valor anterior no es igual\n";
        }


        //copypaste del original
        if(costCurrentSolution < costBestSolution){
                cudaWrapper->AcceptanceBestSolution();
                cudaWrapper->copySolutionToHost(bestSolution, previousSolution);
                costTempSol = calCosto(bestSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
                costBestSolution = costCurrentSolution;
                costPreviousSolution = costCurrentSolution;
                saParams.c_accepta++;
                saParams.count_rechaso = 0;

                cpu = cpu_one_tid_newSolution(r.stu, r.col); //funcion CPU que imita al kernel

                dif = costTempSol - costBestSolution;
                if (dif>0){
                    if (dif > epsilon){
                        cout << "iteracion N " << ciclo+1 << " : Resultado no es igual a CPU, error de "<< dif<< "\n";
                        cout << "accion realizada: -> stu= " << r.stu << " -> col= " << r.col << " costo=" << r.costSolution << "\n";  
                    }
                    
                    if (max_error < dif){
                        max_error = dif;
                        ciclo_max_error = ciclo;
                    }
                }

            }
            else {
                //cout << "iteracion N " << ciclo << " no es mejor solucion que la anterior\n";
                    if(acceptanceCriterion->apply(costPreviousSolution,costCurrentSolution,dist_accepta ) == 1) {
                    //cudaWrapper->AcceptanceSolution();
                    costPreviousSolution = costCurrentSolution;
                    saParams.count_rechaso = 0;
                    saParams.c_accepta++;
                }
                else {
                    saParams.count_rechaso++;
                    
                }
            }
        }
    
    //cout << "Despues: " << bestSolution[r.stu] << "\n";

    double costBestSol = calCosto(bestSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);

    cout << setprecision(16) << "CPU thread: " << cpu.costSolution << "\n";
    cout << setprecision(16) << "Total CPU:  " << costBestSol << "\n";
    cout << setprecision(16) << "Total GPU:  " << costBestSolution << "\n";
    cout << setprecision(16) << "Error maximo: " << max_error << " ciclo: " << ciclo_max_error << "\n";

    delete cudaWrapper;
}

//imitacion del kernel de un movimiento a la vez, pero en CPU para pruebas unitarias 
DataResult SimulatedAnnealing::cpu_one_tid_newSolution(int aluchange, int newSchool) {
    int n_students   = saParams.n_students;
    int n_colegios   = saParams.n_colegios;
    int currentSchool= currentSolution[aluchange];

    DataResult out;
    out.stu = aluchange;
    out.col = newSchool;

    //copia de los valores
    double sumDist     = currentVars[0];
    double totalSesc   = currentVars[1];
    double totalCupo   = currentVars[2];
    double penalty     = currentVars[3];

    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////

    // Distancia
    sumDist -= distMat[aluchange][currentSchool];

    // Seg de la escuela actual
    int totA = aluxcol[currentSchool];
    int vulA = aluVulxCol[currentSchool];
    int novA = totA - vulA;
    totalSesc -= std::fabs((vulA/(double)totalVuln) - (novA/(double)(n_students - totalVuln)));

    // Costocupo escuela actual 
    double p_costCupo = (double)totA/cupoArray[currentSchool];
    totalCupo -= calcCostoCupo(p_costCupo);

    //std::cout << std::setprecision(17) << "CPU Vars: " << totA << " " << cupoArray[currentSchool] << "\n";
    // Seg de la escuela nueva
    int totB = aluxcol[newSchool];
    int vulB = aluVulxCol[newSchool];
    int novB = totB - vulB;
    totalSesc -= std::fabs( (vulB/(double)totalVuln) - (novB/(double)(n_students - totalVuln)) );
    
    //costocupo escuela nueva
    p_costCupo = (double)totB/cupoArray[newSchool];
    totalCupo -= calcCostoCupo(p_costCupo);

    ////////////////////////////////////////////////////////////////
    /////// Realiza Movimiento
    ////////////////////////////////////////////////////////////////

    //elimina el estudiante de la escuela actual
    aluxcol[currentSchool]-=1;
    aluVulxCol[currentSchool]-=alumnosSep[aluchange];
    //Asigna al estudiante a la nueva escuela
    currentSolution[aluchange] = newSchool;
    aluxcol[newSchool]+=1;
    aluVulxCol[newSchool]+=alumnosSep[aluchange];

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////

    //distancia de la nueva escuela
    sumDist += distMat[aluchange][newSchool];

    //seg de la escuela actual
    int totA2 = aluxcol[currentSchool];
    int vulA2 = aluVulxCol[currentSchool];
    int novA2 = totA2 - vulA2;
    totalSesc += std::fabs( (vulA2/(double)totalVuln) - (novA2/(double)(n_students - totalVuln)) );

    //costocupo escuela actual
    p_costCupo = double(totA2)/cupoArray[currentSchool];
    totalCupo += calcCostoCupo(p_costCupo);

    //seg de la escuela antigua
    int totB2 = aluxcol[newSchool];
    int vulB2 = aluVulxCol[newSchool];
    int novB2 = totB2 - vulB2;
    totalSesc += std::fabs( (vulB2/(double)totalVuln) - (novB2/(double)(n_students - totalVuln)) );

    //costcupo escuela antigua
    p_costCupo = double(totB2)/cupoArray[newSchool];
    totalCupo += calcCostoCupo(p_costCupo);


    //Calculo de penalty
    //penalty de la escuela actual
    penalty   -= h_penalty_matrix[aluchange*saParams.n_colegios +currentSchool];

    //penalty de la escuela nueva
    penalty   += h_penalty_matrix[aluchange*saParams.n_colegios +newSchool];

    //std::cout << std::setprecision(17) << "CPU cost 2 B: " << totalCupo << "\n";
    currentVars[0] = sumDist;
    currentVars[1] = totalSesc;
    currentVars[2] = totalCupo;
    currentVars[3] = penalty;

    //std::cout << std::setprecision(17) << "CPU 2: " << currentVars[0] << " " << currentVars[1] << " "<< currentVars[2] << " "<< currentVars[3] << "\n";


    // Combinar con alphas (mismas normalizaciones que el kernel)
    const double var1 = (sumDist / (double)n_students) / saParams.max_dist;
    const double var2 = (totalSesc * 0.5);
    const double var3 = (totalCupo / (double)n_colegios);
    const double var4 = (penalty / n_students);

    out.costSolution = alpha[0]*var1 + alpha[1]*var2 + alpha[2]*var3 + alpha[3]*var4;
    return out;
}







//para contar cuantas asignaciones se hicieron en primera, segunda... N preferencia
int* SimulatedAnnealing::summaryPreferences(const int* currentSolution,
                                const std::vector<Info_alu>& alumnos)
{

    //las primeras max_choices elementos son las preferencias y la ultima son los que estan fuera
    int* asignacion_por_pref = new int[saParams.max_choices + 1]();
    
    for (int i = 0; i < (int)alumnos.size(); ++i) {
        int col_asignado = currentSolution[i]+1; //en el dataset van de 1 a 63, pero en el 
        const auto& A = alumnos[i];
        bool encontrado = false;

        //iteramos por num_ele a cada alumno
        for (int pref = 0; pref < A.num_ele; ++pref) {
            if (col_asignado == A.choices[pref]) {
                asignacion_por_pref[pref]++;
                encontrado = true;
                break;
            }
        }

        if (!encontrado) {
            asignacion_por_pref[saParams.max_choices]++;
        }
    }

    std::cout << "Asignados en preferencia numero ->";
    for (int p = 0; p < saParams.max_choices; ++p){
    std::cout << (p+1) << ":" << asignacion_por_pref[p] << "  ";
}
    std::cout << "99:" << asignacion_por_pref[saParams.max_choices] << "  ";
    std::cout << "\n";

    return asignacion_por_pref; 
}


int* SimulatedAnnealing::summaryCostoCupo(const int* currentSolution,
                                          const std::vector<Info_colegio>& colegios)
{
    int n_colegios = saParams.n_colegios;

    //recuento de alumnos por colegio
    int* ocupados = new int[n_colegios]();
    int* vulnerables = new int[n_colegios]();
    
    // contar ocupación actual
    for (int i = 0; i < saParams.n_students; ++i) {
        ocupados[currentSolution[i]]++;
        vulnerables[currentSolution[i]] += alumnosSep[i];
    }

    //resumen
    //[0] = colegios con cupos
    //[1] = cupos restantes totales
    //[2] = colegios en sobrecupo
    int* resumen = new int[4]();

    std::ofstream fout("summaryCupos.txt");
    if (!fout.is_open()) {
        std::cerr << "Error abriendo summaryCupos.txt\n";
    }
    int capacidad_total = 0;
    for (int j = 0; j < n_colegios; ++j) {
        //int capacidad = (int)std::floor(colegios[j].num_alu * 1.1); 
        int capacidad = colegios[j].cupos;
        capacidad_total += capacidad;
        int ocup = ocupados[j];
        int vul = vulnerables[j];
        int vacantes = capacidad - ocup;

        if (vacantes > 0) {
            resumen[0]++;          // colegios con cupos
            //std::cout << "col con cupo: "<< j<< " | Vacantes: " << vacantes << "\n";
            resumen[1] += vacantes; // cupos restantes totales
        }
        if (ocup > capacidad) {
            resumen[2]++;          // colegios en sobrecupo
        }
        fout << colegios[j].rbd << "\t" << capacidad << "\t" << ocup << "\t" << vul << "\n";
    }
    fout.close();
    resumen[4] = resumen[1] - (capacidad_total-saParams.n_students);
    std::cout << "\n" << "Colegios aun con vacantes: " << resumen[0]
            << " | Vacantes totales restantes: " << resumen[1] << " | colegios en sobrecupo: " << resumen[2] 
            << " | Alumnos en sobrecupo "<< resumen[4] << "\n";
    return resumen;

}


inline uint32_t lottery_simple(int student_index, int rbd, uint32_t seed)
{
    uint32_t x = seed;

    // Mezcla fuerte estilo MurmurHash3
    x ^= (uint32_t)student_index * 0x27d4eb2dU;
    x ^= (uint32_t)rbd           * 0x165667b1U;

    // Finalizador de 32 bits (muy buena avalancha)
    x ^= x >> 15;
    x *= 0x85ebca6bU;
    x ^= x >> 13;
    x *= 0xc2b2ae35U;
    x ^= x >> 16;

    return x;
}


//funcion que replica la asignacion de alumnos del SAE en base a la documentacion oficial (bases de datos, SAE)
void SimulatedAnnealing::asignacionSAE(
    const std::vector<Info_alu>& alumnos,
    const std::vector<Info_colegio>& colegios,
    std::vector<int>& solution
){
    struct Postulante {
        int i;          //numero de estudiante
        int prio;       //prioridad
        uint32_t lotto; //numero obtenido en la loteria aleatoria (lottery_simple)
    };

    
    solution.assign(saParams.n_students, -1);

    //Variables Auxiliares

    //capacidad de cada establecimiento
    std::unordered_map<int,int> capacidad;
    capacidad.reserve(colegios.size()*2);

    //vacantes SEP de cada establecimiento
    std::unordered_map<int,int> cuota_sep;
    cuota_sep.reserve(colegios.size()*2);

    //"admitidos parciales" de cada establecimiento 
    std::unordered_map<int, std::vector<Postulante>> admitidos;
    admitidos.reserve(colegios.size()*2);

    std::vector<int> next_pref(saParams.n_students, 0);
    std::vector<bool> esta_admitido(saParams.n_students, false);
    std::vector<bool> en_nuevos(saParams.n_students, false);

    for (const auto& C : colegios) {
        capacidad[C.rbd] = C.cupos;
        cuota_sep[C.rbd] = C.vacantes_prio;
    }

    //================= FASE 1 =================
    //En la doc del sae, proceso 2017 en el cual se baso, se habla de 2 fases
    //la primera es la tombola tradicional, asignados en funcion de una lista de prioridades y un numero aleatorio
    //la fase 2 entra en accion cuando hay postulantes sin cupo asignado, donde en el proceso inicial
    //se habla de que deben ser asignados a su colegio mas cercano, partiendo del estudiante con el colegio cercano mas lejano


    //en la fase 1 la logica sera de que se movera hasta que no se hallan realizado cambios
    //definiendo un cambio como un estudiante que postulo a un establecimiento y fue rechazado, postulando su siguiente opcion
    //o un estudiante que fue asignado a un establecimiento fue desplazado por un estudiante de mayor prioridad
    bool cambios = true;

    while (cambios) {
        cambios = false;

        //revision de los alumnos que estan siendo asignados "tentativamente"
        std::fill(esta_admitido.begin(), esta_admitido.end(), false);
        for (const auto& kv : admitidos)
            for (const auto& p : kv.second)
                esta_admitido[p.i] = true;

        //lista de postulantes        
        std::unordered_map<int, std::vector<Postulante>> listas;
        listas.reserve(colegios.size()*2);

        //preparacion de la postulacion de cada estudiante
        for (int i = 0; i < saParams.n_students; ++i){

            //aquellos que estan siendo asignados "tentativamente" no realizan postulacion en dicha iteracion
            if (esta_admitido[i]) continue;

            int p = next_pref[i];
            const auto& A = alumnos[i];

            //ademas de los que ya no tienen colegios a los que postular
            if (p >= A.num_ele) continue;

            //el -1, es porque en el dataset los colegios se enumeran de 1 hasta N, en cambio en C
            //se enumera del 0 hasta N-1
            size_t idx = static_cast<size_t>(A.choices[p]) - 1;
            int rbd_pref = colegios[idx].rbd;

            //la mayoria de las prioridades vienen definidas por el dataset (database oficial)
            //sin embargo la prioridad 4 se realiza una mini tombola referente a los alumnos prioritarios
            int prio = A.prioridades[p];

            listas[rbd_pref].push_back({
                i, prio, lottery_simple(i, rbd_pref, saParams.seed)
            });
        }

        //se realiza la postulacion de los establecimientos a su respectiva escuela
        for (auto& kv : listas){

            int rbd = kv.first;
            auto& nuevos_postulantes = kv.second;

            //V es la cantidad de vacantes del colegio==rbd
            //V en el dataset es la capacidad absoluta, o sea es la sumatoria de los cupos normales y SAE
            int V = capacidad[rbd];
            if (V <= 0) continue;

            //cuota SEP
            int cuota = cuota_sep[rbd];

            const auto& antiguos = admitidos[rbd];

            //la tombola en cada iteracion es una pool de estudiantes nuevos y antiguos, donde los nuevos
            //pueden desplazar a un estudiante asignado "tentativamente" (galey shappley)
            std::vector<Postulante> pool = antiguos;
            pool.insert(pool.end(), nuevos_postulantes.begin(), nuevos_postulantes.end());

            //la documentacion dice que los primeros "cuota" elementos (cantidad de cupos SEP del establecimiento)
            //seleccionados aleatoreamente entran con prioridad 4, por lo que se realiza una mini tombola
            //para seleccionar los que entran con una prioridad mas "favorable"
            std::vector<Postulante*> seps;
            for (auto& p : pool)
                if (p.prio == 4)
                    seps.push_back(&p);

            //simplemente se ordena en funcion del numero random obtenido
            std::sort(seps.begin(), seps.end(),
                [](const Postulante* a, const Postulante* b){
                    if (a->lotto != b->lotto) return a->lotto < b->lotto;
                    return a->i < b->i;
                });

            if ((int)seps.size() > cuota){
                for (size_t t = cuota; t < seps.size(); ++t)
                    seps[t]->prio = 7;
            }

            //ahora realizamos la tombola principal, donde seleccionamos los primeros V
            //postulantes y los asignamos tentativamente
            std::sort(pool.begin(), pool.end(),
                [](const Postulante& a, const Postulante& b){
                    if (a.prio != b.prio) return a.prio < b.prio;
                    if (a.lotto != b.lotto) return a.lotto < b.lotto;
                    return a.i < b.i;
                });

            std::vector<Postulante> nuevos_admitidos;
            for (int i = 0; i < std::min(V, (int)pool.size()); ++i)
                nuevos_admitidos.push_back(pool[i]);

            //guardamos registros de los cambios 
            std::fill(en_nuevos.begin(), en_nuevos.end(), false);
            for (auto& p : nuevos_admitidos)
                en_nuevos[p.i] = true;

            //aquellos desplazados y rechazados se aumenta en 1 su puntero
            for (auto& p : antiguos)
                if (!en_nuevos[p.i]){
                    next_pref[p.i]++;
                    cambios = true;
                }
            
            for (auto& p : nuevos_postulantes)
                if (!en_nuevos[p.i]){
                    next_pref[p.i]++;
                    cambios = true;
                }
            
            //actualizamos el vector solucion
            admitidos[rbd] = std::move(nuevos_admitidos);
        }
    }

    //================= CONSTRUIR SOLUCION =================
    for (const auto& kv : admitidos){
        int rbd = kv.first;
        for (const auto& p : kv.second)
            solution[p.i] = rbd;
    }

    //================= FASE 2 =================
    //como mencionamos anteriormente, en la implementacion 2017, en la cual esta basada
    //se aplica una fase de heuristica, donde hay que asignar a los estudiantes a su colegio
    //mas cercano, partiendo del estudiante mas lejano


    //buscamos las vacantes que quedan en el conjunto de colegios que aun no se han utilizado
    std::unordered_map<int,int> vac_by_rbd;
    vac_by_rbd.reserve(colegios.size()*2);

    long long total_vacantes_restantes = 0;

    for (const auto& C : colegios) {
        int usados = 0;
        auto it = admitidos.find(C.rbd);
        if (it != admitidos.end())
            usados = (int)it->second.size();

        int vac = std::max(0, C.cupos - usados);
        vac_by_rbd[C.rbd] = vac;
        total_vacantes_restantes += vac;
    }

    //hash para trabajar con rbd en lugar de indices
    //es para manejarme mejor, me enreda lo de restar -1 a los indices
    std::unordered_map<int,int> rbd2idx;
    rbd2idx.reserve(colegios.size()*2);
    std::vector<int> idx2rbd(colegios.size());

    for (int j = 0; j < (int)colegios.size(); ++j){
        rbd2idx[colegios[j].rbd] = j;
        idx2rbd[j] = colegios[j].rbd;
    }

    //obtenemos los estudiantes que no fueron asignados
    std::vector<int> unassigned;
    unassigned.reserve(saParams.n_students);

    for (int i = 0; i < saParams.n_students; ++i)
        if (solution[i] == -1)
            unassigned.push_back(i);

    //helper que ayuda a rearmar el conjunto de colegios
    //se activa cuando un colegio se llena        
    auto rebuild = [&](){
        std::vector<int> v;
        for (int j = 0; j < (int)colegios.size(); ++j)
            if (vac_by_rbd[idx2rbd[j]] > 0)
                v.push_back(j);
        return v;
    };

    std::vector<int> schools = rebuild();

    //struct para ordenamiento de los postulantes (menor distancia, colegio)
    struct Best { double d; int j; };

    //iteraremos mientras existan estudiantes inasignados y existan vacantes
    while (!unassigned.empty() && total_vacantes_restantes > 0 && !schools.empty()){

        //la posicion en vector marca el estudiante
        std::vector<Best> best(unassigned.size(), {1e18, -1});

        //calculamos distancias hacia el colegio mas cercano de cada estudiante
        for (size_t t = 0; t < unassigned.size(); ++t){
            int i = unassigned[t];
            for (int j : schools){
                int rbd = idx2rbd[j];
                if (vac_by_rbd[rbd] <= 0) continue;
                double d = distMat[i][j];
                if (d < best[t].d){
                    best[t] = {d, j};
                }
            }
        }

        //se tiene que buscar a el estudiante cuyo colegio cercano es el mas lejano de todas las menores distancias estudiante escuela
        int pick = -1;
        double worst = -1;

        for (size_t t = 0; t < unassigned.size(); ++t){
            if (best[t].j == -1) continue;
            if (best[t].d > worst){
                worst = best[t].d;
                pick = (int)t;
            }
        }

        if (pick == -1) break;

        //realizamos la asignacion
        int i = unassigned[pick];
        int j = best[pick].j;
        int rbd = idx2rbd[j];

        solution[i] = rbd;
        vac_by_rbd[rbd]--;
        total_vacantes_restantes--;

        unassigned[pick] = unassigned.back();
        unassigned.pop_back();

        //si el colegio se lleno, se rearma el conjunto de colegios
        if (vac_by_rbd[rbd] == 0)
            schools = rebuild();
    }

    //================= CALCULO FINAL =================
    //rearmamos el vector solucion / calculo de metricas

    //calcular cuantos asignados en cada preferencia
    std::vector<int> asignacion_por_pref(saParams.max_choices + 1, 0);

    for (int i = 0; i < saParams.n_students; ++i){

        int rbd_asignado = solution[i];

        if (rbd_asignado == -1){
            asignacion_por_pref[saParams.max_choices]++;
            continue;
        }

        const auto& A = alumnos[i];
        int pref_final = -1;

        for (int p = 0; p < A.num_ele; ++p){
            size_t idx = (size_t)A.choices[p] - 1;
            if (colegios[idx].rbd == rbd_asignado){
                pref_final = p;
                break;
            }
        }

        if (pref_final != -1)
            asignacion_por_pref[pref_final]++;
        else
            asignacion_por_pref[saParams.max_choices]++;
    }

    //colegios con vacantes
    int colegios_con_vac_final = 0;
    for (int j = 0; j < (int)colegios.size(); ++j){
        int rbd = idx2rbd[j];
        if (vac_by_rbd[rbd] > 0) ++colegios_con_vac_final;
    }
    
    //rearmado del vector solucion en rbd
    std::vector<int> sol_idx(solution.size());
    std::transform(solution.begin(), solution.end(), sol_idx.begin(),
        [&](int rbd){ auto it = rbd2idx.find(rbd); return rbd == -1 ? -1 : (it != rbd2idx.end() ? it->second : -1); });

    solution = sol_idx;    
    
    //mismas metricas que el calculo principal
    double var1 = meanDist(sol_idx.data(),distMat)/saParams.max_dist;
    double var2 = S(sol_idx.data(), alumnosSep, totalVuln);
    double var3 = costCupo(sol_idx.data(),cupoArray);
    double var4 = penaltyParents(sol_idx.data(), h_penalty_matrix)/(saParams.n_students);

    std::cout << "--------------- Resultados Algoritmo SAE ----------------" << "\n";
    std::cout << "Costo solucion*: " << var1*saParams.alpha1 + var2*saParams.alpha2 + var3*saParams.alpha3 + var4*saParams.alpha4 << "\n";
    std::cout << "distancia: " << var1 << "\n";
    std::cout << "Segregacion: " << var2 << "\n";
    std::cout << "CostoCupo: " << var3 << "\n";
    std::cout << "Penalty final: " << var4 << "\n";

    std::cout << "Asignados en preferencia numero ->";
    for (int p = 0; p < saParams.max_choices; ++p){
        std::cout << (p+1) << ":" << asignacion_por_pref[p] << "  ";
    }
    std::cout << "99: " << asignacion_por_pref[saParams.max_choices] << "\n";

    std::cout << "Colegios aun con vacantes: " << colegios_con_vac_final
              << " | Vacantes totales restantes: " << total_vacantes_restantes << "\n";

    std::cout << "--------------- Finalizo con exito ----------------" << "\n";

    std::ofstream out("asignacion_sae.txt");
    for (int i = 0; i < saParams.n_students; ++i) {
        out << i << "," << sol_idx[i] << "\n";
    }

    guardarResultadoSAE("info-registerSAE.txt",
                    var1, var2, var3, var4,
                    asignacion_por_pref,
                    saParams.max_choices);

}





//funcion obsoleta
int SimulatedAnnealing::balanceCostoCupo(
    int* currentSolution,
    const std::vector<Info_alu>& alumnos,
    const std::vector<Info_colegio>& colegios
){

    std::vector<int> capacidad(saParams.n_colegios);
    std::vector<int> ocupados(saParams.n_colegios, 0);

    for (int j = 0; j < saParams.n_colegios; ++j) {
        capacidad[j] = (int)std::floor(colegios[j].num_alu * 1.1);
    }

    for (int i = 0; i < saParams.n_students; ++i) {
        ocupados[currentSolution[i]]++;
    }


    std::unordered_map<int, std::vector<int>> alus_sobrecupo;
    std::unordered_map<int, int> vacantes_col;

    for (int j = 0; j < saParams.n_colegios; ++j) {
        vacantes_col[j] = capacidad[j] - ocupados[j];
        if (vacantes_col[j] < 0) { //colegio en sobrecupo
            int count = 0;
            for (int i = 0; i < saParams.n_students && count < fabs(vacantes_col[j]); ++i) {
                if (currentSolution[i] == j) {
                    alus_sobrecupo[j].push_back(i); //guardar ids de estudiantes
                    count++;
                }
            }
        }
    }

    //fase media heurisitica, robado del algoritmo del SAE
    //solucion temporal mientras veo como diseñar una funcion que penalice con el tiempo pero a un recocido simulado

    //heuristica por penalty
    for (auto& kv : alus_sobrecupo) {
    int colegio_sob = kv.first;
    auto& alumnos_vec = kv.second;

    //revolver el vector para dar aleatoriedad
    std::shuffle(alumnos_vec.begin(), alumnos_vec.end(), mt);

        for (int alu_id : alumnos_vec) {
            if (vacantes_col[colegio_sob] >= 0) { //mientras aun el colegio esta en sobrecupo ver si puede sacar alumnos
                break;
            }
            //recorrer preferencias del alumno
            for (int pref_col : alumnos[alu_id].choices) {
                if (vacantes_col[pref_col-1] > 0) {//choices del dataset va del 1 al 63, pero en el arreglo esta del 0 al 62
                    
                    //cout << "el alumno "<< alu_id+1 << " estaba en " <<colegio_sob+1<<" y tiene preferencia en: " << pref_col<< "\n";
                    currentSolution[alu_id] = pref_col-1;
                    vacantes_col[pref_col-1]--;   // ocupa vacante
                    vacantes_col[colegio_sob]++; // libera en sobrecupo
                    break; // ya lo movimos, no seguimos buscando
                }
            }
        }
    }


    //heuristica por distancia
    std::vector<int> sobrantes;
    sobrantes.reserve(1024);
    for (const auto& kv : alus_sobrecupo) {
        const auto& vec = kv.second;
        sobrantes.insert(sobrantes.end(), vec.begin(), vec.end());
    }

    //Helper: construir índices de colegios con vacantes (>0)
    auto rebuild_schools_with_vac = [&]() {
        std::vector<int> idx;
        idx.reserve(saParams.n_colegios);
        for (int j = 0; j < saParams.n_colegios; ++j)
            if (vacantes_col[j] > 0) idx.push_back(j);
        return idx;
    };
    std::vector<int> schools_with_vac_idx = rebuild_schools_with_vac();

    //Bucle principal
    while (!sobrantes.empty() && !schools_with_vac_idx.empty()) {
        // PRUNEA alumnos cuyo colegio de origen ya no está en sobrecupo
        // (balance >= 0). Esto evita mover alumnos innecesariamente.
        {
            size_t w = 0;
            for (size_t t = 0; t < sobrantes.size(); ++t) {
                int i = sobrantes[t];
                int col_origen = currentSolution[i];
                if (vacantes_col[col_origen] < 0) { // aún en sobrecupo
                    sobrantes[w++] = i;
                }
            }
            sobrantes.resize(w);
            if (sobrantes.empty()) break;
        }

        //para cada alumno sobrante, calcular su colegio más cercano con cupo
        struct Best { double d; int j; };
        std::vector<Best> best(sobrantes.size(), {std::numeric_limits<double>::infinity(), -1});
        int factibles = 0;

        for (size_t t = 0; t < sobrantes.size(); ++t) {
            int i = sobrantes[t];
            int col_origen = currentSolution[i];

            double bestd = std::numeric_limits<double>::infinity();
            int bestj = -1;

            for (int j : schools_with_vac_idx) {
                if (j == col_origen) continue;          //excluir el mismo colegio de origen por seguridad
                if (vacantes_col[j] <= 0) continue;     
                double d = distMat[i][j];
                if (d < bestd) { bestd = d; bestj = j; }
            }

            best[t] = {bestd, bestj};
            if (bestj != -1) ++factibles;
        }

        if (factibles == 0) break; // nadie tiene dónde ir

        //elegir el alumno con peor "mejor distancia" (mayor best.d)
        int pick_t = -1;
        double worst = -1.0;
        int tie_i = std::numeric_limits<int>::max();
        for (size_t t = 0; t < sobrantes.size(); ++t) {
            if (best[t].j == -1) continue;
            int i = sobrantes[t];
            if (best[t].d > worst || (best[t].d == worst && i < tie_i)) {
                worst = best[t].d;
                tie_i = i;
                pick_t = (int)t;
            }
        }
        if (pick_t == -1) break;

        //Reasignar ese alumno
        int i_pick   = sobrantes[pick_t];
        int j_pick   = best[pick_t].j;
        int col_orig = currentSolution[i_pick];

        //mover
        //cout << "el estudiante "<<i_pick<< " se movio a "<<j_pick<< "\n";
        currentSolution[i_pick] = j_pick;
        //actualizar balances
        vacantes_col[j_pick]--;      // consume un cupo en destino
        vacantes_col[col_orig]++;    // libera presión en origen

        //quitar de "sobrantes"
        sobrantes[pick_t] = sobrantes.back();
        sobrantes.pop_back();

        // 3.e) Si el destino se quedó sin vacantes, reconstruir la lista de colegios con cupo
        if (vacantes_col[j_pick] == 0) {
            schools_with_vac_idx = rebuild_schools_with_vac();
        }
        // PD: si algun colegio pasó de 0 a >0 (por liberar origen),
        // no es obligatorio reconstruir aquí; en la siguiente iteración
        // lo considerara el rebuild cuando algún colegio llegue a 0
    }

    double var1 = meanDist(bestSolution, distMat)/saParams.max_dist;
    double var2 = S(bestSolution, alumnosSep, totalVuln);
    double var3 = costCupo(bestSolution, cupoArray);
    double var4 = penaltyParents(bestSolution,h_penalty_matrix)/saParams.n_students;
    double costSolution21 = alpha[0]*var1 + alpha[1]*var2 + alpha[2]*var3 + alpha[3]*var4;
    cout << "----- Fase 2.1-----\n";
    cout << "Costo solucion: " << costSolution21 << "\n";
    cout << "distancia: " << var1 << "\n"; //lo normalice
    cout << "Segregación: " << var2 << "\n";
    cout << "CostoCupo: " << var3 << "\n";
    cout << "Penalty final: " << var4 << "\n";
    int* summaryPrefs = summaryPreferences(bestSolution, dataSet->students);
    summaryCostoCupo(currentSolution, colegios);

    return summaryPrefs[saParams.max_choices];
}









double SimulatedAnnealing::runCPU() {
    CUDAWrapper* cudaWrapper = new CUDAWrapper(cuParams, saParams, mt);
    inicializationValues(cudaWrapper);
    cudaWrapper->memInit(previousSolution,
        bestSolution,
        currentSolution,
        cupoArray,
        alumnosSep,
        totalVuln,
        aluxcol,
        aluVulxCol,
        matrestest,
        alpha,
        currentVars);
    
    int n_students   = saParams.n_students;
    int n_colegios   = saParams.n_colegios;
    DataResult solution;
    auto start = std::chrono::high_resolution_clock::now();

    cout << "--------------- Primeros datos -------------\n";
    cout << "Primer costo de solución: " << costBestSolution << "\n";
    cout << "Primer distancia: " << meanDist(currentSolution, distMat)/saParams.max_dist << "\n"; //lo normalize
    cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    cout << "Primer CostoCupo: " << costCupo(currentSolution, cupoArray) << "\n";
    cout << "Penalty inicial: " << penaltyParents(currentSolution,h_penalty_matrix)/saParams.n_students << "\n\n";

    std::ofstream tiempoCPU;
    tiempoCPU.open("../../save/tiempo_CPU.txt", std::ios::app);

    while(saParams.temp > saParams.min_temp){
        auto start_shuffle = std::chrono::high_resolution_clock::now();
        memcpy(currentSolution, previousSolution, n_students * sizeof(int));
        memcpy(aluxcol, previousAluxCol, n_colegios * sizeof(int));
        memcpy(aluVulxCol, previousAluVulxCol, n_colegios * sizeof(int));
        memcpy(currentVars, previousVars, 4 * sizeof(double));

        shuffle(saParams.shuffle_student, saParams.max_changes_students, dist);
        shuffle(saParams.shuffle_colegios, saParams.max_changes_school, dist2);

        auto end_shuffle = std::chrono::high_resolution_clock::now();
        double time_taken_shuffle = std::chrono::duration_cast<std::chrono::nanoseconds>(end_shuffle - start_shuffle).count();

        auto start_it = std::chrono::high_resolution_clock::now();

        int newSchool = saParams.shuffle_colegios[0];
        double cost_solution;
        double min_cost = 1.0;

        //verion CPU primer kernel, mover a todos los estudiantes a una misma escuela, uno a la vez, y comparar el mejor Z
        for (int i=0; i<saParams.n_students; i++){
            int aluchange = saParams.shuffle_student[i];
            int currentSchool= currentSolution[aluchange];
            cost_solution = 0.0;

            double sumDist     = currentVars[0];
            double totalSesc   = currentVars[1];
            double totalCupo   = currentVars[2];
            double penalty     = currentVars[3];

            ////////////////////////////////////////////////////////////////
            /////// Descuenta antes de mover
            ////////////////////////////////////////////////////////////////

            // Distancia
            sumDist -= distMat[aluchange][currentSchool];

            // Seg de la escuela actual
            int totA = aluxcol[currentSchool];
            int vulA = aluVulxCol[currentSchool];
            int novA = totA - vulA;
            totalSesc -= std::fabs((vulA/(double)totalVuln) - (novA/(double)(n_students - totalVuln)));

            // Costocupo escuela actual 
            double p_costCupo = (double)totA/cupoArray[currentSchool];
            totalCupo -= calcCostoCupo(p_costCupo);

            //std::cout << std::setprecision(17) << "CPU Vars: " << totA << " " << cupoArray[currentSchool] << "\n";
            // Seg de la escuela nueva
            int totB = aluxcol[newSchool];
            int vulB = aluVulxCol[newSchool];
            int novB = totB - vulB;
            totalSesc -= std::fabs( (vulB/(double)totalVuln) - (novB/(double)(n_students - totalVuln)) );
            
            //costocupo escuela nueva
            p_costCupo = (double)totB/cupoArray[newSchool];
            totalCupo -= calcCostoCupo(p_costCupo);

            ////////////////////////////////////////////////////////////////
            ////// Calculó despues de mover
            //////////////////////////////////////////////////////////////

            //distancia de la nueva escuela
            sumDist += distMat[aluchange][newSchool];

            //seg de la escuela actual
            int totA2 = aluxcol[currentSchool] -1;
            int vulA2 = aluVulxCol[currentSchool] -alumnosSep[aluchange];
            int novA2 = totA2 - vulA2;
            totalSesc += std::fabs( (vulA2/(double)totalVuln) - (novA2/(double)(n_students - totalVuln)) );

            //costocupo escuela actual
            p_costCupo = double(totA2)/cupoArray[currentSchool];
            totalCupo += calcCostoCupo(p_costCupo);

            //seg de la escuela antigua
            int totB2 = aluxcol[newSchool] +1;
            int vulB2 = aluVulxCol[newSchool] +alumnosSep[aluchange];
            int novB2 = totB2 - vulB2;
            totalSesc += std::fabs( (vulB2/(double)totalVuln) - (novB2/(double)(n_students - totalVuln)) );

            //costcupo escuela antigua
            p_costCupo = double(totB2)/cupoArray[newSchool];
            totalCupo += calcCostoCupo(p_costCupo);


            //Calculo de penalty
            //penalty de la escuela actual
            penalty   -= h_penalty_matrix[aluchange*saParams.n_colegios +currentSchool];

            //penalty de la escuela nueva
            penalty   += h_penalty_matrix[aluchange*saParams.n_colegios +newSchool];


            // Combinar con alphas (mismas normalizaciones que el kernel)
            cost_solution += alpha[0]*(sumDist / (double)n_students) / saParams.max_dist;
            cost_solution += alpha[1]*(totalSesc * 0.5);
            cost_solution += alpha[2]*(totalCupo / (double)n_colegios);
            cost_solution += alpha[3]*(penalty / n_students);

            if (cost_solution < min_cost){
                min_cost = cost_solution;
                solution.costSolution = cost_solution;
                solution.stu = aluchange;
                solution.col = newSchool;
            }


        }

        auto end_it = std::chrono::high_resolution_clock::now();
        double time_taken_it = std::chrono::duration_cast<std::chrono::nanoseconds>(end_it - start_it).count();

        auto start_onesolution = std::chrono::high_resolution_clock::now();

        //verion CPU del segundo kernel, realizar el mejor movimiento del paso anterior
        DataResult movimiento = cpu_one_tid_newSolution(solution.stu, solution.col);
        costCurrentSolution = movimiento.costSolution;

        auto end_onesolution = std::chrono::high_resolution_clock::now();
        double time_taken_onesolution = std::chrono::duration_cast<std::chrono::nanoseconds>(end_onesolution - start_onesolution).count();

        auto start_post = std::chrono::high_resolution_clock::now();

        //verificar error
        if(costCurrentSolution<0.00 || isnan(costCurrentSolution)){
            cout << "error" << endl;
            cout << saParams.count << endl;
            std::cout << saParams.shuffle_colegios[cuParams.selectThread] << "\n";
            std::cout << saParams.shuffle_student[cuParams.selectBlock] << "\n";
            std::cout << "distancia: " << meanDist(currentSolution,distMat) << "\n";
            std::cout << "Segregación: " << S(currentSolution,alumnosSep, totalVuln) << "\n";
            std::cout << "CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";
            std::cout << "Penalty: " << penaltyParents(currentSolution,h_penalty_matrix) << "\n";
            std::cout << costCurrentSolution;
            exit(1);
        }

        //copypaste del ciclo original
        if(costCurrentSolution < costBestSolution){
            //cudaWrapper->AcceptanceBestSolution();
            memcpy(bestSolution, currentSolution, n_students * sizeof(int));
            memcpy(previousSolution, currentSolution, n_students * sizeof(int));
            memcpy(previousAluxCol, aluxcol, n_colegios * sizeof(int));
            memcpy(previousAluVulxCol, aluVulxCol, n_colegios * sizeof(int));
            memcpy(previousVars, currentVars, 4 * sizeof(double));
            memcpy(bestVars, currentVars, 4 * sizeof(double));
            costBestSolution = costCurrentSolution;
            costPreviousSolution = costCurrentSolution;

            saParams.c_accepta++;
            saParams.count_rechaso = 0;
            //cout << costCurrentSolution << " | " << saParams.count <<  endl;

        }
        else {
            if(acceptanceCriterion->apply(costPreviousSolution,costCurrentSolution,dist_accepta ) == 1) {
                //cudaWrapper->AcceptanceSolution();
                memcpy(previousSolution, currentSolution, n_students * sizeof(int));
                memcpy(previousAluxCol, aluxcol, n_colegios * sizeof(int));
                memcpy(previousAluVulxCol, aluVulxCol, n_colegios * sizeof(int));
                memcpy(previousVars, currentVars, 4 * sizeof(double));
                costPreviousSolution = costCurrentSolution;

                saParams.count_rechaso = 0;
                saParams.c_accepta++;
            }
            else {
                saParams.count_rechaso++;
            }
        }

        //supongo que esto hace bajar la temperatura
        if(lengthTemperature->apply()){
            coolingScheme->apply();
        }
        
        reheatingMethod->apply();
        saParams.count_trials++;
        saParams.count++;
        
        auto end_post = std::chrono::high_resolution_clock::now();
        double time_taken_post = std::chrono::duration_cast<std::chrono::nanoseconds>(end_post - start_post).count();
        tiempoCPU << std::setprecision(10)<<time_taken_shuffle << "," << time_taken_it<< "," << time_taken_onesolution<< ","<<time_taken_post << "\n";

    }

    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    cout << "--------------- Resultado Final ----------------" << "\n";
    cout << "Numero de Ciclos: " << saParams.count << "\n";
    cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
    cout << "Costo de la mejor solución: " << costBestSolution << "\n";
    cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
    cout << "Tiempo de ejecución de SA: " << time_taken << "\n";
    cout << "distancia: " << meanDist(bestSolution, distMat)/saParams.max_dist << "\n"; //lo normalice
    cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << "CostoCupo: " << costCupo(bestSolution, cupoArray) << "\n";
    cout << "Penalty final: " << penaltyParents(bestSolution,h_penalty_matrix)/saParams.n_students << "\n";

    int* summaryPrefs = summaryPreferences(bestSolution, dataSet->students);
    int* data_costocupo = summaryCostoCupo(bestSolution, dataSet->colegios);
    cout << "--------------- Finalizo con exito ----------------" << "\n";
    
    return costBestSolution;
}


void SimulatedAnnealing::runSAE(){
    CUDAWrapper* cudaWrapper = new CUDAWrapper(cuParams, saParams, mt);
    inicializationValues(cudaWrapper);
    cudaWrapper->memInit(previousSolution,
        bestSolution,
        currentSolution,
        cupoArray,
        alumnosSep,
        totalVuln,
        aluxcol,
        aluVulxCol,
        matrestest,
        alpha,
        currentVars);
    

    std::vector<int> solution;
    asignacionSAE(dataSet->students, dataSet->colegios, solution);
    summaryCostoCupo(solution.data(), dataSet->colegios);

}



void SimulatedAnnealing::guardarResultadoSAE(const std::string& filename,
                         double var1, double var2, double var3, double var4,
                         const std::vector<int>& asignacion_por_pref,
                         int max_choices)
{
    bool archivo_existe = std::filesystem::exists(filename);

    std::ofstream out(filename, std::ios::app); // append mode
    if (!out.is_open()) {
        std::cerr << "Error al abrir " << filename << "\n";
        return;
    }

    // ----------------------------------------------------
    // Escribir encabezado SOLO la primera vez
    // ----------------------------------------------------
    if (!archivo_existe) {
        out << "seed,distancia,segregacion,costocupo,penalty,Z";
        for (int p = 0; p < max_choices; ++p) {
            out << ",pref_" << (p + 1);
        }
        out << ",unassigned\n";  // preferencia 99
    }

    // ----------------------------------------------------
    // Escribir los valores de una ejecución del algoritmo
    // ----------------------------------------------------
    double cost_solution = saParams.alpha1*var1 + saParams.alpha2*var2 + saParams.alpha3*var3 +saParams.alpha4*var4; 
    out << saParams.seed << "," << var1 << "," << var2 << "," << var3 << "," << var4 << "," << cost_solution;

    for (int p = 0; p < max_choices; ++p) {
        out << "," << asignacion_por_pref[p];
    }

    // El último índice es el de no asignados
    out << "," << asignacion_por_pref[max_choices] << "\n";
}