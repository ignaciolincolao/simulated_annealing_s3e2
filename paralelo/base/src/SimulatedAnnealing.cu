#include <SimulatedAnnealing.cuh>
#include <CUDAWrapper.cuh>

#include <limits>
#include <iomanip> //para imprimir mas decimales

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
        choices_parents,
        currentVars);

    cout << "--------------- Primeros datos -------------\n";
    cout << "Primer costo de solución: " << costBestSolution << "\n";
    cout << "Primer distancia: " << meanDist(currentSolution, distMat) << "\n";
    cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    cout << "Primer CostoCupo: " << costCupo(currentSolution, cupoArray) << "\n";
    cout << "Penalty inicial: " << penaltyParents(currentSolution)/saParams.p_weight << "\n\n";
#if SAVE_DATA
    #if ENABLE_OPEN_RECORD_INFO
    recordManager->openRecordInfo();
    recordManager->SaveInfoInit(costBestSolution,
        meanDist(currentSolution, distMat),
        S(currentSolution, alumnosSep, totalVuln),
        costCupo(currentSolution, cupoArray),
        penaltyParents(currentSolution));
    recordManager->closeRecordInfo();
    #endif
    #ifdef ENABLE_OPEN_RECORD_GRAPHICS
    recordManager->openRecordGraphics();
    recordManager->SaveGraphicsInit(meanDist(currentSolution, distMat),
    S(currentSolution, alumnosSep, totalVuln),
    costCupo(currentSolution, cupoArray),
    costCurrentSolution,
    penaltyParents(currentSolution));
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
    while(saParams.temp > saParams.min_temp){
        ///////////////////////////////////////////////////
        /// Copia Solución Anterior a la actual
        ///////////////////////////////////////////////////
        cudaWrapper->memCopyPrevToCurrent();
        ///////////////////////////////////////////////////
        ///  Selecciona aleatoria mente a los alumnos
        ///////////////////////////////////////////////////
        shuffle(saParams.shuffle_student, saParams.max_changes_students, dist);
        shuffle(saParams.shuffle_colegios, saParams.max_changes_school, dist2);
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
            //cout << costCurrentSolution << " | " << saParams.count << " | " << id_select <<  endl;

#if SAVE_DATA
    #ifdef ENABLE_OPEN_RECORD_REGISTER       
            cudaWrapper->copySolutionToHost(bestSolution, previousSolution);
            recordManager->vector_costCurrentSolution.emplace_back(costBestSolution);
            recordManager->vector_meanDist.emplace_back(meanDist(bestSolution, distMat));
            recordManager->vector_segregation.emplace_back(S(bestSolution, alumnosSep, totalVuln));
            recordManager->vector_costoCupo.emplace_back(costCupo(bestSolution, cupoArray));
            recordManager->vector_penalty.emplace_back(penaltyParents(bestSolution));
            recordManager->vector_temp.emplace_back(saParams.temp);
            recordManager->vector_count.emplace_back(saParams.count);
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
    cout << "distancia: " << meanDist(bestSolution, distMat) << "\n";
    cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << "CostoCupo: " << costCupo(bestSolution, cupoArray) << "\n";
    cout << "Penalty final: " << penaltyParents(bestSolution)/saParams.p_weight << "\n";
    cout << "--------------- Finalizo con exito ----------------" << "\n";
    

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
        penaltyParents(bestSolution));
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
        penaltyParents(bestSolution),
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
        bestSolution
    );
    recordManager->closeRecordRegister();
    #endif
    #ifdef ENABLE_OPEN_RECORD_MOVE_SOLUTION
    recordManager->openRecordMoveSolution();
    recordManager->AllMovementFinish();
    recordManager->closeRecordMoveSolution();
    #endif
    #ifdef ENABLE_OPEN_RECORD_SIMCE_UPDATE
    recordManager->openRecordInfoSimce();
    recordManager->simceScoreUpdate(bestSolution, dataSet->ptr_students,dataSet->ptr_colegios);
    recordManager->closeRecordInfoSimce();
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
    saParams.p_weight = saParams.n_students*500000.0;
    //cout << fixed << setprecision(70) << endl;
    //srand(time(NULL));


    ///////////////////////////////////////////////////
    /// Inicializa Variables y arreglos
    ///////////////////////////////////////////////////

    aluxcol= (int *)malloc(sizeof(int)*saParams.n_colegios);
    aluVulxCol = (int *)malloc(sizeof(int)*saParams.n_colegios);
    choices_parents = (uint8_t *)malloc(5 * saParams.n_students);
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

    //borrable???
    //creo que aqui se calculo de forma normal, pero mas adelante se sobreescribe el valor con uno normalizado
    //entonces calCosto no es la version CPU
    costBestSolution = calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;

    //fin borrable ------

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
    currentVars[3] = penaltyParents(currentSolution);
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
    //es 50000.0 o 500000.0, sigo analizando estaba en 50000.0
    var4 = currentVars[3] /(saParams.n_students*500000.0);
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

    //creo que al oscar le falto incluir una division "/(saParams.n_students*50000.0)"
    //no tiene sentido todo lo demas es como promedio y metio un total que rompe todo, mas encima porque metio un 50000 hardcode
    //que significa
    //ahora entiendo porque piden tanto comentar los codigos 💀💀
    //me imagino que lo que intento hacer fue normalizar para que siempre de valores entre 0 y 1
    //nota mental agregar ese 50000 a sa.params (estaba en 50000.0)
    double var4 = penaltyParents(currentSolution)/(saParams.n_students*500000.0);
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
        totalcostCupo+=(double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2));
    }
    return totalcostCupo/saParams.n_colegios;
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
        totalcostCupo+= (double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2));
    }
    return totalcostCupo;
}
///////////////////////////////////////////////////
/// Genera una nueva solución en donde asigna a un estudiante a una escuela
/// aleatoriamente
///////////////////////////////////////////////////

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

        cupoArray[x] = ptr_colegios->num_alu+ ((int)((ptr_colegios->num_alu*10)/100));
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
        //agregar las preferencias de los padres (creo que esto es unicamente para la prueba CPU)
        for (std::size_t i = 0; i < 5; i++){
            choices_parents[x * 5 + i] = students[x].choices[i];
        }
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


//version CPU con la misma logica de calPenalty del GPU
//double calcPenaltyCPU(int currentCollege, const uint8_t choices[5]) {
//    double weights[6] = {500000, 0, 100, 200, 300, 400};
//    uint8_t index = 0;
//    for (size_t i = 1; i < 6; i++) {
//        index += (currentCollege == choices[i - 1]) * i;
//    }
//    return weights[index];
//}

/*
void SimulatedAnnealing::simceScoreUpdate(int *bestSolution, Info_alu *ptr_students, Info_colegio *ptr_colegios) {

    std::vector<std::array<double, 4>> schoolScores(saParams.n_colegios, {0.0, 0.0, 0.0, 0.0});
    //vector que almacena los nuevos puntajes simce de cada colegio (rbd, pmat, plen, numero de estudiantes en el colegio)

    for (int x = 0; x < saParams.n_students; x++) {
        int schoolID = bestSolution[x];
        double simceMath = ptr_students[x].pmat;
        double simceLanguage = ptr_students[x].plen;

        schoolScores[schoolID][0] = ptr_colegios[schoolID].rbd; 
        schoolScores[schoolID][1] += simceMath;                //guardamos la suma de los puntajes mat y len
        schoolScores[schoolID][2] += simceLanguage;             
        schoolScores[schoolID][3] += 1;                        //y el numero de alumnos en el colegio
    }

    std::vector<std::array<double, 3>> results; // Almacena {RBD, delta_math, delta_language}
    //vector que almacenara la informacion util (RBD del colegio, variacion en puntaje mat, variacion len)

    for (int i = 0; i < saParams.n_colegios; i++) {
        const auto& score = schoolScores[i];
        if (score[3] > 0) { //descartamos los que se quedaron sin alumnos matriculados despues del cambio
            double avgMath = score[1] / score[3]; //calculamos promedios
            double avgLanguage = score[2] / score[3];

            //calcular variacion entre antes y depsues del cambio
            double deltaMath = avgMath - ptr_colegios[i].pmat;
            double deltaLanguage = avgLanguage - ptr_colegios[i].plen;

            results.push_back({score[0], deltaMath, deltaLanguage}); //utilizamos el RBD real del colegio
            //PD: pusieron que en bestsolution apunta a la posicion en el arreglo de info_col en lugar del rbd real.
        }
    }
    //guardamos los resultados como un txt
    //no me funciono lo de las flag
    const std::string outputFile = "../../save/resultados_simce.csv";
    //PD: el compilador no se que le pasa que no quiere crear carpetas y no reconoce mkdir
    //estoy creando la carpeta a mano, y tira todo a la carpeta de release/debug
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: No se pudo crear o abrir el archivo: " << outputFile << "\n";
        return;
    }

    outFile << "rbd,delta_pmat,delta_plen\n";
    

    for (const auto& res : results) {
        outFile << static_cast<int>(res[0]) << "," // rbd
                << res[1] << ","                   // delta_pmat
                << res[2] << "\n";                 // delta_plen
    }
    outFile.close();
    std::cout << "Archivo de actualizacion de puntajes guardado en: " << outputFile << "\n";

    //analisis:
    //se observa que elimina a la gran cantidad de los colegios, esto es porque en el grafico de R se observa que
    //las preferencias se concentran en un grupo cerrado de colegios, y hay colegios que no tienen ninguna preferencia
    //ademas como al ubicar a un individuo en un colegio que no selecciono (dandole un valor de 5000), el sistema prefiere
    //matricular a un individuo con sobrecupo respecto a un colegio que no sea del agrado de los padres.

    //no he tocado la funcion de penalty del oscar
}
*/

void SimulatedAnnealing::ValidateGPU(){
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
        choices_parents,
        currentVars);

    //error max detectado de 2e-16 (ciclo 54)
    int n_ciclos = 12000;
    std::vector<DataResult> sol;
    int id_select= 0;
    DataResult cpu;
    double costTempSol;
    double dif;
    double max_error = -1;
    int ciclo_max_error = -1;

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
        //cout << "accion: -> stu= " << r.stu << " -> col= " << r.col << " costo=" << r.costSolution << "\n";
        

        //2nd kernel (se cambia a un alu de col)
        cudaWrapper->newSolutionUpdate(costCurrentSolution, id_select);
        
        //copypaste del original
        if(costCurrentSolution < costBestSolution){
                cudaWrapper->AcceptanceBestSolution();
                cudaWrapper->copySolutionToHost(bestSolution, previousSolution);
                costTempSol = calCosto(bestSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
                costBestSolution = costCurrentSolution;
                costPreviousSolution = costCurrentSolution;
                saParams.c_accepta++;
                saParams.count_rechaso = 0;

                cpu = cpu_one_tid_newSolution(r.stu); //funcion CPU que imita al kernel

                dif = costTempSol - costBestSolution;
                if (dif>0){
                    cout << "iteracion N " << ciclo+1 << " : Resultado no es igual a CPU, error de "<< dif<< "\n";
                    cout << "accion realizada: -> stu= " << r.stu << " -> col= " << r.col << " costo=" << r.costSolution << "\n";
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
    cout << setprecision(16) << "Total CPU: " << costBestSol << "\n";
    cout << setprecision(16) << "Total GPU: " << costBestSolution << "\n";
    cout << setprecision(16) << "Error maximo: " << max_error << " ciclo: " << ciclo_max_error << "\n";



    //comprobacion del primer kernel
    /*
    int mismatches = 0;
    double eps = 1e-9;
    const int N = (int)sol.size();
    for (int tid = 0; tid < N; ++tid) {
        DataResult cpu = cpu_one_tid_newSolution(tid);
        const DataResult& gpu = sol[tid];

        bool ok_stu = (cpu.stu == gpu.stu);
        bool ok_col = (cpu.col == gpu.col);

        bool ok_cost;
        const double SENT = static_cast<double>(0xffffffffffffffffULL);
        if (cpu.costSolution >= 1e18 || gpu.costSolution >= 1e18) {
            ok_cost = ( (cpu.costSolution >= 1e18) && (gpu.costSolution >= 1e18) );
        } else {
            ok_cost = std::fabs(cpu.costSolution - gpu.costSolution) <= eps;
        }

        if (!(ok_stu && ok_col && ok_cost)) {
            if (mismatches < 10) {
                std::cout << "Mismatch tid " << tid
                          << " CPU(stu=" << cpu.stu << ", col=" << cpu.col
                          << ", cost=" << cpu.costSolution << ") vs "
                          << "GPU(stu=" << gpu.stu << ", col=" << gpu.col
                          << ", cost=" << gpu.costSolution << ")\n";
            }
            ++mismatches;
        }
    }

    std::cout << "[validateGPU_newSolution] "
              << (mismatches ? "Fallos=" : "OK. fallos=")
              << mismatches << " de " << N << "\n";

    */
    delete cudaWrapper;
}


//imitacion de kernel 1 en CPU (robado de GPT)

static inline double calcPenaltyCPU(int currentCollege, const uint8_t ch[5]) {
    const double w[6] = {500000, 0, 100, 200, 300, 400};
    uint8_t idx = 0;
    for (int i = 1; i <= 5; ++i) idx += (currentCollege == ch[i-1]) * i;
    return w[idx];
}

// Mismo cálculo que newSolution_kernel, pero leyendo los miembros de 'this'
DataResult SimulatedAnnealing::cpu_one_tid_newSolution(int tid) const {
    const int n_students   = saParams.n_students;
    const int n_colegios   = saParams.n_colegios;
    const int aluchange    = tid;
    const int newSchool    = saParams.shuffle_colegios[0]; // igual que el kernel actual
    //cout << tid << " " << aluchange << " "<< saParams.shuffle_colegios[0] << "\n";
    const int currentSchool= currentSolution[aluchange];

    DataResult out; // ctor pone cost=0, col=-1, stu=-1
    out.stu = aluchange;
    out.col = newSchool;

    // Sentinela si no hay movimiento (igual que kernel)
    if (newSchool == currentSchool) {
        out.costSolution = static_cast<double>(0xffffffffffffffffULL);
        return out;
    }

    // Snapshot de agregados
    double sumDist     = currentVars[0];
    double totalSesc   = currentVars[1];
    double totalCupo   = currentVars[2];
    double penalty     = currentVars[3];
    //std::cout << std::setprecision(17) << "CPU: " << currentVars[0] << " " << currentVars[1] << " "<< currentVars[2] << " "<< currentVars[3] << "\n";

    // choices del alumno
    uint8_t ch[5] = {
        choices_parents[aluchange*5+0],
        choices_parents[aluchange*5+1],
        choices_parents[aluchange*5+2],
        choices_parents[aluchange*5+3],
        choices_parents[aluchange*5+4],
    };

    // Matriz de distancias en HOST es contigua: stride = n_colegios
    const size_t stride = static_cast<size_t>(n_colegios);

    // --- Descontar antes de mover ---
    sumDist -= matrestest[static_cast<size_t>(aluchange) * stride + currentSchool];

    int totA = aluxcol[currentSchool];
    int vulA = aluVulxCol[currentSchool];
    int novA = totA - vulA;
    totalSesc -= std::fabs( (vulA/(double)totalVuln) - (novA/(double)(n_students - totalVuln)) );
    penalty   -= calcPenaltyCPU(currentSchool, ch);
    totalCupo -= (double)totA * std::fabs((double)cupoArray[currentSchool] - totA)
                 / std::pow(cupoArray[currentSchool]*0.5, 2);

    //std::cout << std::setprecision(17) << "CPU Vars: " << totA << " " << cupoArray[currentSchool] << "\n";

    int totB = aluxcol[newSchool];
    int vulB = aluVulxCol[newSchool];
    int novB = totB - vulB;
    totalSesc -= std::fabs( (vulB/(double)totalVuln) - (novB/(double)(n_students - totalVuln)) );
    
    totalCupo -= (double)totB * std::fabs((double)cupoArray[newSchool] - totB)
                 / std::pow(cupoArray[newSchool]*0.5, 2);

    //std::cout << std::setprecision(17) << "CPU cost 1: " << totalCupo << "\n";
    penalty   += calcPenaltyCPU(newSchool, ch);


    //ELimina el estudiante de la escuela actual
    aluxcol[currentSchool]-=1;
    aluVulxCol[currentSchool]-=alumnosSep[aluchange];
    //Asigna al estudiante a la nueva escuela
    currentSolution[aluchange] = newSchool;
    aluxcol[newSchool]+=1;
    aluVulxCol[newSchool]+=alumnosSep[aluchange];

    // --- Sumar después de mover ---
    sumDist += matrestest[static_cast<size_t>(aluchange) * stride + newSchool];
    int totA2 = aluxcol[currentSchool];
    int vulA2 = aluVulxCol[currentSchool];
    int novA2 = totA2 - vulA2;
    totalSesc += std::fabs( (vulA2/(double)totalVuln) - (novA2/(double)(n_students - totalVuln)) );

    //std::cout << std::setprecision(17) << "CPU Vars 2: " << totA2 << " " << cupoArray[currentSchool] << "\n";
    totalCupo += (double)totA2 * std::fabs((double)cupoArray[currentSchool] - totA2)
                 / std::pow(cupoArray[currentSchool]*0.5, 2);



    int totB2 = aluxcol[newSchool];
    int vulB2 = aluVulxCol[newSchool];
    int novB2 = totB2 - vulB2;
    totalSesc += std::fabs( (vulB2/(double)totalVuln) - (novB2/(double)(n_students - totalVuln)) );

    //std::cout << std::setprecision(17) << "CPU cost 2 A: " << totalCupo << "\n";
    totalCupo += (double)totB2 * std::fabs((double)cupoArray[newSchool] - totB2)
                 / std::pow(cupoArray[newSchool]*0.5, 2);

    //std::cout << std::setprecision(17) << "CPU cost 2 B: " << totalCupo << "\n";
    currentVars[0] = sumDist;
    currentVars[1] = totalSesc;
    currentVars[2] = totalCupo;
    currentVars[3] = penalty;

    //std::cout << std::setprecision(17) << "CPU 2: " << currentVars[0] << " " << currentVars[1] << " "<< currentVars[2] << " "<< currentVars[3] << "\n";


    // Combinar con alphas (mismas normalizaciones que el kernel)
    const double weight_n_students = n_students * 500000.0;     // igual que memInit
    const double var1 = (sumDist / (double)n_students) / saParams.max_dist;
    const double var2 = (totalSesc * 0.5);
    const double var3 = (totalCupo / (double)n_colegios);
    const double var4 = (penalty / weight_n_students);

    out.costSolution = alpha[0]*var1 + alpha[1]*var2 + alpha[2]*var3 + alpha[3]*var4;
    return out;
}

