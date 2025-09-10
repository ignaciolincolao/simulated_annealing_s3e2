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
            recordManager->vector_penalty.emplace_back(penaltyParents(bestSolution,h_penalty_matrix));
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
    cout << "Penalty final: " << penaltyParents(bestSolution,h_penalty_matrix)/saParams.n_students << "\n";
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
        bestSolution
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
    //es 50000.0 o 500000.0, sigo analizando estaba en 50000.0
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

/*
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
        //printf("Penalty %f, alu: %d, col: %d", h_penalty_matrix[i*saParams.n_colegios +currentSolution[i]], i, currentSolution[i]);
        penalty += h_penalty_matrix[i*saParams.n_colegios +currentSolution[i]];
    }

    return penalty;
}


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
    

    std::vector<int> num_ele_vec;
    num_ele_vec.reserve(dataSet->students.size());

    for (const auto& alu : dataSet->students) {
        num_ele_vec.push_back(alu.num_ele);
    }

    std::vector<int> choices_matrix;
    choices_matrix.reserve(saParams.n_students * saParams.max_choices);

for (const auto& alu : dataSet->students) {
    // solo tomo los primeros saParams.max_choices elementos
    for (int j = 0; j < saParams.max_choices; ++j) {
        choices_matrix.push_back(static_cast<int>(alu.choices[j]));
    }
}

    std::vector<float> penalty_matrix(saParams.n_students * saParams.n_colegios);

    cudaWrapper->compute_penalty_matrix(choices_matrix.data(),
                            num_ele_vec.data(),
                            penalty_matrix.data());
    


    //error max detectado de 2e-16 (ciclo 54)
    int n_ciclos = 10;
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
        //cout << "accion: -> stu= " << r.stu << " -> col= " << r.col << " costo=" << r.costSolution << "\n";
        
        //2nd kernel (se cambia a un alu de col)
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

                cpu = cpu_one_tid_newSolution(r.stu); //funcion CPU que imita al kernel

                dif = costTempSol - costBestSolution;
                if (dif>0){
                    if (dif > epsilon){
                        //cout << "iteracion N " << ciclo+1 << " : Resultado no es igual a CPU, error de "<< dif<< "\n";
                        //cout << "accion realizada: -> stu= " << r.stu << " -> col= " << r.col << " costo=" << r.costSolution << "\n";  
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













/**
 * FUNCIÓN HOST MEJORADA: Gestión completa con parámetros flexibles
 *
 * Esta función te permite experimentar fácilmente con diferentes configuraciones
 * del modelo sin recompilar el kernel. Piensa en ella como el "panel de control"
 * de tu sistema de penalizaciones.
 *
 * PARÁMETROS NUEVOS Y SUS EFECTOS:
 * @param max_pref_penalty: Controla la discontinuidad entre preferencias vs no-preferencias
 *                          - Valores bajos (0.2-0.4): Gran salto, enfatiza importancia de estar en la lista
 *                          - Valores altos (0.6-0.8): Salto moderado, más "perdón" para asignaciones fuera de preferencias
 * @param alpha: Controla la curvatura de preferencias dentro de la lista
 *               - Valores bajos (0.5): Preferencias más "planas", diferencia sutil entre 1ra y última opción
 *               - Valores altos (2.0): Preferencias muy jerarquizadas, primeras opciones mucho mejor valoradas

void SimulatedAnnealing::compute_penalty_matrix(
    int* h_preferences_matrix,        // Entrada: matriz de preferencias en CPU
    int* h_num_preferences,           // Entrada: número de preferencias por estudiante
    float* h_penalty_matrix,          // Salida: matriz de penalizaciones en CPU
    int num_students,
    int num_schools,
    int max_preferences_per_student,  // Típicamente igual a num_schools, pero puede ser menor
    float alpha,               // Curvatura exponencial (recomendado: 0.5 - 2.0)
    float max_pref_penalty    // Penalización máxima para preferencias (recomendado: 0.3 - 0.8)
) {
   
    
   
    // Transferir datos de entrada a GPU
    cudaMemcpy(d_preferences, h_preferences_matrix, prefs_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_num_preferences, h_num_preferences, counts_size, cudaMemcpyHostToDevice);
   
    // Configuración de grid: balance entre paralelismo y cache efficiency
    // 16×16 = 256 threads por block es óptimo para la mayoría de GPUs modernas
    dim3 block_size(16, 16);    
    dim3 grid_size(
        (num_schools + block_size.x - 1) / block_size.x,      // Ceil division para cobertura completa
        (num_students + block_size.y - 1) / block_size.y
    );
   
    // Información sobre configuración de ejecución
    int total_blocks = grid_size.x * grid_size.y;
    int total_threads = total_blocks * block_size.x * block_size.y;
    printf("Configuración de ejecución: %d blocks, %d threads total\n", total_blocks, total_threads);
   
    // Ejecutar kernel con parámetros flexibles
    compute_preference_penalty_matrix<<<grid_size, block_size>>>(
        d_preferences, d_num_preferences, d_penalty_matrix,
        num_students, num_schools, max_preferences_per_student,
        alpha, max_pref_penalty  // Los nuevos parámetros configurables
    );
   
    // Sincronizar y verificar errores detalladamente
    error = cudaDeviceSynchronize();
    if (error != cudaSuccess) {
        printf("Error ejecutando kernel: %s\n", cudaGetErrorString(error));
       
        // Cleanup en caso de error
        cudaFree(d_preferences);
        cudaFree(d_num_preferences);
        cudaFree(d_penalty_matrix);
        return;
    }
   
    cudaMemcpy(h_penalty_matrix, d_penalty_matrix, matrix_size, cudaMemcpyDeviceToHost);

    printf("penalty:  %f %f %f \n", h_penalty_matrix[num_schools*1+60], h_penalty_matrix[num_schools*1+59], h_penalty_matrix[num_schools*1+58]);
    printf("Kernel ejecutado exitosamente. Matriz de penalizaciones calculada.\n");

    // === NUEVO: escribir matriz a archivo txt ===
    std::ofstream fout("penalty_matrix.txt");
    if (!fout) {
        std::cerr << "Error: no se pudo abrir penalty_matrix.txt para escribir\n";
    } else {
        fout << std::fixed << std::setprecision(4);
        for (int i = 0; i < num_students; ++i) {
            for (int j = 0; j < num_schools; ++j) {
                fout << h_penalty_matrix[i * num_schools + j];
                if (j < num_schools - 1) fout << ",";
            }
            fout << "\n";  // salto de línea por estudiante
        }
    }
    fout.close();
    std::cout << "Matriz de penalizaciones escrita en penalty_matrix.txt\n";

    // Cleanup: liberar toda la memoria GPU
    cudaFree(d_preferences);
    cudaFree(d_num_preferences);
    cudaFree(d_penalty_matrix);
}

/**
 * UTILIDAD MEJORADA: Validación comprehensiva con diagnósticos detallados
 *
 * Esta función no solo verifica que los resultados sean correctos, sino que también
 * te ayuda a entender cómo se están comportando tus parámetros en la práctica.
 * Es como un "chequeo médico" para tu matriz de penalizaciones.
 *
 * @param penalty_matrix: La matriz calculada por el kernel
 * @param preferences_matrix: Datos originales para cross-validation
 * @param num_preferences: Número de preferencias por estudiante
 * @param num_students, num_schools: Dimensiones
 * @param max_pref_penalty: El valor máximo usado en el cálculo
 * @param verbose: Si imprimir estadísticas detalladas
 *
 * @return true si todos los valores son válidos, false si hay problemas
 
bool validate_penalty_matrix(
    float* penalty_matrix,
    int* preferences_matrix,
    int* num_preferences,
    int num_students,
    int num_schools,
    int max_preferences_per_student,
    float max_pref_penalty,
    bool verbose = false
) {
    int invalid_count = 0;
    int preference_penalties = 0;      // Conteo de penalizaciones para colegios preferidos
    int non_preference_penalties = 0;  // Conteo de penalizaciones para colegios NO preferidos
    float sum_pref_penalties = 0.0f;   // Para calcular promedio de penalizaciones de preferencias
   
    for (int student = 0; student < num_students; student++) {
        int student_prefs = num_preferences[student];
       
        for (int school = 0; school < num_schools; school++) {
            float penalty = penalty_matrix[student * num_schools + school];
           
            // Determinar si este colegio está en las preferencias de este estudiante
            bool is_preferred = false;
            for (int p = 0; p < student_prefs; p++) {
                if (preferences_matrix[student * max_preferences_per_student + p] == school) {
                    is_preferred = true;
                    break;
                }
            }
           
            if (is_preferred) {
                // Verificar rango válido para preferencias: [0.0, max_pref_penalty]
                if (penalty < 0.0f || penalty > max_pref_penalty + 1e-5f) {  // Pequeña tolerancia para precision
                    invalid_count++;
                    if (verbose) {
                        printf("Error: Estudiante %d, Colegio %d (preferido): penalización %.6f fuera de rango [0, %.3f]\n",
                               student, school, penalty, max_pref_penalty);
                    }
                } else {
                    preference_penalties++;
                    sum_pref_penalties += penalty;
                }
            } else {
                // Verificar valor exacto para no-preferencias: debe ser 1.0
                if (fabsf(penalty - 1.0f) > 1e-5f) {  // Tolerancia para float precision
                    invalid_count++;
                    if (verbose) {
                        printf("Error: Estudiante %d, Colegio %d (NO preferido): penalización %.6f, debería ser 1.0\n",
                               student, school, penalty);
                    }
                } else {
                    non_preference_penalties++;
                }
            }
        }
    }
   
    // Calcular estadísticas útiles para entender el comportamiento del modelo
    if (verbose) {
        float avg_pref_penalty = (preference_penalties > 0) ?
                                 sum_pref_penalties / preference_penalties : 0.0f;
       
        printf("\n=== ESTADÍSTICAS DE VALIDACIÓN ===\n");
        printf("Total de asignaciones evaluadas: %d\n", num_students * num_schools);
        printf("Penalizaciones para colegios preferidos: %d (promedio: %.4f)\n",
               preference_penalties, avg_pref_penalty);
        printf("Penalizaciones para colegios NO preferidos: %d\n", non_preference_penalties);
        printf("Errores encontrados: %d\n", invalid_count);
       
        // Análisis de la distribución de preferencias por estudiante
        int min_prefs = num_schools + 1, max_prefs = 0;
        float sum_prefs = 0.0f;
        for (int i = 0; i < num_students; i++) {
            int prefs = num_preferences[i];
            min_prefs = (prefs < min_prefs) ? prefs : min_prefs;
            max_prefs = (prefs > max_prefs) ? prefs : max_prefs;
            sum_prefs += prefs;
        }
        float avg_prefs = sum_prefs / num_students;
       
        printf("\n=== ANÁLISIS DE PREFERENCIAS ===\n");
        printf("Preferencias por estudiante - Mínimo: %d, Máximo: %d, Promedio: %.2f\n",
               min_prefs, max_prefs, avg_prefs);
        printf("Esto significa que el %.1f%% de las evaluaciones son para colegios preferidos\n",
               100.0f * preference_penalties / (num_students * num_schools));
       
        if (avg_prefs / num_schools < 0.3f) {
            printf("INSIGHT: Los estudiantes son selectivos (promedio %.1f preferencias de %d colegios)\n",
                   avg_prefs, num_schools);
            printf("         La discontinuidad en %.3f vs 1.0 será muy efectiva\n", max_pref_penalty);
        } else {
            printf("INSIGHT: Los estudiantes son flexibles (promedio %.1f preferencias de %d colegios)\n",
                   avg_prefs, num_schools);
            printf("         Considera ajustar max_pref_penalty para mayor diferenciación\n");
        }
    }
   
    if (invalid_count > 0) {
        printf("Advertencia: %d valores de penalización fuera de los rangos esperados\n", invalid_count);
        return false;
    }
   
    return true;
}
*/