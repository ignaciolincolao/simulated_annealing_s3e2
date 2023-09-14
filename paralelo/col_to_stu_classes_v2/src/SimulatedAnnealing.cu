#include <SimulatedAnnealing.cuh>
#include <CUDAWrapper.cuh>

#include <limits>

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
    }
    
/*
static std::mutex addInfoMutex;

static void addInfoToSave(RecordManager *recordManager,
                          double costCurrentSolution,
                          double meanDist,
                          double S,
                          double costCupo,
                          SimulatedParams *saParams)
{
    std::lock_guard<std::mutex> lock(addInfoMutex);
    recordManager->vector_costCurrentSolution.push_back(costCurrentSolution);
    recordManager->vector_meanDist.push_back(meanDist);
    recordManager->vector_segregation.push_back(S);
    recordManager->vector_costoCupo.push_back(costCupo);
    recordManager->vector_temp.push_back(saParams->temp);
    recordManager->vector_count.push_back(saParams->count);
}
*/
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
    cout << "Primer distancia: " << meanDist(currentSolution, distMat) << "\n";
    cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    cout << "Primer CostoCupo: " << costCupo(currentSolution, cupoArray) << "\n\n";
#if SAVE_DATA
    recordManager->openRecordInfo();
    recordManager->openRecordGraphics();
    recordManager->openRecordGraphicsBestSolution();

    recordManager->SaveInfoInit(costBestSolution,
                                meanDist(currentSolution, distMat),
                                S(currentSolution, alumnosSep, totalVuln),
                                costCupo(currentSolution, cupoArray));

    recordManager->SaveGraphicsInit(meanDist(currentSolution, distMat),
                                    S(currentSolution, alumnosSep, totalVuln),
                                    costCupo(currentSolution, cupoArray),
                                    costCurrentSolution);

    recordManager->SaveGraphicsBestSolution(currentSolution);


    recordManager->closeRecordInfo();
    recordManager->closeRecordGraphics();
    recordManager->closeRecordGraphicsBestSolution();
#endif
    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////
    auto start = std::chrono::high_resolution_clock::now();
    ///////////////////////////////////////////////////
    /// Comienza a ejecutarse el algoritmo de SA
    ///////////////////////////////////////////////////
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
        

    
        ///////////////////////////////////////////////////
        ///  Metodo Nuevo
        //////////////////////////////////////////////////
        cudaWrapper->getCurrentSolutionGpuToHost(costCurrentSolution);  
        if(costCurrentSolution >= costPreviousSolution){
            if(acceptanceCriterionApply() == 1){
                cudaWrapper->newSolutionRandomSelection(dist,
                    dist2);
            }
        }
        ///////////////////////////////////////////////////
        ///  Actualiza la nueva solución en la GPU
        //////////////////////////////////////////////////
        cudaWrapper->newSolutionUpdate(costCurrentSolution);
            
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
        
        
        ///////////////////////////////////////////////////
        /// 
        //////////////////////////////////////////////////
        if(costCurrentSolution < costBestSolution){
            cudaWrapper->AcceptanceBestSolution();
            costBestSolution = costCurrentSolution;
            costPreviousSolution = costCurrentSolution;
            saParams.c_accepta++;
            saParams.count_rechaso = 0;

            // futures.push_back(std::async(std::launch::async,
            //                   addInfoToSave,
            //                   recordManager,
            //                   costCurrentSolution,
            //                   meanDist(currentSolution, distMat),
            //                   S(currentSolution, alumnosSep, totalVuln),
            //                   costCupo(currentSolution, cupoArray),
            //                   &saParams
            //                 ));
#if SAVE_DATA
            cudaWrapper->copySolutionToHost(bestSolution, previousSolution);
            recordManager->vector_costCurrentSolution.emplace_back(costBestSolution);
            recordManager->vector_meanDist.emplace_back(meanDist(bestSolution, distMat));
            recordManager->vector_segregation.emplace_back(S(bestSolution, alumnosSep, totalVuln));
            recordManager->vector_costoCupo.emplace_back(costCupo(bestSolution, cupoArray));
            recordManager->vector_temp.emplace_back(saParams.temp);
            recordManager->vector_count.emplace_back(saParams.count);
#endif
        }
        else {
            if(acceptanceCriterion->apply(costPreviousSolution,costCurrentSolution,dist_accepta ) == 1) {

                cudaWrapper->AcceptanceSolution();
                costPreviousSolution = costCurrentSolution;

                saParams.count_rechaso = 0;
                saParams.c_accepta++;
            }
            else {
                saParams.count_rechaso++;
                
            }
        }

        if(lengthTemperature->apply()){
            coolingScheme->apply();
        }
        reheatingMethod->apply();
        

        
        //cout << costCurrentSolution << costPreviousSolution << "| |" << saParams.temp << "| |" << saParams.count<< endl;
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
    cout << "--------------- Finalizo con exito ----------------" << "\n";

#if SAVE_DATA
    recordManager->openRecordInfo();
    recordManager->openRecordGraphics();
    recordManager->openRecordGraphicsBestSolution();
    recordManager->openRecordRegister();

    recordManager->SaveInfoFinish(costPreviousSolution,
                                  costBestSolution,
                                  costCurrentSolution,
                                  time_taken,
                                  meanDist(bestSolution, distMat),
                                  S(bestSolution, alumnosSep, totalVuln),
                                  costCupo(bestSolution, cupoArray));

    recordManager->SaveGraphicsFinish();

    recordManager->SaveGraphicsBestSolution(bestSolution);
    recordManager->SaveInfoRegister(
        time_taken,
        costBestSolution,
        meanDist(bestSolution, distMat),
        S(bestSolution, alumnosSep, totalVuln),
        costCupo(bestSolution, cupoArray),
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
        cuParams.n_thread
    );

    recordManager->closeRecordInfo();
    recordManager->closeRecordGraphics();
    recordManager->closeRecordGraphicsBestSolution();
    recordManager->closeRecordRegister();
#endif

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
    costBestSolution = calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;

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
    previousVars[0] = currentVars[0];
    previousVars[1] = currentVars[1];
    previousVars[2] = currentVars[2];
    
    double var1,var2,var3;
    var1 = (currentVars[0]/saParams.n_students);
    var1= (var1/saParams.max_dist);
    //cout << var1 << "\n";
    var2 = (currentVars[1]/2.0);
    //cout << var2 << "\n";
    var3 = (currentVars[2] /saParams.n_colegios);
    costBestSolution = (double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;
    auto start_compare = std::chrono::high_resolution_clock::now();
    auto end_compare = std::chrono::high_resolution_clock::now();
    double time_taken_v1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_compare - start_compare).count();
    //double time_taken_v2;
    //double vector_time1 =0;
    // double vector_time2=0;
    saParams.max_changes_students = min(cuParams.n_thread*cuParams.n_block, saParams.n_students);
    saParams.max_changes_school = min(cuParams.n_block, saParams.n_colegios);

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
    return (double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
}

///////////////////////////////////////////////////
/// Distancia promedio que recorreran los estudiantes
///////////////////////////////////////////////////
double SimulatedAnnealing::meanDist(const int *currentSolution, double  **distMat){
    double sumDist=0.0;
    for(int i=0;i<saParams.n_students;i++){
        sumDist+=round_n(distMat[i][currentSolution[i]]); // distMat[estudiante][escuela]
    }
    //cout << "meanDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << saParams.n_students << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist/saParams.n_students;
}

double SimulatedAnnealing::sumDist(const int *currentSolution, double  **distMat){
    double sumDist=0.0;
    for(int i=0;i<saParams.n_students;i++){
        sumDist+=round_n(distMat[i][currentSolution[i]]); // distMat[estudiante][escuela]
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
            totalSesc+=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(saParams.n_students-totalVuln))));
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
            totalSesc+=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(saParams.n_students-totalVuln))));
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
    // double a = 0.0;
    for(int j=0;j<saParams.n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<saParams.n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        totalcostCupo+=round_n((double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2)));
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
        totalcostCupo+= round_n((double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2)));
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
    for(int x=0; x<3; x++){
        sumaAlpha +=alpha[x];
    }
    for(int x=0; x<3; x++){
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


SimulatedAnnealing::~SimulatedAnnealing() {

}

