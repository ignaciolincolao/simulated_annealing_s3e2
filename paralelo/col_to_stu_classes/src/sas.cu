#include <sas-old.cuh>
#include <kernel.cuh>
#include <ExplorationCriterion.hpp>
#include <AcceptanceCriterion.hpp>
#include <TemperatureLength.hpp>
#include <ReheatingMethods.hpp>
#include <CoolingScheme.hpp>

#include <limits>

#define DECIMAL 16

#include <assert.h>
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
///////////////////////////////////////////////////
/// Variables constantes CUDA
///////////////////////////////////////////////////



typedef std::numeric_limits<double> dbl;

///////////////////////////////////////////////////
/// Variables globales.
///////////////////////////////////////////////////




double sasFunc() {
    int x = 0, z = 0;
    int totalVuln = 0;
    cout.precision(dbl::max_digits10);
    //cout << fixed << setprecision(70) << endl;
    //srand(time(NULL));
    ///////////////////////////////////////////////////
    /// Genera archivo de almacenamiento de datos
    ///////////////////////////////////////////////////

    /*
    * Prepara archivos para guardar los datos
    */
    
    ofstream info;
    string infotxt = ruta_save + prefijo_save +"-info.txt"; 
    info.open(infotxt);
    /*
    * 
    */
    ofstream info_test;
    string nameinfo_test = ruta_save + prefijo_save+"-info-test.txt"; 
    info_test.open(nameinfo_test);
    /*
    * Genera los archivos que contienen información de los estados de estudiantes y escuelas durante
    * la ejecución del algoritmo
    */
    ofstream info_graficos;
    string name_info_graficos = ruta_save + prefijo_save +"-info-graficos.txt";
    info_graficos.open(name_info_graficos);
    ///////////////////////////////////////////////////
    /// Datos colegios
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura Info_colegio
    ///////////////////////////////////////////////////
    Info_colegio *ptr_colegios;
    vector<Info_colegio> colegios;
    getDataSchool(colegios);
    ptr_colegios = colegios.data();
    n_colegios = colegios.size();

    ///////////////////////////////////////////////////
    /// Datos Alumnos
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura info_student
    ///////////////////////////////////////////////////
    Info_alu *ptr_students;
    vector<Info_alu> students;
    getDataStudents(students,totalVuln);
    ptr_students = students.data();
    n_students = students.size();

    ///////////////////////////////////////////////////
    /// Inicializa Variables y arreglos
    ///////////////////////////////////////////////////


    int aluVulxCol[n_colegios], aluxcol[n_colegios];
    int previousAluxCol[n_colegios];
    int previousAluVulxCol[n_colegios];
    int bestAluxCol[n_colegios];
    int bestAluVulxCol[n_colegios];

    int *previousSolution = nullptr;
    int *bestSolution = nullptr;
    int *currentSolution = nullptr;
    double **distMat = nullptr;
    int *cupoArray = nullptr;
    int *alumnosSep = nullptr;

    
    double  costBestSolution,
        costPreviousSolution,
        costCurrentSolution,
        *ptr_alpha = &alpha[0];
    
    int count = 0;

    cudaMallocHost((void**)&previousSolution, sizeof(int)*n_students);
    cudaMallocHost((void**)&bestSolution, sizeof(int)*n_students);
    cudaMallocHost((void**)&currentSolution, sizeof(int)*n_students);
    cudaMallocHost((void**)&cupoArray, sizeof(int)*n_colegios);
    cudaMallocHost((void**)&alumnosSep, sizeof(int)*n_students);
    /*
    previousSolution = (int *)malloc(sizeof(int)*n_students);
    bestSolution=(int *)malloc(sizeof(int)*n_students);
    currentSolution=(int *)malloc(sizeof(int)*n_students);
    cupoArray=(int *)malloc(sizeof(int)*n_colegios);
    alumnosSep = (int *)malloc( sizeof(int)*n_students);
    */
    distMat=(double **)malloc(sizeof(double)*n_students);
    for(x=0; x < n_students; x++) {
        distMat[ x ]=(double *)malloc(sizeof(double)*n_colegios);
    }

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
                students,
                colegios);
    assignSchoolToArray(previousSolution, bestSolution, currentSolution, ptr_colegios, ptr_students, cupoArray);
    calcDist(ptr_colegios, ptr_students, distMat);
    max_dist = getMaxDistance(distMat);
    normalizedAlpha(alpha);

    ///////////////////////////////////////////////////
    /// Registro de datos
    ///////////////////////////////////////////////////
    costBestSolution = calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;
    
    
    cout << "--------------- Primeros datos -------------" << "\n";
    cout << "Primer costo de solución: " << costBestSolution << "\n";
    cout << "Primer distancia: " << meanDist(currentSolution,distMat) << "\n";
    cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    cout << "Primer CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";

    info << "--------------- Primeros datos -------------" << "\n";
    info << "Primer costo de solución: " << costBestSolution << "\n";
    info << "Primer distancia: " << meanDist(currentSolution,distMat) << "\n";
    info << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    info << "Primer CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";


    ///////////////////////////////////////////////////
    /// Generación de archivos que almacenan información de los graficos
    ///////////////////////////////////////////////////

    info_graficos << setprecision(13);
    info_graficos << count << "," 
                << meanDist(currentSolution,distMat)/max_dist << "," // Distancia promedio recorrida por los estudiantes normalizada
                << meanDist(currentSolution,distMat) << "," // Distancia promedio recorrida por los estudiantes
                << S(currentSolution, alumnosSep, totalVuln) << "," // Indice de duncan
                << costCupo(currentSolution,cupoArray) << "," // Costo cupo de las escuelas
                << costCurrentSolution << "," // Solución actual
                << temp << setprecision(13) << "\n"; // Temperatura actual

    count++;
    ///////////////////////////////////////////////////
    /// Genera arreglos que contendran valores del 0 hasta n_students y n_colegios
    ///////////////////////////////////////////////////

    int *shuffle_student, *shuffle_colegios;
    cudaMallocHost((void**)&shuffle_student, sizeof(int)*n_students);
    cudaMallocHost((void**)&shuffle_colegios, sizeof(int)*n_colegios);
    for (int i = 0; i < n_students; i++) {
        shuffle_student[i] = i;
    }
    for (int i=0; i < n_colegios; i++){
        shuffle_colegios[i]=i;
    }
    ///////////////////////////////////////////////////
    /// Posicion estudiantes
    ///////////////////////////////////////////////////

    ofstream info_graficos_bestSolution;
    string name_info_graficos_bestSolution = ruta_save + prefijo_save +"-info-graficos_bestSolution.txt"; // concatenar
    info_graficos_bestSolution.open(name_info_graficos_bestSolution);
    for(x = 0; x < n_students; x++){
        info_graficos_bestSolution << currentSolution[x] << ",";
    }
    info_graficos_bestSolution << "\n";

    ///////////////////////////////////////////////////
    /// Genera distribuciones para seleccionar un estudiante y una escuela al azar
    ///////////////////////////////////////////////////

    dist = uniform_int_distribution<int>(0, n_students-1);
    dist2 = uniform_int_distribution<int>(0, n_colegios-1);

    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////
    auto start = std::chrono::high_resolution_clock::now();
    ///////////////////////////////////////////////////
    /// Comienza a ejecutarse el algoritmo de SA
    ///////////////////////////////////////////////////


    vector<double> vector_costCurrentSolution;
    vector<double> vector_meanDist;
    vector<double> vector_segregation;
    vector<double> vector_costoCupo;
    vector<double> vector_temp;
    vector<int> vector_count;


    std::vector<double> vector_historyCostSolution;
    std::vector<double> vector_historyTemp;
    std::vector<double> vector_historymeanDist;
    std::vector<double> vector_historymeanDistNorm;
    std::vector<double> vector_historySegregation;
    std::vector<double> vector_historycostoCupo;
    std::vector<bool> vector_historyAcceptSolution;
    std::vector<int> vector_historyAsign;
    std::vector<std::tuple <int,int>> vector_historyMove;


    
    int count_rechaso=0;
    // int reheating = 0;
    int c_accepta = 0;
    int c_cooling_temperature = 0;
    // int valmaxheating=n_colegios;
    // int count_reheating = 0;
    // double bestTemp = 0;
    double k_reheating_init = k_reheating;
    double temp_init = temp;
    int count_trials = 0;
    float len1_init = len1;
    float len2_init = len2;
    double len3_init = len3;
    double len4_init = len4;

    ////////////////////////////////////////////////////////////////////////
    // VARIABLES DE PRUEBA
    ////////////////////////////////////////////////////////////////////////

    // double costCurrentSolutionV2 = costCurrentSolution;
    double *currentVars;
    cudaMallocHost( (void**)&currentVars,3 * sizeof(double));
    double *previousVars;
    cudaMallocHost( (void**)&previousVars,3 * sizeof(double)); 
    double *bestVars;
    cudaMallocHost( (void**)&bestVars,3 * sizeof(double)); 

    currentVars[0] = sumDist(currentSolution,distMat);
    currentVars[1] = sumS(currentSolution, alumnosSep, totalVuln);
    currentVars[2] = sumCostCupo(currentSolution,cupoArray);
    previousVars[0] = currentVars[0];
    previousVars[1] = currentVars[1];
    previousVars[2] = currentVars[2];
    double var1,var2,var3;
    cout << costBestSolution << endl;
    var1 = (currentVars[0]/n_students);
    var1= (var1/max_dist);
    //cout << var1 << "\n";
    var2 = (currentVars[1]/2.0);
    //cout << var2 << "\n";
    var3 = (currentVars[2] /n_colegios);
    costBestSolution = (double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
    cout << costBestSolution << endl;
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;
    auto start_compare = std::chrono::high_resolution_clock::now();
    auto end_compare = std::chrono::high_resolution_clock::now();
    double time_taken_v1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_compare - start_compare).count();
    //double time_taken_v2;
    double vector_time1 =0;
    // double vector_time2=0;

    /////////////////////////////////////////////////////////////////////////
    // VARIABLES PARA CUDA
    ////////////////////////////////////////////////////////////////////////

    double *d_distMat; /// clon de matriz de distancia
    int *d_currentSolution, *d_bestSolution, *d_previousSolution;
    int *d_alumnosSep; // Array que contendra a los estudiantes vulnerables
    int *d_cupoArray;
    ///////////////
    double *d_array_current_Solution;
    int *d_array_current_Solution_alu;
    int *d_array_current_Solution_col;
    ///////////////
    int *d_aluxcol,*d_previousAluxcol;
    int *d_aluVulxCol,*d_previousAluVulxCol;
    int *d_shuffle_students;
    int *d_shuffle_colegios;
    double *d_currentVars, *d_bestVars, *d_previousVars;
    double *d_costPreviousSolution, *d_costBestSolution, *d_costCurrentSolution;

    int max_changes_students = min(n_thread*n_block, n_students);
    int max_changes_school = min(n_block, n_colegios);



    cudaMalloc((void **) &d_array_current_Solution, n_block * sizeof(double));
    cudaMalloc((void **) &d_costCurrentSolution, 1 * sizeof(double));
    cudaMalloc((void **) &d_costBestSolution, 1 * sizeof(double));
    cudaMalloc((void **) &d_costPreviousSolution, 1 * sizeof(double));
    cudaMalloc((void **) &d_currentVars, 3 * sizeof(double));
    cudaMalloc((void **) &d_bestVars, 3 * sizeof(double));
    cudaMalloc((void **) &d_previousVars, 3 * sizeof(double));
    cudaMalloc((void **) &d_array_current_Solution_alu, n_block * sizeof(int)); 
    cudaMalloc((void **) &d_array_current_Solution_col, n_block * sizeof(int));
    cudaMalloc((void **) &d_shuffle_colegios, max_changes_school  * sizeof(int));
    cudaMalloc((void **) &d_shuffle_students, max_changes_students * sizeof(int));
    cudaMalloc((void **) &d_aluxcol,n_colegios * sizeof(int));
    cudaMalloc((void **) &d_previousAluxcol,n_colegios * sizeof(int));
    cudaMalloc((void **) &d_aluVulxCol,n_colegios * sizeof(int));
    cudaMalloc((void **) &d_previousAluVulxCol,n_colegios * sizeof(int));
    cudaMalloc((void **) &d_currentSolution, n_students * sizeof(int));  // Solución actual
    cudaMalloc((void **) &d_bestSolution, n_students * sizeof(int));
    cudaMalloc((void **) &d_previousSolution, n_students * sizeof(int));
    cudaMalloc((void **) &d_alumnosSep, n_students * sizeof(int)); // arreglo que contiene la id de cada usuario vulnerable
    cudaMalloc((void **) &d_cupoArray, n_colegios * sizeof(int));


    double *matrestest;
    cudaMallocHost( (void**)&matrestest,sizeof(double) * n_students * n_colegios); 
    double *array_costCurrentSolution = (double *) malloc(sizeof(double) * n_block * n_thread);
    for (x = 0; x < n_students; x++) {
        for (z = 0; z < n_colegios; z++) {
            matrestest[n_colegios * x + z] = distMat[x][z];
        }
    }
    for (x = 0; x < n_block; x++){
        for (z = 0; z < n_thread; z++){
            array_costCurrentSolution[n_thread * x + z] = 0.0;
        }
    }



    ////////////////////////////////////////////////////
    /////// Stream 
    ///////////////////////////////////////////////////
    int deviceId;
    int numberOfSMs;
    cudaDeviceProp deviceProp;
    cudaGetDevice(&deviceId);
    cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId); // Calcula el numero de SMstream 
    cudaGetDeviceProperties(&deviceProp, 0);
    int threadsPerBlock = 256;
    int numberOfBlocks = 32 * numberOfSMs;
    int NUM_STREAMS = 10;
    int nWarp = deviceProp.warpSize;
    cudaStream_t streams[NUM_STREAMS];
    for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamCreate(&streams[i]); }

   

    ///////////////////////////////////////////////////
    /// Valores que nunca van a cambiar
    //////////////////////////////////////////////////////




    size_t pitch;
    cudaMallocPitch(&d_distMat,
                    &pitch,
                    n_colegios * sizeof(double),
                    n_students); // Reserva memoria para la matriz de distancia


    gpuErrchk( cudaMemcpyToSymbolAsync( d_alpha, alpha,  3 * sizeof(double),0,cudaMemcpyHostToDevice,streams[2]));
    

    gpuErrchk( cudaMemcpyToSymbolAsync( d_n_students, &n_students, sizeof(int),0,cudaMemcpyHostToDevice,streams[3]));
    gpuErrchk( cudaMemcpyToSymbolAsync( d_n_colegios, &n_colegios, sizeof(int),0,cudaMemcpyHostToDevice,streams[4]));
    gpuErrchk( cudaMemcpyToSymbolAsync( d_max_dist, &max_dist, sizeof(double),0,cudaMemcpyHostToDevice,streams[5]));
    gpuErrchk( cudaMemcpyToSymbolAsync( d_totalVuln, &totalVuln, sizeof(int),0,cudaMemcpyHostToDevice,streams[6]));


    size_t h_pitchBytes = n_colegios * sizeof(double);
    cudaMemcpy2DAsync(d_distMat,
                 pitch,
                 matrestest,
                 h_pitchBytes,
                 n_colegios * sizeof(double),
                 n_students,
                 cudaMemcpyHostToDevice,
                 streams[3]);


    ///////////////////////////////////////////////////
    /// Inicializa las distribuciónes
    ///////////////////////////////////////////////////

    dist = std::uniform_int_distribution<int>(0, n_students-1);
    dist2 = std::uniform_int_distribution<int>(0, n_colegios-1);

    ///////////////////////////////////////////////////
    /// Contador de tiempo de ejecución en cuda
    ///////////////////////////////////////////////////

    cudaEvent_t start_cuda;
    cudaEvent_t stop_cuda;
    cudaEventCreate(&start_cuda);
    cudaEventCreate(&stop_cuda);
    // float elapsedTime;
    // double timeCuda = 0.0;
    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////
    //int vef_count = 0;


    cudaMemcpyAsync(d_currentSolution, currentSolution, n_students * sizeof(int), cudaMemcpyHostToDevice,streams[2]);
    cudaMemcpyAsync(d_previousSolution, currentSolution, n_students * sizeof(int), cudaMemcpyHostToDevice,streams[3]);
    cudaMemcpyAsync(d_bestSolution, currentSolution, n_students * sizeof(int), cudaMemcpyHostToDevice,streams[4]);
    cudaMemcpyAsync(d_aluxcol, aluxcol, n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[5]);
    cudaMemcpyAsync(d_previousAluxcol, aluxcol, n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[6]);
    cudaMemcpyAsync(d_aluVulxCol, aluVulxCol, n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[7]);
    cudaMemcpyAsync(d_previousAluVulxCol, aluVulxCol, n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[8]);
    cudaMemcpyAsync(d_currentVars, currentVars, 3 * sizeof(double), cudaMemcpyHostToDevice,streams[9]);
    cudaMemcpyAsync(d_previousVars, currentVars, 3 * sizeof(double), cudaMemcpyHostToDevice,streams[0]);
    cudaMemcpyAsync(d_bestVars, currentVars, 3 * sizeof(double), cudaMemcpyHostToDevice,streams[1]);
    cudaMemcpyAsync(d_alumnosSep, alumnosSep, n_students * sizeof(int), cudaMemcpyHostToDevice,streams[2]);
    cudaMemcpyAsync(d_cupoArray, cupoArray, n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[3]);

    cudaError_t errSync  = cudaGetLastError();
    cudaError_t errAsync = cudaDeviceSynchronize();
    
    if (errSync != cudaSuccess) 
        printf("0 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("0 Async kernel error: %s\n", cudaGetErrorString(errAsync));
    ///////////////////////////// Incorporar para acceder mas rapido al costCurrentSolution
    //int deviceId;
    //cudaGetDevice(&deviceId);                                         // The ID of the currently active GPU device.
    //cudaMemPrefetchAsync(pointerToSomeUMData, size, deviceId); 

    CoolingScheme cooling = CoolingScheme(&temp, coolingRate);
    //Reheating reheating = Reheating(&temp, &k_reheating, &n_reheating);

    while(cooling.getTemp() > min_temp){

        copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[0]>>>(d_currentSolution, d_previousSolution,n_students);
        copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[1]>>>(d_aluxcol, d_previousAluxcol,n_colegios);
        copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[2]>>>(d_aluVulxCol, d_previousAluVulxCol,n_colegios);
        copyVars<<<1,3,0,streams[3]>>>(d_currentVars, d_previousVars);
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("1 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("1 Async kernel error: %s\n", cudaGetErrorString(errAsync));
        //for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamSynchronize(streams[i]); }

        /*
        memcpy(currentSolution,previousSolution,sizeof(int)*n_students);
        memcpy(aluxcol,previousAluxCol,sizeof(int)*n_colegios);
        memcpy(aluVulxCol,previousAluVulxCol,sizeof(int)*n_colegios);
        memcpy(currentVars,previousVars,sizeof(double)*3);
        */

        ///////////////////////////////////////////////////
        ///  Selecciona aleatoria mente a los alumnos
        ///////////////////////////////////////////////////

        shuffle(shuffle_student, max_changes_students, dist);
        shuffle(shuffle_colegios, max_changes_school, dist2);
        ///////////////////////////////////////////////////
        /// Actualiza la memoria en CUDA
        ///////////////////////////////////////////////////
        /*
        cudaMemcpy(d_currentSolution, currentSolution, n_students * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_aluxcol, aluxcol, n_colegios * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_aluVulxCol, aluVulxCol, n_colegios * sizeof(int), cudaMemcpyHostToDevice);
        */



        ///////////////////////////////////////////////////
        ///  Envia datos a GPU
        ///////////////////////////////////////////////////

        
        cudaMemcpyAsync(d_shuffle_students, shuffle_student, max_changes_students* sizeof(int), cudaMemcpyHostToDevice,streams[0]);
        cudaMemcpyAsync(d_shuffle_colegios, shuffle_colegios, max_changes_school * sizeof(int), cudaMemcpyHostToDevice,streams[1]);
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("2 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("2 Async kernel error: %s\n", cudaGetErrorString(errAsync));

        ///////////////////////////////////////////////////
        ///  Ejecuta los kernel
        //////////////////////////////////////////////////
        //cout << (n_block/nWarp+1) << endl;
        newSolution_kernel<<<n_block,n_thread,
        (n_thread/nWarp+1) * sizeof(double)+ (n_thread/nWarp+1)* sizeof(int) + (n_thread/nWarp+1)* sizeof(int)>>>(
                        d_array_current_Solution,
                                d_array_current_Solution_alu,
                                d_array_current_Solution_col,
                                d_cupoArray,
                                d_alumnosSep,
                                d_aluxcol,
                                d_aluVulxCol,
                                d_currentSolution,
                                d_distMat,
                                d_shuffle_students,
                                d_shuffle_colegios,
                                d_currentVars,
                                pitch);
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("3 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("3 Async kernel error: %s\n", cudaGetErrorString(errAsync));
        reduce_block_kernel<<<1,n_block,
        (n_block/nWarp+1)* sizeof(double)+ (n_block/nWarp+1)* sizeof(int)+ (n_block/nWarp+1)* sizeof(int)>>>(d_array_current_Solution,
                d_array_current_Solution_alu,
                d_array_current_Solution_col);
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("4 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("4 Async kernel error: %s\n", cudaGetErrorString(errAsync));

        /********************************
        /* Metodo Nuevo
        */
        cudaMemcpy(&costCurrentSolution,&d_array_current_Solution[0], sizeof(double),cudaMemcpyDeviceToHost);
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("5 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("5 Async kernel error: %s\n", cudaGetErrorString(errAsync));

        if(costCurrentSolution >= costPreviousSolution){
            if(metropolisAC1(costPreviousSolution,costCurrentSolution)==1){
                selectThread = dist(mt);
                selectBlock = dist2(mt);
                cudaMemcpy(&d_array_current_Solution_alu[0],&selectThread, sizeof(int),cudaMemcpyHostToDevice);
                cudaMemcpy(d_array_current_Solution_col,&selectBlock, sizeof(int),cudaMemcpyHostToDevice);
                errSync  = cudaGetLastError();
                errAsync = cudaDeviceSynchronize();
                if (errSync != cudaSuccess) 
                printf("6 Sync kernel error: %s\n", cudaGetErrorString(errSync));
                if (errAsync != cudaSuccess)
                printf("6 Async kernel error: %s\n", cudaGetErrorString(errAsync));
                //cout << "son iguales" << endl;
            }
        }

        /* best
        if(costCurrentSolution >= costPreviousSolution){
            if(costCurrentSolution > costPreviousSolution){
                selectThread = dist(mt);
                selectBlock = dist2(mt);
                cudaMemcpy(&d_array_current_Solution_alu[0],&selectThread, sizeof(int),cudaMemcpyHostToDevice);
                cudaMemcpy(d_array_current_Solution_col,&selectBlock, sizeof(int),cudaMemcpyHostToDevice);
                errSync  = cudaGetLastError();
                errAsync = cudaDeviceSynchronize();
                if (errSync != cudaSuccess) 
                printf("6 Sync kernel error: %s\n", cudaGetErrorString(errSync));
                if (errAsync != cudaSuccess)
                printf("6 Async kernel error: %s\n", cudaGetErrorString(errAsync));
                //cout << "son iguales" << endl;
                }
            else{
                if(metropolisAC1(costPreviousSolution,costCurrentSolution)==1){
                    reduce_block_max<<<1,n_block,
                    (n_block/nWarp+1)* sizeof(double)+ (n_block/nWarp+1)* sizeof(int)+ (n_block/nWarp+1)* sizeof(int)>>>(d_array_current_Solution,
                    d_array_current_Solution_alu,
                    d_array_current_Solution_col);
                    errSync  = cudaGetLastError();
                    errAsync = cudaDeviceSynchronize();
                    if (errSync != cudaSuccess) 
                    printf("6 Sync kernel error: %s\n", cudaGetErrorString(errSync));
                    if (errAsync != cudaSuccess)
                    printf("6 Async kernel error: %s\n", cudaGetErrorString(errAsync));
                }
            }
            
        }
        */
        calculateSolution<<<1,1>>>(d_array_current_Solution,
            d_array_current_Solution_alu,
            d_array_current_Solution_col,
            d_cupoArray,
            d_alumnosSep,
            d_aluxcol,
            d_aluVulxCol,
            d_currentSolution,
            d_distMat,
            pitch,
            d_currentVars,
            d_costCurrentSolution);
        cudaMemcpy(&costCurrentSolution,&d_array_current_Solution[0], sizeof(double),cudaMemcpyDeviceToHost);
        //cudaMemcpyAsync(&selectThread,&d_array_current_Solution_alu[0], sizeof(int),cudaMemcpyDeviceToHost,streams[1]);
        //cudaMemcpyAsync(&selectBlock,d_array_current_Solution_col, sizeof(int),cudaMemcpyDeviceToHost,streams[2]);
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("8 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("8 Async kernel error: %s\n", cudaGetErrorString(errAsync));

        
        /********************************
        // Metodo antiguo
        */
        //cout << endl;
        /*
        calculateSolution<<<1,1>>>(d_array_current_Solution,
                    d_array_current_Solution_alu,
                    d_array_current_Solution_col,
                    d_cupoArray,
                    d_alumnosSep,
                    d_aluxcol,
                    d_aluVulxCol,
                    d_currentSolution,
                    d_distMat,
                    pitch,
                    d_currentVars,
                    d_costCurrentSolution);
        cudaMemcpy(&costCurrentSolution,&d_array_current_Solution[0], sizeof(double),cudaMemcpyDeviceToHost);
        //cudaMemcpyAsync(&selectThread,&d_array_current_Solution_alu[0], sizeof(int),cudaMemcpyDeviceToHost,streams[1]);
        //cudaMemcpyAsync(&selectBlock,d_array_current_Solution_col, sizeof(int),cudaMemcpyDeviceToHost,streams[2]);
        cudaDeviceSynchronize();
        */
        /*
        if(costCurrentSolution > costBestSolution){
            if(metropolisAC1(costPreviousSolution,costCurrentSolution)==1){
            reduce_block_max<<<1,n_block,
            (n_block/nWarp+1)* sizeof(double)+ (n_block/nWarp+1)* sizeof(int)+ (n_block/nWarp+1)* sizeof(int)>>>(d_array_current_Solution,
            d_array_current_Solution_alu,
            d_array_current_Solution_col);
            calculateSolution<<<1,1>>>(d_array_current_Solution,
                d_array_current_Solution_alu,
                d_array_current_Solution_col,
                d_cupoArray,
                d_alumnosSep,
                d_aluxcol,
                d_aluVulxCol,
                d_currentSolution,
                d_distMat,
                pitch,
                d_currentVars,
                d_costCurrentSolution);
            cudaMemcpy(&costCurrentSolution,&d_array_current_Solution[0], sizeof(double),cudaMemcpyDeviceToHost);
            //cudaMemcpyAsync(&selectThread,&d_array_current_Solution_alu[0], sizeof(int),cudaMemcpyDeviceToHost,streams[1]);
            //cudaMemcpyAsync(&selectBlock,d_array_current_Solution_col, sizeof(int),cudaMemcpyDeviceToHost,streams[2]);
            cudaDeviceSynchronize();
            }
        }*/
        
        //exit(0);
        ///////////////////////////////////////////////////
        ///  Actualizo datos basicos
        ///////////////////////////////////////////////////
    
        /*
        vector_historyAsign.push_back(currentSolution[shuffle_student[selectBlock]]);           
        aluxcol[currentSolution[shuffle_student[selectBlock]]]-=1; ///
        aluVulxCol[currentSolution[shuffle_student[selectBlock]]]-=alumnosSep[shuffle_student[selectBlock]]; ///
        aluxcol[shuffle_colegios[selectThread]]+=1; ///
        aluVulxCol[shuffle_colegios[selectThread]]+=alumnosSep[shuffle_student[selectBlock]]; ///
        currentSolution[shuffle_student[selectBlock]] = shuffle_colegios[selectThread]; ///
        */

        ///////////////////////////////////////////////////
        /// Salida en caso de error
        ///////////////////////////////////////////////////
        //std::cout << costCurrentSolution << "\n";
        //std::cout << selectThread << "\n";
        //std::cout << selectBlock << "\n";
        
        if(costCurrentSolution<0.00 || isnan(costCurrentSolution)){
            std::cout << shuffle_colegios[selectThread] << "\n";
            std::cout << shuffle_student[selectBlock] << "\n";
            std::cout << "distancia: " << meanDist(currentSolution,distMat) << "\n";
            std::cout << "Segregación: " << S(currentSolution,alumnosSep, totalVuln) << "\n";
            std::cout << "CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";
            std::cout << costCurrentSolution;
            exit(1);
        }
        
        
        
        if(costCurrentSolution < costBestSolution){
            
            copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[0]>>>(d_bestSolution, d_currentSolution,n_students);
            copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[1]>>>(d_previousSolution, d_currentSolution,n_students);
            copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[2]>>>(d_previousAluxcol, d_aluxcol,n_colegios);
            copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[3]>>>(d_previousAluVulxCol, d_aluVulxCol,n_colegios);
            copyVars<<<1,3,0,streams[4]>>>(d_previousVars, d_currentVars);
            copyVars<<<1,3,0,streams[5]>>>(d_bestVars, d_currentVars);
            copyCost<<<1,1,0,streams[6]>>>(d_costBestSolution,d_costCurrentSolution);
            copyCost<<<1,1,0,streams[7]>>>(d_costPreviousSolution,d_costCurrentSolution);
            //for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamSynchronize(streams[i]); }
            errSync  = cudaGetLastError();
            errAsync = cudaDeviceSynchronize();
            if (errSync != cudaSuccess) 
            printf("9 Sync kernel error: %s\n", cudaGetErrorString(errSync));
            if (errAsync != cudaSuccess)
            printf("9 Async kernel error: %s\n", cudaGetErrorString(errAsync));
            /*
            memcpy(bestSolution,currentSolution,sizeof(int)*n_students);
            memcpy(previousSolution,currentSolution,sizeof(int)*n_students);
            memcpy(previousAluxCol,aluxcol,sizeof(int)*n_colegios);
            memcpy(previousAluVulxCol,aluVulxCol,sizeof(int)*n_colegios);
            memcpy(previousVars,currentVars,sizeof(double)*3);
            memcpy(bestVars,currentVars,sizeof(double)*3);
            */
            
            costBestSolution = costCurrentSolution;
            costPreviousSolution = costCurrentSolution;
            //cout << costBestSolution << "| |" << temp << "| |" << count<< endl;
            /*
            vector_costCurrentSolution.push_back(costCurrentSolution);
            vector_meanDist.push_back(meanDist(currentSolution,distMat));
            vector_segregation.push_back(S(currentSolution, alumnosSep, totalVuln));
            vector_costoCupo.push_back(costCupo(currentSolution,cupoArray));
            vector_temp.push_back(temp);
            vector_count.push_back(count);
            */
            c_accepta++;
            count_rechaso = 0;
        }
        // En el caso que el la solución actual sea mas alta intenta aceptar una peor solución en base
        // a la función acepta
        else{
            if(metropolisAC1(costPreviousSolution,costCurrentSolution) == 1) {

                copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[0]>>>(d_previousSolution, d_currentSolution,n_students);
                copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[1]>>>(d_previousAluxcol, d_aluxcol,n_colegios);
                copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[2]>>>(d_previousAluVulxCol, d_aluVulxCol,n_colegios);
                copyVars<<<1,3,0,streams[3]>>>(d_previousVars, d_currentVars);
                copyCost<<<1,1,0,streams[4]>>>(d_costPreviousSolution,d_costCurrentSolution);
                //for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamSynchronize(streams[i]); }
                errSync = cudaGetLastError();
                errAsync = cudaDeviceSynchronize();
                if (errSync != cudaSuccess) 
                printf("10 Sync kernel error: %s\n", cudaGetErrorString(errSync));
                if (errAsync != cudaSuccess)
                printf("10 Async kernel error: %s\n", cudaGetErrorString(errAsync));

                /*
                memcpy(previousSolution,currentSolution,sizeof(int)*n_students);
                memcpy(previousAluxCol,aluxcol,sizeof(int)*n_colegios);
                memcpy(previousAluVulxCol,aluVulxCol,sizeof(int)*n_colegios);
                memcpy(previousVars,currentVars,sizeof(double)*3);
                */
                costPreviousSolution = costCurrentSolution;

                count_rechaso = 0;
                c_accepta++;
            }
            else{
                count_rechaso++;
                
            }
        }

        if(temperatureTL7(c_cooling_temperature, c_accepta, len1, len2, n_colegios, count)){
        //if(temperatureTL8(temp, c_cooling_temperature, count_trials, len1, len2, coolingRate)){
        //if(temperatureTL9(temp, c_cooling_temperature, count_trials, len3, len4, coolingRate)){
        //if(temperatureTL11(temp, c_cooling_temperature, count_trials, len3, len4, coolingRate)){
            cooling.CS2();
            cout << cooling.getTemp() << "\n";
        }

        //reheating.TR11(temp, k_reheating, n_reheating, count_rechaso);
        //reheating.TR12(temp, k_reheating, n_reheating, count);
        //reheating.TR13(temp, k_reheating, n_reheating, c_cooling_temperature);
        //reheating.TR14(temp, k_reheating, k_reheating_init, n_reheating, count_rechaso, e_const);
        
        
        
        ///////////////////////////////////////////////////
        /// History
        ///////////////////////////////////////////////////
        /*
        vector_historyCostSolution.push_back(costCurrentSolution);
        vector_historyTemp.push_back(temp);
        vector_historymeanDist.push_back(meanDist(currentSolution,distMat));
        vector_historymeanDistNorm.push_back(meanDist(currentSolution,distMat)/max_dist);
        vector_historySegregation.push_back(S(currentSolution, alumnosSep, totalVuln));
        vector_historycostoCupo.push_back(costCupo(currentSolution,cupoArray));
        if(count_rechaso==0){
            vector_historyAcceptSolution.push_back(1);
        }
        else{
            vector_historyAcceptSolution.push_back(0);
        }
        vector_historyMove.push_back(std::tuple<int,int>(shuffle_colegios[selectThread],shuffle_student[selectBlock]));     
        */
        
        //cout << costCurrentSolution << costPreviousSolution << "| |" << temp << "| |" << count<< endl;
        errSync  = cudaGetLastError();
        errAsync = cudaDeviceSynchronize();
        if (errSync != cudaSuccess) 
        printf("6 Sync kernel error: %s\n", cudaGetErrorString(errSync));
        if (errAsync != cudaSuccess)
        printf("6 Async kernel error: %s\n", cudaGetErrorString(errAsync));
        count_trials++;
        count++;
    }
    cudaMemcpyAsync(bestSolution, d_bestSolution, n_students * sizeof(int), cudaMemcpyDeviceToHost,streams[0]);
    cudaMemcpyAsync(previousSolution, d_previousSolution, n_students * sizeof(int), cudaMemcpyDeviceToHost,streams[1]);
    ///////////////////////////////////////////////////
    /// Obtiene el tiempo de ejecución
    ///////////////////////////////////////////////////
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    for(x=0;x<n_students;x++){
        info_graficos_bestSolution << bestSolution[x] << ",";
    }
    

    for(x=0; x<vector_count.size(); x++){
        info_graficos << vector_count.at(x) << "," 
                    << vector_meanDist.at(x)/max_dist << "," // Distancia promedio recorrida por los estudiantes normalizada
                    << vector_meanDist.at(x) << "," 
                    << vector_segregation.at(x) << "," 
                    << vector_costoCupo.at(x) << "," 
                    << vector_costCurrentSolution.at(x) << "," 
                    << fixed << vector_temp.at(x) << setprecision(13) << "\n";
    }




    ///////////////////////////////////////////////////
    /// Almacenamiento de datos
    ///////////////////////////////////////////////////
    cout.precision(dbl::max_digits10);
    cout << "--------------- Resultado Final ----------------" << "\n";
    cout << "Numero de Ciclos " << count << "\n";
    cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
    cout << "Costo de la mejor solución: " << costBestSolution << "\n";
    cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
    cout << "Tiempo de ejecución de SA: " << time_taken << "\n";
    cout << "distancia: " << meanDist(bestSolution,distMat) << "\n";
    cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << "CostoCupo: " << costCupo(bestSolution,cupoArray) << "\n";

    cout << "Cal costo " << calCosto(bestSolution, distMat, ptr_alpha, alumnosSep, totalVuln, cupoArray) << endl;
    cout << "Costo de: " << costBestSolution << "\n";


    //cout << fixed << setprecision(70) << endl;
    cout << sumDist(bestSolution,distMat) << "\n";
    cout << bestVars[0] << endl;
    cout << sumS(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << bestVars[1] << endl;
    cout << sumCostCupo(bestSolution,cupoArray) << "\n";
    cout << bestVars[2] << endl;
    cout << "Tiempo de ejecución de SA get_result: " << vector_time1 << "\n";

    cout << "--------------- Finalizo con exito ----------------" << "\n";


    info << "--------------- Resultado Final ----------------" << "\n";
    info << "Numero de Ciclos " << count << "\n";
    info << "Costo de la solución previa: " << costPreviousSolution << "\n";
    info << "Costo de la mejor solución: " << costBestSolution << "\n";
    info << "Costo de la solución actual: " << costCurrentSolution << "\n";
    info << "Tiempo de ejecución de SA: " << time_taken << "\n";
    info << "distancia: " << meanDist(bestSolution,distMat) << "\n";
    info << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    info << "CostoCupo: " << costCupo(bestSolution,cupoArray) << "\n";
    info << "--------------- Finalizo con exito ----------------" << "\n";


    info_test << fixed << time_taken << setprecision(9) << "," 
            << costBestSolution << "," 
            << meanDist(bestSolution,distMat)/max_dist 
            << "," << meanDist(bestSolution,distMat) 
            << "," << S(bestSolution, alumnosSep, totalVuln) 
            << "," << costCupo(bestSolution,cupoArray) 
            << "," << count 
            << "," << fixed << temp_init << setprecision(13) 
            << "," << fixed << cooling.getTemp() << setprecision(13) 
            << "," << min_temp 
            << "," << seed
            << "," << alpha1 
            << "," << alpha2 
            << "," << alpha3 
            << "," << alpha[0]
            << "," << alpha[1]
            << "," << alpha[2]
            << "," << coolingRate 
            << "," << k_reheating_init 
            << "," << e_const
            << "," << n_reheating
            << "," << len1_init
            << "," << len2_init
            << "," << len3_init
            << "," << len4_init
            << "," << len1
            << "," << len2
            << "," << len3
            << "," << len4
            << "," << Th
            << "," << n_block 
            << "," << n_thread 
            << ","<< name_exp << "\n";

    info_graficos_bestSolution.close();
    cout << ".";
    info_graficos.close();
    cout << ".";
    info_test.close();
    info.close();
    cout << ".\n";
    cout << " Archivos Guardado" << "\n";



    for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamDestroy(streams[i]); }
    cudaFree(d_array_current_Solution);
    cudaFree(d_costCurrentSolution);
    cudaFree(d_costBestSolution);
    cudaFree(d_costPreviousSolution);
    cudaFree(d_currentVars);
    cudaFree(d_bestVars);
    cudaFree(d_previousVars);
    cudaFree(d_array_current_Solution_alu);
    cudaFree(d_array_current_Solution_col);
    cudaFree(d_shuffle_colegios);
    cudaFree(d_shuffle_students);
    cudaFree(d_aluxcol);
    cudaFree(d_previousAluxcol);
    cudaFree(d_aluVulxCol);
    cudaFree(d_previousAluVulxCol);
    cudaFree(d_currentSolution);
    cudaFree(d_bestSolution);
    cudaFree(d_previousSolution);
    cudaFree(d_alumnosSep);
    cudaFree(d_cupoArray);
    cudaFree(d_distMat);
    cudaFree(d_alpha);

    cudaEventDestroy(start_cuda);
    cudaEventDestroy(stop_cuda);

    return (costBestSolution);

}

///////////////////////////////////////////////////
///////////////////////////////////////////////////


///////////////////////////////////////////////////
/// Calcula el costo
///////////////////////////////////////////////////
double calCosto(int currentSolution[], double **distMat, const double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]){
    double var1 = meanDist(currentSolution,distMat)/max_dist;
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
double meanDist(const int currentSolution[], double  **distMat){
    double sumDist=0.0;
    for(int i=0;i<n_students;i++){
        sumDist+=round_n(distMat[i][currentSolution[i]]); // distMat[estudiante][escuela]
    }
    //cout << "meanDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << n_students << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist/n_students;
}

double sumDist(const int currentSolution[], double  **distMat){
    double sumDist=0.0;
    for(int i=0;i<n_students;i++){
        sumDist+=round_n(distMat[i][currentSolution[i]]); // distMat[estudiante][escuela]
    }
    //cout << "sumDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << n_students << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist;
}


///////////////////////////////////////////////////
/// Calcula segregación por duncan
///////////////////////////////////////////////////

double S(const int currentSolution[],const int alumnosSep[], int totalVuln){
    double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
        }
    }
    return totalSesc/2.0;
}

double sumS(const int currentSolution[],const int alumnosSep[], int totalVuln){
    double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
        }
    }
    return totalSesc;
}


///////////////////////////////////////////////////
/// Calcula el costo de tener los estudiantes en las escuelas
///////////////////////////////////////////////////

double costCupo(int currentSolution[],int cupoArray[]){
    double totalcostCupo = 0.0;
    int totalAluCol = 0;
    // double a = 0.0;
    for(int j=0;j<n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        totalcostCupo+=round_n((double)totalAluCol*fabs(((double)cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2)));
    }
    return totalcostCupo/n_colegios;
}



double sumCostCupo(int currentSolution[],int cupoArray[]){
    double totalcostCupo = 0.0;
    int totalAluCol = 0;
    for(int j=0;j<n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<n_students; i++){
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

void newSolution(int currentSolution[],const int previousSolution[]){
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> dist(0, n_students);
    random_device rd2;
    mt19937 mt2(rd2());
    uniform_int_distribution<int> dist2(0, n_colegios);
    int selectStudent=dist(mt);
    int selectSchool = dist2(mt2);
    for(int x=0; x<n_students; x++){
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
void assignSchoolToArray(int previousSolution[], int bestSolution[], int currentSolution[], Info_colegio *ptr_colegios, Info_alu *ptr_students, int cupoArray[]){
    Info_alu *ptr_aux = ptr_students;
    for(int x=0;x < n_colegios;x++){
        for(int y=0; y < n_students; y++){
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
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, double **distMat){
    Info_colegio *ptr_aux = ptr_colegios;
    for(int x=0;x < n_students ;x++){
        for(int y=0; y < n_colegios; y++){
            distMat[x][y] = sqrt( pow((ptr_students->latitude - ptr_colegios->latitude),2)+pow((ptr_students->longitude - ptr_colegios->longitude),2))/1000;
            ptr_colegios++;

        }
        ptr_colegios = ptr_aux;
        ptr_students++;
    }
}

///////////////////////////////////////////////////
/// newSolution_v2, tiene como entrada la información de los estado actual de la solución, y alcula de inmediato la
/// distancia promedio, el costocupo y segregación total.
///////////////////////////////////////////////////
double newSolution_v2(int n_students,int n_colegios,int totalVuln,int aluxcol[],int aluVulxCol[],int cupoArray[],double **distMat, int currentSolution[], const double ptr_alpha[]){
    double sumDist=0;
    // double mean=0.0;
    double totalcostCupo = 0.0;
    double totalSesc = 0.0;
    int aluVulCol, aluNoVulCol,totalAluCol;
    for(int i=0;i<n_students;i++){
        sumDist+=distMat[i][currentSolution[i]]; // distMat[estudiante][escuela]
    }
    for(int n=0; n<n_colegios; n++) {
        totalAluCol = aluxcol[n];
        aluVulCol = aluVulxCol[n];
        aluNoVulCol =totalAluCol - aluVulCol;
        // Calcula el costo cupo
        totalcostCupo+=totalAluCol*fabs((cupoArray[n]-totalAluCol)/pow(((double)cupoArray[n]/2),2));
        // Calcula el total sesc
        totalSesc+=((double)1/2)*fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln)));
    }
    totalcostCupo = totalcostCupo/n_colegios;
    double var1 = (sumDist/(double)n_students)/max_dist;
    //cout << var1 << "\n";
    double var2 = totalSesc;
    //cout << var2 << "\n";
    double var3 = totalcostCupo;
    //cout << var3 << "\n";
    return (double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
}

void shuffle(int *values, const int max_change, uniform_int_distribution<int> distri) {
    int randvalue1,randvalue2,tem_value;
    for (int i = 0; i<max_change; i++) {
        randvalue1 = distri(mt);
        randvalue2 = i;
        tem_value = values[randvalue1];
        values[randvalue1] = values[randvalue2];
        values[randvalue2] = tem_value;
    }
}

void getDataSchool(vector<Info_colegio> &colegios){
    string line_colegios;
    ifstream info_school("colegios_utm.txt"); // concatenar
    int cx = 0;
    while (getline(info_school, line_colegios)) {
        stringstream linestream(line_colegios);
        string data;
        colegios.push_back(Info_colegio());
        getline(linestream, data, ',');
        colegios[cx].rbd = stoi(data);
        getline(linestream, data, ',');
        colegios[cx].latitude = stod(data);
        getline(linestream, data, ',');
        colegios[cx].longitude = stod(data);
        getline(linestream, data, ',');
        colegios[cx].num_alu = stoi(data);
        getline(linestream, data, ',');
        colegios[cx].prioritario = stoi(data);
        cx++;
    }
    info_school.close();
}

void getDataStudents(vector<Info_alu> &students, int &totalVuln)
{
    string line_student;
    ifstream info_student("alumnos_utm.txt"); // concatenar
    int cx = 0;
    while (getline(info_student, line_student)) {
        stringstream linestream(line_student);
        string data;
        students.push_back(Info_alu());
        getline(linestream, data, ',');
        students[cx].rbd = stoi(data);
        getline(linestream, data, ',');
        students[cx].latitude = stod(data);
        getline(linestream, data, ',');
        students[cx].longitude = stod(data);
        getline(linestream, data, ',');
        students[cx].sep = stoi(data);
        if (students[cx].sep == 1) {
            totalVuln++;
        }
        cx++;

    }
    info_student.close();
}

////////////////////////////////////////////////
////// Obtiene la maxima distancia que un estudiante podria llegar a recorrer
///////////////////////////////////////////////////
double getMaxDistance(double **distMat){
    double max = 0;
    for(int i=0;i<n_students;i++){
        for(int x=0;x<n_colegios;x++){
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
void normalizedAlpha(double alpha[3])
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
void initializeArray(int *aluxcol, int *previousAluxCol, int *bestAluxCol, int *aluVulxCol, int *previousAluVulxCol, int *bestAluVulxCol, int *alumnosSep, vector<Info_alu> &students,vector<Info_colegio> &colegios)
{
    for(int x = 0; x < n_colegios; x++){
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
    for(int x=0; x < n_students; x++) {
        alumnosSep[x] = students[x].sep;
    }
}


double round_n(double x)
{
    double digits = pow(10.0, DECIMAL);
    return trunc(x * digits) / digits;
}


