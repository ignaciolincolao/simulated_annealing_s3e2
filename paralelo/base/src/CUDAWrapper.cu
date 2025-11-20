#include <CUDAWrapper.cuh>
#include <kernel.cuh>

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

CUDAWrapper::CUDAWrapper(CUDAParams& cuParams,SimulatedParams& saParams, mt19937& mt) 
    :  cuParams(cuParams), saParams(saParams), mt(mt)
{
    cudaDeviceProp deviceProp;
    cudaGetDevice(&deviceId);
    cudaDeviceGetAttribute(&numberOfSMs, cudaDevAttrMultiProcessorCount, deviceId); // Calcula el numero de SMstream 
    cudaGetDeviceProperties(&deviceProp, 0);
    threadsPerBlock = 256;
    numberOfBlocks = 32 * numberOfSMs;
    nWarp = deviceProp.warpSize;
    streams = new cudaStream_t[NUM_STREAMS]; 
    for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamCreate(&streams[i]); }
    cudaEventCreate(&start_cuda);
    cudaEventCreate(&stop_cuda);
}
CUDAWrapper::~CUDAWrapper(){
    
    for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamDestroy(streams[i]); }
    cudaFree(d_array_current_Solution);
    cudaFree(d_costCurrentSolution);
    cudaFree(d_costBestSolution);
    cudaFree(d_costPreviousSolution);
    cudaFree(d_currentVars);
    cudaFree(d_bestVars);
    cudaFree(d_previousVars);
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

    cudaFree(d_prevMove); //pruebas unitarias
    cudaFree(d_costPrevSolUnitTest);

    //variables de update de penalty
    cudaFree(d_penalty_matrix);
    cudaFree(d_matrix_solution);
    
    cudaEventDestroy(start_cuda);
    cudaEventDestroy(stop_cuda);
}

void CUDAWrapper::memInit(
    int*& previousSolution,
    int*& bestSolution,
    int*& currentSolution,
    int*& cupoArray,
    int*& alumnosSep,
    int& totalVuln,
    int*& aluxcol,
    int*& aluVulxCol,
    double*& matrestest,
    double*& alpha,
    double*& currentVars
    
){
    cudaMalloc(&d_array_current_Solution, (cuParams.n_block*cuParams.n_thread) * sizeof(DataResult));
    cudaMalloc((void **) &d_costCurrentSolution, 1 * sizeof(double));
    cudaMalloc((void **) &d_costBestSolution, 1 * sizeof(double));
    cudaMalloc((void **) &d_costPreviousSolution, 1 * sizeof(double));
    cudaMalloc((void **) &d_currentVars, 4 * sizeof(double));
    cudaMalloc((void **) &d_bestVars, 4 * sizeof(double));
    cudaMalloc((void **) &d_previousVars, 4 * sizeof(double));
    cudaMalloc((void **) &d_shuffle_colegios, saParams.max_changes_school  * sizeof(int));
    cudaMalloc((void **) &d_shuffle_students, saParams.max_changes_students * sizeof(int));
    cudaMalloc((void **) &d_aluxcol,saParams.n_colegios * sizeof(int));
    cudaMalloc((void **) &d_previousAluxcol,saParams.n_colegios * sizeof(int));
    cudaMalloc((void **) &d_aluVulxCol,saParams.n_colegios * sizeof(int));
    cudaMalloc((void **) &d_previousAluVulxCol,saParams.n_colegios * sizeof(int));
    cudaMalloc((void **) &d_currentSolution, saParams.n_students * sizeof(int));  // Solución actual
    cudaMalloc((void **) &d_bestSolution, saParams.n_students * sizeof(int));
    cudaMalloc((void **) &d_previousSolution, saParams.n_students * sizeof(int));
    cudaMalloc((void **) &d_alumnosSep, saParams.n_students * sizeof(int)); // arreglo que contiene la id de cada usuario vulnerable
    cudaMalloc((void **) &d_cupoArray, saParams.n_colegios * sizeof(int));

    cudaMalloc((void **) &d_prevMove, 2 * sizeof(int)); //guardar movimiento anterior para pruebas unitarias
    cudaMalloc((void **) &d_costPrevSolUnitTest, 1 * sizeof(double)); //guardar el costo del movimiento anterior realizado con el nuevo kernel

    cudaMalloc((void **) &d_matrix_solution, cuParams.n_block*cuParams.n_thread*sizeof(GPU_move));
    ///////////////////////////////////////////////////
    /// Genera arreglos que contendran valores del 0 hasta saParams.n_students y saParams.n_colegios
    ///////////////////////////////////////////////////


    for (int i = 0; i < saParams.n_students; i++) {
        saParams.shuffle_student[i] = i;
    }
    for (int i=0; i < saParams.n_colegios; i++){
        saParams.shuffle_colegios[i]=i;
    }
    cudaMallocPitch(&d_distMat,
                    &pitch,
                    saParams.n_colegios * sizeof(double),
                    saParams.n_students); // Reserva memoria para la matriz de distancia


    gpuErrchk( cudaMemcpyToSymbolAsync( d_alpha, alpha,  4 * sizeof(double),0,cudaMemcpyHostToDevice,streams[2]));
    

    gpuErrchk( cudaMemcpyToSymbolAsync( d_n_students, &saParams.n_students, sizeof(int),0,cudaMemcpyHostToDevice,streams[3]));
    gpuErrchk( cudaMemcpyToSymbolAsync( d_n_colegios, &saParams.n_colegios, sizeof(int),0,cudaMemcpyHostToDevice,streams[4]));
    gpuErrchk( cudaMemcpyToSymbolAsync( d_max_dist, &saParams.max_dist, sizeof(double),0,cudaMemcpyHostToDevice,streams[5]));
    gpuErrchk( cudaMemcpyToSymbolAsync( d_totalVuln, &totalVuln, sizeof(int),0,cudaMemcpyHostToDevice,streams[6]));

    size_t h_pitchBytes = saParams.n_colegios * sizeof(double);
    cudaMemcpy2DAsync(d_distMat,
                 pitch,
                 matrestest,
                 h_pitchBytes,
                 saParams.n_colegios * sizeof(double),
                 saParams.n_students,
                 cudaMemcpyHostToDevice,
                 streams[3]);
    
    cudaMemcpyAsync(d_currentSolution, currentSolution, saParams.n_students * sizeof(int), cudaMemcpyHostToDevice,streams[2]);
    cudaMemcpyAsync(d_previousSolution, currentSolution, saParams.n_students * sizeof(int), cudaMemcpyHostToDevice,streams[3]);
    cudaMemcpyAsync(d_bestSolution, currentSolution, saParams.n_students * sizeof(int), cudaMemcpyHostToDevice,streams[4]);
    cudaMemcpyAsync(d_aluxcol, aluxcol, saParams.n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[5]);
    cudaMemcpyAsync(d_previousAluxcol, aluxcol, saParams.n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[6]);
    cudaMemcpyAsync(d_aluVulxCol, aluVulxCol, saParams.n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[7]);
    cudaMemcpyAsync(d_previousAluVulxCol, aluVulxCol,saParams.n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[8]);
    cudaMemcpyAsync(d_currentVars, currentVars, 4 * sizeof(double), cudaMemcpyHostToDevice,streams[9]);
    cudaMemcpyAsync(d_previousVars, currentVars, 4 * sizeof(double), cudaMemcpyHostToDevice,streams[0]);
    cudaMemcpyAsync(d_bestVars, currentVars, 4 * sizeof(double), cudaMemcpyHostToDevice,streams[1]);
    cudaMemcpyAsync(d_alumnosSep, alumnosSep, saParams.n_students * sizeof(int), cudaMemcpyHostToDevice,streams[2]);
    cudaMemcpyAsync(d_cupoArray, cupoArray, saParams.n_colegios * sizeof(int), cudaMemcpyHostToDevice,streams[3]);

    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    
    //error aca
    if (errSync != cudaSuccess) 
        printf("0 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("0 Async kernel error: %s\n", cudaGetErrorString(errAsync));

#ifdef ENABLE_GPU_RECORD_TIME
    time_GPU.open("../../save/tiempo_GPU.txt", std::ios::app);
#endif

}

void CUDAWrapper::memCopyPrevToCurrent(){
    copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[0]>>>(d_currentSolution, d_previousSolution,saParams.n_students);
    copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[1]>>>(d_aluxcol, d_previousAluxcol,saParams.n_colegios);
    copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[2]>>>(d_aluVulxCol, d_previousAluVulxCol,saParams.n_colegios);
    copyVars<<<1,4,0,streams[3]>>>(d_currentVars, d_previousVars);

    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
    printf("1 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
    printf("1 Async kernel error: %s\n", cudaGetErrorString(errAsync));
}

void CUDAWrapper::uploadCurrentMemorySolution(){
    cudaMemcpyAsync(d_shuffle_students, saParams.shuffle_student, saParams.max_changes_students* sizeof(int), cudaMemcpyHostToDevice,streams[0]);
    cudaMemcpyAsync(d_shuffle_colegios, saParams.shuffle_colegios, saParams.max_changes_school * sizeof(int), cudaMemcpyHostToDevice,streams[1]);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
    printf("2 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
    printf("2 Async kernel error: %s\n", cudaGetErrorString(errAsync));
};

void CUDAWrapper::AcceptanceBestSolution(){
    copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[0]>>>(d_bestSolution, d_currentSolution,saParams.n_students);
    copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[1]>>>(d_previousSolution, d_currentSolution,saParams.n_students);
    copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[2]>>>(d_previousAluxcol, d_aluxcol,saParams.n_colegios);
    copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[3]>>>(d_previousAluVulxCol, d_aluVulxCol,saParams.n_colegios);
    copyVars<<<1,4,0,streams[4]>>>(d_previousVars, d_currentVars);
    copyVars<<<1,4,0,streams[5]>>>(d_bestVars, d_currentVars);
    copyCost<<<1,1,0,streams[6]>>>(d_costBestSolution,d_costCurrentSolution);
    copyCost<<<1,1,0,streams[7]>>>(d_costPreviousSolution,d_costCurrentSolution);
    //for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamSynchronize(streams[i]); }
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
    printf("9 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
    printf("9 Async kernel error: %s\n", cudaGetErrorString(errAsync));
}

void CUDAWrapper::AcceptanceSolution(){
    copyMemSolution<<<numberOfBlocks,threadsPerBlock,0,streams[0]>>>(d_previousSolution, d_currentSolution,saParams.n_students);
    copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[1]>>>(d_previousAluxcol, d_aluxcol,saParams.n_colegios);
    copyMemCol<<<numberOfBlocks,threadsPerBlock,0,streams[2]>>>(d_previousAluVulxCol, d_aluVulxCol,saParams.n_colegios);
    copyVars<<<1,4,0,streams[3]>>>(d_previousVars, d_currentVars);
    copyCost<<<1,1,0,streams[4]>>>(d_costPreviousSolution,d_costCurrentSolution);
    //for (int i = 0; i < NUM_STREAMS; ++i) { cudaStreamSynchronize(streams[i]); }
    errSync = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
        printf("10 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("10 Async kernel error: %s\n", cudaGetErrorString(errAsync));

}

void CUDAWrapper::newSolution(){

#ifdef ENABLE_GPU_RECORD_TIME
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
#endif
    newSolution_kernel<<<cuParams.n_block, cuParams.n_thread>>>(
                                d_array_current_Solution,
                                d_cupoArray,
                                d_alumnosSep,
                                d_aluxcol,
                                d_aluVulxCol,
                                d_currentSolution,
                                d_distMat,
                                d_shuffle_students,
                                d_shuffle_colegios,
                                d_currentVars,
                                pitch,
                                d_penalty_matrix,
                                d_matrix_solution
                            );
#ifdef ENABLE_GPU_RECORD_TIME
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);

    time_GPU << std::setprecision(10)<<int(ms*1e6) <<",";
    time_GPU.flush();
#endif
    CUDAWrapper::synchronizeBucle();
#ifdef ENABLE_GPU_RECORD_TIME
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
#endif

}

void CUDAWrapper::find_minimum(){
#ifdef ENABLE_GPU_RECORD_TIME
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
#endif
    reduce_kernel<<<cuParams.n_block,cuParams.n_thread,sizeof(DataResult) *(cuParams.n_thread/32)>>>(d_array_current_Solution,cuParams.n_block*cuParams.n_thread);
    
    reduce_kernel_update<<<1,((cuParams.n_block+32-1)/32)*32,sizeof(DataResult) *((((cuParams.n_block+32-1)/32)*32)/32)>>>(
        d_array_current_Solution,
        cuParams.n_block,
        d_matrix_solution,
        d_currentVars,
        d_aluxcol,
        d_aluVulxCol,
        d_currentSolution,
        d_costCurrentSolution
    );


#ifdef ENABLE_GPU_RECORD_TIME
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    time_GPU << std::setprecision(10)<<int((ms)*1e6) <<",";
    time_GPU.flush();
#endif

    CUDAWrapper::synchronizeBucle();

#ifdef ENABLE_GPU_RECORD_TIME
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
#endif
}

/*
void CUDAWrapper::newSolutionRandomSelection(uniform_int_distribution<int> dist,
    uniform_int_distribution<int> dist2)
{

    cuParams.selectThread = dist(mt);
    cuParams.selectBlock = dist2(mt);
    cudaMemcpy(&d_array_current_Solution_alu[0],&cuParams.selectThread, sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_array_current_Solution_col,&cuParams.selectBlock, sizeof(int),cudaMemcpyHostToDevice);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
    printf("6 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
    printf("6 Async kernel error: %s\n", cudaGetErrorString(errAsync));


}
*/


void CUDAWrapper::newSolutionUpdate(double& costCurrentSolution, int id_select)
{   
#ifdef ENABLE_GPU_RECORD_TIME
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
#endif
/*
    calculateSolution<<<1,1>>>(d_array_current_Solution,
    d_cupoArray,
    d_alumnosSep,
    d_aluxcol,
    d_aluVulxCol,
    d_currentSolution,
    d_distMat,
    pitch,
    d_currentVars,
    d_costCurrentSolution,
    id_select,
    d_prevMove, //para pruebas unitarias
    d_penalty_matrix
);

*/
    getCurrentSolutionGpuToHost(costCurrentSolution);

#ifdef ENABLE_GPU_RECORD_TIME
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);

    time_GPU << std::setprecision(10)<<int(ms*1e6) <<",";
    time_GPU.flush();
#endif

    synchronizeBucle();

#ifdef ENABLE_GPU_RECORD_TIME 
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
#endif
}



void CUDAWrapper::getCurrentSolutionGpuToHost(double& costCurrentSolution)
{
    cudaMemcpy(&costCurrentSolution,&d_costCurrentSolution[0], sizeof(double),cudaMemcpyDeviceToHost);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
        printf("5 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("5 Async kernel error: %s\n", cudaGetErrorString(errAsync));
}



void CUDAWrapper::synchronizeBucle(){
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
        printf("6 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("6 Async kernel error: %s\n", cudaGetErrorString(errAsync));
}

void CUDAWrapper::copySolutionToHost(int* bestSolution, int* previousSolution){
    cudaMemcpyAsync(bestSolution, d_bestSolution, saParams.n_students * sizeof(int), cudaMemcpyDeviceToHost,streams[0]);
    cudaMemcpyAsync(previousSolution, d_previousSolution, saParams.n_students * sizeof(int), cudaMemcpyDeviceToHost,streams[1]);
    CUDAWrapper::synchronizeBucle();
}


std::tuple<int,int> CUDAWrapper::getMovementDeviceToHost(int idx){
    int alu;
    int col;
    cudaMemcpy(&alu,&d_array_current_Solution[idx].stu, sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(&col, &d_array_current_Solution[idx].col, sizeof(int),cudaMemcpyDeviceToHost);
    CUDAWrapper::synchronizeBucle();
    return std::make_tuple(alu,col);
}

void CUDAWrapper::sortSolutions(){
    int size  = cuParams.n_block * cuParams.n_thread;
    thrust::sort(thrust::device, d_array_current_Solution, d_array_current_Solution + size);
    CUDAWrapper::synchronizeBucle();
}

/*
void CUDAWrapper::UpdateSelectionDeviceToHost(int*&  currentSolution){
    cudaMemcpy(&cuParams.selectThread,&d_array_current_Solution_alu[0], sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(&cuParams.selectBlock, d_array_current_Solution_col, sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(currentSolution,d_currentSolution, saParams.n_students * sizeof(int), cudaMemcpyDeviceToHost);
}
*/
/*
void CUDAWrapper::mallocHostInit(double* currentVars,double *previousVars,double* bestVars){

}
*/

void CUDAWrapper::mallocHost(
    int*& previousSolution,
    int*&  bestSolution,
    int*&  currentSolution,
    int*&  cupoArray,
    int*&  alumnosSep,
    double*& matrestest,
    double*& currentVars,
    double*& previousVars,
    double*& bestVars){
    cudaMallocHost((void**)&previousSolution, sizeof(int)*saParams.n_students);
    cudaMallocHost((void**)&bestSolution, sizeof(int)*saParams.n_students);
    cudaMallocHost((void**)&currentSolution, sizeof(int)*saParams.n_students);
    cudaMallocHost((void**)&cupoArray, sizeof(int)*saParams.n_colegios);
    cudaMallocHost((void**)&alumnosSep, sizeof(int)*saParams.n_students);
    cudaMallocHost((void**)&matrestest,sizeof(double) * saParams.n_students * saParams.n_colegios); 
    cudaMallocHost((void**)&saParams.shuffle_student, sizeof(int)*saParams.n_students);
    cudaMallocHost((void**)&saParams.shuffle_colegios, sizeof(int)*saParams.n_colegios);
    cudaMallocHost((void**)&currentVars,4 * sizeof(double));
    cudaMallocHost((void**)&previousVars,4 * sizeof(double)); 
    cudaMallocHost((void**)&bestVars,4 * sizeof(double)); 
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
}

void CUDAWrapper::UpdateCurrentVarsHostToGPU(double*& currentVars){
    cudaMemcpy(d_currentVars, currentVars, 4 * sizeof(double),cudaMemcpyHostToDevice);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    
    if (errSync != cudaSuccess) 
        printf("0 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("0 Async kernel error: %s\n", cudaGetErrorString(errAsync));
}


void CUDAWrapper::getSolution(std::vector<DataResult>& out) {
    const int N = cuParams.n_block * cuParams.n_thread;
    out.resize(N);
    cudaMemcpy(out.data(), d_array_current_Solution,
               N * sizeof(DataResult), cudaMemcpyDeviceToHost);
}

void CUDAWrapper::previousSolution(int id_select)
{
        calculatePreviousSolution<<<1,1>>>(d_cupoArray,
        d_alumnosSep,
        d_aluxcol,
        d_aluVulxCol,
        d_currentSolution,
        d_distMat,
        pitch,
        d_currentVars,
        d_costPrevSolUnitTest,
        id_select,
        d_prevMove, //para pruebas unitarias
        d_penalty_matrix
);
        synchronizeBucle();
}

void CUDAWrapper::getPreviousSolutionUnitTest(double& costPrevSolUnitTest)
{
    cudaMemcpy(&costPrevSolUnitTest,&d_costPrevSolUnitTest[0], sizeof(double),cudaMemcpyDeviceToHost);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
        printf("5 Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
        printf("5 Async kernel error: %s\n", cudaGetErrorString(errAsync));
}


void CUDAWrapper::compute_penalty_matrix(
    int* h_preferences_matrix,        //matriz de preferencias en CPU
    int* h_num_preferences,           //número de preferencias por estudiante
    float* h_penalty_matrix          //matriz de penalizaciones en CPU (salida)
) {
    //necesitamos hacer esta reservera al inicio, ya que el metodo global se ejecuta cuando ya se hizo el primer computo CPU
    cudaMalloc((void **) &d_preferences_matrix, saParams.n_students * saParams.max_choices * sizeof(int)); //matriz de preferencias
    cudaMalloc((void **) &d_num_preferences,  saParams.n_students * sizeof(int)); //vector que indica cuantas preferencias escogio cada estudiante
    cudaMalloc((void **) &d_penalty_matrix, saParams.n_students * saParams.n_colegios * sizeof(float)); //matriz de penalidades calculadas
    
    //transferir datos a CPU
    cudaMemcpy(d_preferences_matrix, h_preferences_matrix, saParams.n_students * saParams.max_choices * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_num_preferences, h_num_preferences, saParams.n_students * sizeof(int), cudaMemcpyHostToDevice);
   
    // Configuración de grid: balance entre paralelismo y cache efficiency
    // 16×16 = 256 threads por block es óptimo para la mayoría de GPUs modernas
    dim3 block_size(16, 16);    
    dim3 grid_size(
        (saParams.n_colegios + block_size.x - 1) / block_size.x,      // Ceil division para cobertura completa
        (saParams.n_students + block_size.y - 1) / block_size.y
    );
   
    // Ejecutar kernel con parámetros flexibles
    compute_preference_penalty_matrix<<<grid_size, block_size>>>(
        d_preferences_matrix, d_num_preferences, d_penalty_matrix,
        saParams.n_students, saParams.n_colegios, saParams.max_choices,
        saParams.penalty_curve, saParams.penalty_max_pref  // Los nuevos parámetros configurables
    );
   
    cudaError_t error;
    error = cudaDeviceSynchronize();
    if (error != cudaSuccess) {
        printf("Error al calcular la matriz de penalidades: %s\n", cudaGetErrorString(error));
       
        //Cleanup en caso de error
        cudaFree(d_preferences_matrix);
        cudaFree(d_num_preferences);
        cudaFree(d_penalty_matrix);
        return;
    }
   
    //para debug guardar la matriz de penalty como un txt
    cudaMemcpy(h_penalty_matrix, d_penalty_matrix,  saParams.n_students * saParams.n_colegios * sizeof(float), cudaMemcpyDeviceToHost);

    std::ofstream fout("penalty_matrix.txt");
    if (!fout) {
        std::cerr << "Error: no se pudo abrir penalty_matrix.txt para escribir\n";
    } else {
        fout << std::fixed << std::setprecision(4);
        for (int i = 0; i < saParams.n_students; ++i) {
            for (int j = 0; j < saParams.n_colegios; ++j) {
                fout << h_penalty_matrix[i * saParams.n_colegios + j];
                if (j < saParams.n_colegios - 1) fout << ",";
            }
            fout << "\n";
        }
    }
    fout.close();

    //liberar lo que no se utilizara
    cudaFree(d_preferences_matrix);
    cudaFree(d_num_preferences);
}

