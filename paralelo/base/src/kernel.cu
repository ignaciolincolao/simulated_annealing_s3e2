#include <kernel.cuh>

__constant__ double d_alpha[4];
__constant__ int d_n_students;
__constant__ int d_n_colegios;
__constant__ double d_max_dist;
__constant__ int d_totalVuln;

__global__ void newSolution_kernel(
    DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    const int* __restrict__ d_aluxcol,
    const int* __restrict__ d_aluVulxCol,
    const int* __restrict__ d_currentSolution,
    const double* __restrict__ d_distMat,
    const int* __restrict__ d_shuffle_students,
    const int* __restrict__ d_shuffle_colegios,
    const double* __restrict__ d_currentVars,
    size_t pitch,
    const float * __restrict__ d_penalty_matrix //matriz de penalidades

) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int     aluchange,
            newSchool,
            aluVulCol = 0,
            aluNoVulCol = 0,
            totalAluCol = 0,
            currentSchool,
            col_solution,
            alu_solution;
    double  totalcostCupo= 0.0,
            totalSesc= 0.0,
            sumDist = 0.0,
            penalty = 0.0,
            cost_solution;
    aluchange = d_shuffle_students[tid%d_n_students];
    
    newSchool = d_shuffle_colegios[0];
    currentSchool = d_currentSolution[aluchange];
    col_solution = newSchool;
    alu_solution = aluchange;
    d_array_current_Solution[tid].col = col_solution;
    d_array_current_Solution[tid].stu = alu_solution;
    
    sumDist = d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    penalty = d_currentVars[3];
    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia
    sumDist -= d_distMat[aluchange * pitch / sizeof(double) + currentSchool];
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc -= fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln)));
    // costcupo escuela actual 
    //penalty -= calcPenalty(currentSchool, choices);
    penalty -= d_penalty_matrix[aluchange * d_n_colegios + currentSchool];
    

    double p_costCupo = double(totalAluCol)/d_cupoArray[currentSchool];
    totalcostCupo -= calcCostoCupo(p_costCupo);
    //totalcostCupo -= (double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2);

    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc -= fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln)));

    // costcupo escuela nueva
    p_costCupo = double(totalAluCol)/d_cupoArray[newSchool];
    totalcostCupo -= calcCostoCupo(p_costCupo);
    //totalcostCupo -= (double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2);

    //penalty += calcPenalty(newSchool, choices);
    penalty += d_penalty_matrix[aluchange * d_n_colegios + newSchool];

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////
    sumDist += d_distMat[aluchange * pitch / sizeof(double) + newSchool];
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool]-1;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluVulCol -= d_alumnosSep[aluchange];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln)));

    // costcupo escuela actual
    p_costCupo = double(totalAluCol)/d_cupoArray[currentSchool];
    totalcostCupo += calcCostoCupo(p_costCupo);
    //totalcostCupo += (double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2);
    
    // seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool] + 1;
    aluVulCol = d_aluVulxCol[newSchool];
    aluVulCol += d_alumnosSep[aluchange];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln)));

    // costcupo escuela antigua
    p_costCupo = double(totalAluCol)/d_cupoArray[newSchool];
    totalcostCupo += calcCostoCupo(p_costCupo);
    //totalcostCupo += ((double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2));

    cost_solution = d_alpha[0] * (sumDist / (d_n_students * d_max_dist));
    cost_solution += d_alpha[1] * (totalSesc * 0.5);
    cost_solution += d_alpha[2] * (totalcostCupo / d_n_colegios);
    cost_solution += d_alpha[3] * (penalty/d_n_students);
    //printf("alpha 1 valor: %f\n", d_alpha[0]);
    //printf("alpha 2 valor: %f\n", d_alpha[1]);
    //printf("alpha 3 valor: %f\n", d_alpha[2]);
    //printf("alpha 4 valor: %f\n", d_alpha[3]);
    d_array_current_Solution[tid].costSolution  = (newSchool != currentSchool) * cost_solution + (double)(0xffffffffffffffff) * (newSchool == currentSchool);
}



__global__ void reduce_kernel(DataResult *d_array_current_Solution, int N){
    #define FULL_MASK 0xFFFFFFFF
    extern __shared__ DataResult sharedMem[];
    DataResult* solutions = (DataResult*)sharedMem;
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    DataResult val;
    val.costSolution=(double)(0xffffffffffffffff);
    val.col = -1;
    val.stu = -1;

    if (idx < N) {
        val.costSolution = d_array_current_Solution[idx].costSolution;
        val.col = d_array_current_Solution[idx].col;
        val.stu = d_array_current_Solution[idx].stu;
    }
    int warpID = threadIdx.x>>5 ;
    int threadWarp = threadIdx.x & 31;
    __syncthreads();

    // Reducción a nivel de warp, cada warp encontrara al mejor y lo dejara en la memoria compartida
    for (int salto=16; salto>0; salto>>=1){ // salto>>=1 es igual a salto/2 
        double neighbour_solution = __shfl_down_sync(FULL_MASK,val.costSolution,salto);
        int col = __shfl_down_sync(FULL_MASK,val.col,salto);
        int stu = __shfl_down_sync(FULL_MASK,val.stu,salto);
        if(neighbour_solution < val.costSolution){
            val.costSolution = neighbour_solution;
            val.col = col;
            val.stu = stu;
        }
    }

    __syncthreads();

    if (threadWarp==0){
        solutions[warpID].costSolution = val.costSolution;
        solutions[warpID].col = val.col;
        solutions[warpID].stu = val.stu;
    }
    __syncthreads();
    // Reducción entre los mejores de los warps
    DataResult val2;
    val2.costSolution=(double)(0xffffffffffffffff);
    val2.col = -1;
    val2.stu = -1;
    if(warpID == 0){
        val = (threadIdx.x < blockDim.x/32)?solutions[threadWarp]:val2;
        for(int salto=16; salto>0; salto>>=1){
            double neighbour_solution = __shfl_down_sync(FULL_MASK,val.costSolution,salto);
            int a1 = __shfl_down_sync(FULL_MASK,val.col,salto);
            int a2 = __shfl_down_sync(FULL_MASK,val.stu,salto);
            if(neighbour_solution < val.costSolution){
                val.costSolution = neighbour_solution;
                val.col = a1;
                val.stu = a2;
            }
        }
        __syncthreads();
        if(threadWarp==0){
            d_array_current_Solution[blockIdx.x].costSolution = val.costSolution;
            d_array_current_Solution[blockIdx.x].col = val.col;
            d_array_current_Solution[blockIdx.x].stu = val.stu;
        }
    }
}



__global__ void calculateSolution(
    DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    double *d_currentVars,
    double *d_costCurrentSolution,
    int id_select,
    int * d_prevMove, //para pruebas unitarias
    const float * __restrict__ d_penalty_matrix //matriz de penalidades
){

    int aluchange,
    colchange,
    newSchool,
    aluVulCol= 0,
    aluNoVulCol= 0,
    totalAluCol= 0,
    currentSchool;

    double  totalcostCupo= 0.0,
            totalSesc= 0.0,
            penalty = 0.0,
            sumDist = 0.0,
            var1,
            var2,
            var3,
            var4;
    /// Inicializa arrays

    aluchange = d_array_current_Solution[id_select].stu;
    colchange = d_array_current_Solution[id_select].col;
    currentSchool = d_currentSolution[aluchange];

    //movimiento anterior para pruebas unitarias
    d_prevMove[0] = aluchange;
    d_prevMove[1] = d_currentSolution[aluchange];

    newSchool = colchange;

    sumDist= d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    penalty = d_currentVars[3];

    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////

    // Distancia
    sumDist-=d_distMat[aluchange * pitch / sizeof(double) + currentSchool];

    // Seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
  
    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    // Costocupo escuela actual 
    double p_costCupo = double(totalAluCol)/d_cupoArray[currentSchool];
    totalcostCupo -= calcCostoCupo(p_costCupo);

    //totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    

    // Seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    
    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costocupo escuela nueva
    p_costCupo = double(totalAluCol)/d_cupoArray[newSchool];
    totalcostCupo -= calcCostoCupo(p_costCupo);

    //totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2);
    
    ////////////////////////////////////////////////////////////////
    /////// Realiza Movimiento
    ////////////////////////////////////////////////////////////////

    //ELimina el estudiante de la escuela actual
    d_aluxcol[currentSchool]-=1;
    d_aluVulxCol[currentSchool]-=d_alumnosSep[aluchange];

    //Asigna al estudiante a la nueva escuela
    d_currentSolution[aluchange] = newSchool;
    d_aluxcol[newSchool]+=1;
    d_aluVulxCol[newSchool]+=d_alumnosSep[aluchange];


    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////

    //distancia de la nueva escuela
    sumDist+=d_distMat[aluchange * pitch / sizeof(double) + newSchool];
    
    //seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Original 2 ac: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, aluchange);

    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));
    
    //costocupo escuela actual
    p_costCupo = double(totalAluCol)/d_cupoArray[currentSchool];
    totalcostCupo += calcCostoCupo(p_costCupo);
    //totalcostCupo+=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);

    //seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool];
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Original 2 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);

    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costcupo escuela antigua
    p_costCupo = double(totalAluCol)/d_cupoArray[newSchool];
    totalcostCupo += calcCostoCupo(p_costCupo);
    //totalcostCupo+=((double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2));
    

    //Calculo de penalty
    //penalty de la escuela actual
    //penalty -= calcPenalty(currentSchool, choices);
    penalty -= d_penalty_matrix[aluchange * d_n_colegios + currentSchool];

    //printf("Prev, alu: %d, col: %d, penalty: %f \n", aluchange,currentSchool,  d_penalty_matrix[aluchange * d_n_colegios + currentSchool]);
    //penalty de la escuela nueva
    //penalty += calcPenalty(newSchool, choices);
    penalty += d_penalty_matrix[aluchange * d_n_colegios + newSchool];
    //printf("Sig, alu: %d, col: %d, penalty: %f \n", aluchange,newSchool,  d_penalty_matrix[aluchange * d_n_colegios + newSchool]);

    //actualizar las sumatorias guardadas
    d_currentVars[0] = sumDist;
    d_currentVars[1] = totalSesc;
    d_currentVars[2] = totalcostCupo;
    d_currentVars[3] = penalty;

    var1 = (sumDist/d_n_students);
    var1= (var1/d_max_dist);

    var2 = (totalSesc*0.5);
    var3 = (totalcostCupo /d_n_colegios);
    var4 = (penalty / d_n_students);


    //printf("Original var1: %f, var2: %f, var3: %f, var4: %f \n", var1,var2,var3,var4);


    d_costCurrentSolution[0] = (double)((d_alpha[0] * var1) + (d_alpha[1] * var2) + (d_alpha[2] * var3) + (d_alpha[3] * var4));
}


__global__ void copyMemSolution(
    int *solution,
    int *new_solution,
    int N){
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < N; i += stride){
        solution[i] = new_solution[i];
    }
}
__global__ void copyMemCol(
    int *col,
    int *new_col,
    int N){
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < N; i += stride){
        col[i] = new_col[i];
    }
}
__global__ void copyVars(
    double *var,
    double *new_var){

    var[threadIdx.x] = new_var[threadIdx.x];
}

__global__ void copyCost(
    double *costCurrentSolution,
    double *new_costCurrentSolution
    ){

        costCurrentSolution[0] = new_costCurrentSolution[0];

    }


inline __device__ double calcPenalty(double currentCollege, uint8_t choices[5]) {
    double weights[6] = {500000, 0, 100, 200, 300, 400};
    uint8_t index = 0;
    for (size_t i = 1; i < 6; i++)
        index += (currentCollege == choices[i - 1]) * i;

    return weights[index];
}

inline __device__ double calcCostoCupo(double p_costCupo) {
    double r = 0.6;
    double s = 6.0;
    int l_izq = p_costCupo <= 0.5 && p_costCupo >= 0.0; 
    int l_der = p_costCupo > 0.5 && p_costCupo <= 1.0; 

    double penalty_sobrecupo = calcCostoCupo_sobrecupo(p_costCupo);
    double costCupoEscuela =  l_izq*(pow(2.0, r)*pow(p_costCupo, r)) + l_der*(pow(2.0, s) * pow(1.0 - p_costCupo, s))+penalty_sobrecupo;
    //printf("costo: %f, porcentaje %f, l_der: %f \n",costCupoEscuela,p_costCupo,l_der*(pow(2.0, s) * pow(1.0 - p_costCupo, s)));
    return costCupoEscuela;
}

inline __device__ double calcCostoCupo_sobrecupo(double p_costCupo) {
    const double k = 10.0;  //controla la pendiente del crecimiento exponencial
    const double S = 1.0;   //factor de escala 

    //penalidad, 0 si no hay sobrecupo, 1 creciendo exponencialmente si hay sobrecupo
    double over = fmax(0.0, p_costCupo - 1.0);
    double costCupoEscuela  = S * (exp(k * over) - 1.0);

    return costCupoEscuela;
}


__global__ void calculatePreviousSolution(
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    double *d_currentVars,
    double *d_costPrevSolUnitTest,
    int id_select,
    int * d_prevMove, //para pruebas unitarias
    const float * __restrict__ d_penalty_matrix //matriz de penalidades
){

    int aluchange,
    colchange,
    newSchool,
    aluVulCol= 0,
    aluNoVulCol= 0,
    totalAluCol= 0,
    currentSchool;

    double  totalcostCupo= 0.0,
            totalSesc= 0.0,
            penalty = 0.0,
            sumDist = 0.0,
            var1,
            var2,
            var3,
            var4;
    /// Inicializa arrays

    aluchange = d_prevMove[0];
    colchange = d_prevMove[1];
    currentSchool = d_currentSolution[aluchange];
    newSchool = colchange;

    //printf("Previous-> alu: %d, old col: %d, new col: %d \n", aluchange, colchange, currentSchool);

    sumDist= d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    penalty = d_currentVars[3];

    //printf("GPU: %lf |%lf |%lf |%lf \n",sumDist,totalSesc,totalcostCupo,penalty);
    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////

    // Distancia
    sumDist-=d_distMat[aluchange * pitch / sizeof(double) + currentSchool];

    // Seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 2 ac: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, aluchange);

    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    // Costocupo escuela actual 
    double p_costCupo = double(totalAluCol)/d_cupoArray[currentSchool];
    totalcostCupo -= calcCostoCupo(p_costCupo);
    //totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    

    //seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 2 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);
    
    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costocupo escuela nueva
    p_costCupo = double(totalAluCol)/d_cupoArray[newSchool];
    totalcostCupo -= calcCostoCupo(p_costCupo);
    //totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2);
    


    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////

    // Distancia de la nueva escuela
    sumDist+=d_distMat[aluchange * pitch / sizeof(double) + newSchool];
    
    //seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool] -1;
    aluVulCol = d_aluVulxCol[currentSchool] -d_alumnosSep[aluchange];
    aluNoVulCol = totalAluCol - aluVulCol;
    //printf("Previous 1 ac: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, aluchange);

    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));
    
    //costocupo escuela actual
    p_costCupo = double(totalAluCol)/d_cupoArray[currentSchool];
    totalcostCupo += calcCostoCupo(p_costCupo);
    //totalcostCupo+=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    

    //seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool] +1 ;
    aluVulCol = d_aluVulxCol[newSchool] +d_alumnosSep[aluchange];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 1 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n\n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);


    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costcupo escuela antigua
    p_costCupo = double(totalAluCol)/d_cupoArray[newSchool];
    totalcostCupo += calcCostoCupo(p_costCupo);
    //totalcostCupo+=((double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2));


    //Calculo de penalty
    //penalty de la escuela actual
    penalty -= d_penalty_matrix[aluchange * d_n_colegios + currentSchool];

    //penalty de la escuela nueva
    penalty += d_penalty_matrix[aluchange * d_n_colegios + newSchool]; 



    var1 = (sumDist/d_n_students);
    var1 = (var1/d_max_dist);

    var2 = (totalSesc*0.5);
    var3 = (totalcostCupo /d_n_colegios);
    var4 = (penalty / d_n_students);

    //printf("Previous var1: %d, var2: %d, var3: %d, var4: %d \n", var1,var2,var3,var4);
    
    d_costPrevSolUnitTest[0] = (double)((d_alpha[0] * var1) + (d_alpha[1] * var2) + (d_alpha[2] * var3) + (d_alpha[3] * var4));

}








/**
 * KERNEL DE PENALIZACIONES POR PREFERENCIAS CON DISCONTINUIDAD PARAMETRIZABLE
 *
 * Calcula la matriz completa P[estudiante][colegio] de penalizaciones basada en:
 * - Teoría de Satisficing: diferencia cualitativa entre opciones consideradas vs no consideradas
 * - Función exponencial parametrizable: [0.0, max_pref_penalty] para preferencias, 1.0 para no-preferencias
 * - Optimización sin bifurcaciones: máximo paralelismo GPU evitando divergencia de warps
 * - Manejo flexible de preferencias variables: cada estudiante puede tener 1 a num_schools preferencias
 *
 * PARÁMETROS:
 * @param preferences_matrix: Matriz de preferencias [num_students × max_preferences_per_student]
 *                           Cada fila es el vector de preferencias de un estudiante (índices de colegios)
 * @param num_preferences: Vector con número de preferencias válidas por estudiante (rango: [1, num_schools])
 * @param penalty_matrix: Matriz de salida [num_students × num_schools]
 * @param num_students, num_schools: Dimensiones del problema
 * @param max_preferences_per_student: Tamaño máximo del vector de preferencias (típicamente = num_schools)
 * @param alpha: Parámetro de decay exponencial (controla curvatura, recomendado: 1.0f)
 * @param max_pref_penalty: Penalización máxima para preferencias declaradas (recomendado: 0.3-0.8)
 *
 * FUNCIÓN MATEMÁTICA GENERALIZADA:
 * penalty = max_pref_penalty × (1 - exp(-α(r-1)))  si colegio está en preferencias (rank r)
 * penalty = 1.0                                    si colegio NO está en preferencias
 *
 * INTERPRETACIÓN DE PARÁMETROS:
 * - alpha bajo (0.5): preferencias más uniformes, diferencia sutil entre 1ra y última opción
 * - alpha alto (2.0): preferencias muy jerarquizadas, 1ra opción mucho mejor que 2da
 * - max_pref_penalty bajo (0.3): gran discontinuidad entre preferencias vs no-preferencias  
 * - max_pref_penalty alto (0.8): discontinuidad más sutil, "zona gris" más amplia
 */
__global__ void compute_preference_penalty_matrix(
    int* preferences_matrix,
    int* num_preferences,
    float* d_penalty_matrix,
    int num_students,
    int num_schools,
    int max_preferences_per_student,
    float alpha,
    float max_pref_penalty
) {
    // Thread mapping: cada thread procesa un par (estudiante, colegio)
    int student_id = blockIdx.y * blockDim.y + threadIdx.y;
    int school_id = blockIdx.x * blockDim.x + threadIdx.x;
   
    // Bounds check: única bifurcación necesaria para memory safety
    if (student_id >= num_students || school_id >= num_schools) return;
   
    // Obtener número de preferencias para este estudiante específico
    // CLAVE: cada estudiante puede tener diferente cantidad (1 a num_schools)
    int student_prefs = num_preferences[student_id];
   
    // BÚSQUEDA SIN BIFURCACIONES: encontrar rank del colegio en preferencias de este estudiante
    // Inicializar con valor que indica "no encontrado"
    int found_rank = max_preferences_per_student + 1;
   
    // Loop desenrollado para evitar divergencia de warps
    // IMPORTANTE: iteramos hasta max_preferences_per_student, no student_prefs
    // Esto asegura que todos los threads ejecuten el mismo número de iteraciones
    #pragma unroll 8  // Optimiza casos comunes (≤8 preferencias)
    for (int p = 0; p < max_preferences_per_student; p++) {
        // Leer preferencia en posición p para este estudiante
        //preferences_matrix da un numero entre 1 - 63, pero necesitamos que este entre 0 -62
        int pref_school = preferences_matrix[student_id * max_preferences_per_student + p]-1;
       
        // Condiciones como factores multiplicativos (0.0 o 1.0):
        // is_valid: esta posición p contiene una preferencia válida para este estudiante
        // is_match: el colegio en esta posición coincide con el school_id que estamos evaluando
        float is_valid = (float)(p < student_prefs);  // Solo las primeras student_prefs posiciones son válidas
        float is_match = (float)(pref_school == school_id);
        float found_here = is_valid * is_match;  // Multiplicación actúa como AND lógico
       
        // Update condicional sin if/else: actualizar found_rank solo si encontramos match
        // Técnica elegante: usar aritmética para hacer update selectivo sin divergencia
        int new_rank = p + 1;  // Ranks empiezan en 1 (primera preferencia = rank 1)
        found_rank = (int)(found_here * (float)new_rank + (1.0f - found_here) * (float)found_rank);
       
        // Nota: continuamos el loop completo aunque hayamos encontrado match
        // Esto mantiene sincronización perfecta entre todos los threads del warp
    }
   
    // CÁLCULO DE PENALIZACIÓN CON PARÁMETROS FLEXIBLES
    // Convertir "encontrado vs no encontrado" en factor multiplicativo
    float was_found = (float)(found_rank <= max_preferences_per_student);
   
    // Penalización exponencial parametrizable para preferencias declaradas
    // NUEVA FÓRMULA: penalty = max_pref_penalty × (1 - exp(-α(r-1)))
    // Donde r es el rank encontrado (1 = primera preferencia)
    //float pref_penalty = max_pref_penalty * (1.0f - expf(-alpha * powf(fmaxf(0.0f, (float)(found_rank - 1)), 1)));
    float pref_penalty = max_pref_penalty * (1.0f - expf(-alpha * fmaxf(0.0f, (float)(found_rank - 1))));
   
    // Penalización para colegios no considerados: salto discontinuo a 1.0
    // Esto preserva la discontinuidad teórica entre "considerado" vs "no considerado"
    float no_pref_penalty = 1.0f;
   
    // Combinar ambos casos usando aritmética pura (sin bifurcaciones)
    // Esta técnica asegura que todos los threads ejecuten las mismas operaciones
    float final_penalty = was_found * pref_penalty + (1.0f - was_found) * no_pref_penalty;
   
    // Guardar resultado en matriz global
    // Layout: penalty_matrix[estudiante * num_schools + colegio]
    d_penalty_matrix[student_id * num_schools + school_id] = final_penalty;
    
    /*
    if (student_id == 1){
        printf("colegio: %d (%d) (%d), penalty: %f \n", school_id,(student_id * num_schools + school_id), found_rank, final_penalty);
    }
    */
}