#include <kernel.cuh>

__constant__ double d_alpha[4];
__constant__ int d_n_students;
__constant__ int d_n_colegios;
__constant__ double d_max_dist;
__constant__ int d_totalVuln;
__constant__ double d_weight_n_students;

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
    const uint8_t *__restrict__ d_choices,
    size_t pitch) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int aluchange,
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
    
    //arreglo de preferencias de estudiantes (robado del oscar, todavia no se bien lo que hace)
    uint8_t choices[5] = {
        d_choices[aluchange * 5 + 0],
        d_choices[aluchange * 5 + 1],
        d_choices[aluchange * 5 + 2],
        d_choices[aluchange * 5 + 3],
        d_choices[aluchange * 5 + 4],
    };

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
    penalty -= calcPenalty(currentSchool, choices);
    
    totalcostCupo -= (double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2);

    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc -= fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln)));

    // costcupo escuela nueva
    totalcostCupo -= (double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2);

    penalty += calcPenalty(newSchool, choices);

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
    totalcostCupo += (double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2);
    
    // seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool] + 1;
    aluVulCol = d_aluVulxCol[newSchool];
    aluVulCol += d_alumnosSep[aluchange];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln)));

    // costcupo escuela antigua
    totalcostCupo += ((double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2));

    cost_solution = d_alpha[0] * (sumDist / (d_n_students * d_max_dist));
    cost_solution += d_alpha[1] * (totalSesc * 0.5);
    cost_solution += d_alpha[2] * (totalcostCupo / d_n_colegios);
    cost_solution += d_alpha[3] * (penalty / d_weight_n_students);
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
    uint8_t *d_choices,
    double *d_currentVars,
    double *d_costCurrentSolution,
    int id_select,
    int * d_prevMove //para pruebas unitarias
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

    //para pruebas unitarias
    d_prevMove[0] = aluchange;
    d_prevMove[1] = d_currentSolution[aluchange];
    //fin

    newSchool = colchange;

    uint8_t choices[5] = {
        d_choices[aluchange * 5 + 0],
        d_choices[aluchange * 5 + 1],
        d_choices[aluchange * 5 + 2],
        d_choices[aluchange * 5 + 3],
        d_choices[aluchange * 5 + 4],
    };
    
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

    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Original 1 ac: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, aluchange);

    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costocupo escuela actual 
    totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    


    //seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Original 1 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);

    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costocupo escuela nueva
    totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2);
    


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
    totalcostCupo+=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    

    //seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool];
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Original 2 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);

    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costcupo escuela antigua
    totalcostCupo+=((double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2));


    //Calculo de penalty
    //penalty de la escuela actual
    penalty -= calcPenalty(currentSchool, choices);

    //penalty de la escuela nueva
    penalty += calcPenalty(newSchool, choices);


    //actualizar las sumatorias guardadas
    d_currentVars[0] = sumDist;
    d_currentVars[1] = totalSesc;
    d_currentVars[2] = totalcostCupo;
    d_currentVars[3] = penalty;

    var1 = (sumDist/d_n_students);
    var1= (var1/d_max_dist);

    var2 = (totalSesc*0.5);
    var3 = (totalcostCupo /d_n_colegios);
    var4 = (penalty / d_weight_n_students);

    //printf("Original var1: %d, var2: %d, var3: %d, var4: %d \n", var1,var2,var3,var4);

    d_costCurrentSolution[0] = (double)((d_alpha[0] * var1) + (d_alpha[1] * var2) + (d_alpha[2] * var3) + (d_alpha[3] * var4));
    //printf("Original cost: %.16f \n", d_costCurrentSolution[0]);
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








__global__ void calculatePreviousSolution(
    //DataResult *d_array_current_Solution,
    const int* __restrict__ d_cupoArray,
    const int* __restrict__ d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    uint8_t *d_choices,
    double *d_currentVars,
    double *d_costPrevSolUnitTest,
    int id_select,
    int * d_prevMove //para pruebas unitarias
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


    uint8_t choices[5] = {
        d_choices[aluchange * 5 + 0],
        d_choices[aluchange * 5 + 1],
        d_choices[aluchange * 5 + 2],
        d_choices[aluchange * 5 + 3],
        d_choices[aluchange * 5 + 4],
    };
    
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

    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 2 ac: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, aluchange);

    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costocupo escuela actual 
    totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    


    //seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 2 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);
    
    totalSesc-=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costocupo escuela nueva
    totalcostCupo-=(double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2);
    


    ////////////////////////////////////////////////////////////////
    /////// Realiza Movimiento
    ////////////////////////////////////////////////////////////////
/*
    //ELimina el estudiante de la escuela actual
    d_aluxcol[currentSchool]-=1;
    d_aluVulxCol[currentSchool]-=d_alumnosSep[aluchange];

    //Asigna al estudiante a la nueva escuela
    d_currentSolution[aluchange] = newSchool;
    d_aluxcol[newSchool]+=1;
    d_aluVulxCol[newSchool]+=d_alumnosSep[aluchange];
*/

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////

    //distancia de la nueva escuela
    sumDist+=d_distMat[aluchange * pitch / sizeof(double) + newSchool];
    
    //seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool]-1;
    aluVulCol = d_aluVulxCol[currentSchool] -d_alumnosSep[aluchange];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 1 ac: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n",totalAluCol,aluVulCol,aluNoVulCol, aluchange);

    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));
    
    //costocupo escuela actual
    totalcostCupo+=(double)totalAluCol*fabs((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]*0.5),2);
    

    //seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool] +1 ;
    aluVulCol = d_aluVulxCol[newSchool] +d_alumnosSep[aluchange];
    aluNoVulCol =totalAluCol - aluVulCol;
    //printf("Previous 1 an: alucol: %d | aluvulcol: %d | alunovulcol: %d | aluchange: %d \n\n",totalAluCol,aluVulCol,aluNoVulCol, newSchool);


    totalSesc+=fabs((aluVulCol/(double)d_totalVuln)-(aluNoVulCol/(double)(d_n_students-d_totalVuln)));

    //costcupo escuela antigua
    totalcostCupo+=((double)totalAluCol*fabs((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]*0.5),2));


    //Calculo de penalty
    //penalty de la escuela actual
    penalty -= calcPenalty(currentSchool, choices);

    //penalty de la escuela nueva
    penalty += calcPenalty(newSchool, choices);

/*
    //actualizar las sumatorias guardadas
    d_currentVars[0] = sumDist;
    d_currentVars[1] = totalSesc;
    d_currentVars[2] = totalcostCupo;
    d_currentVars[3] = penalty;
*/

    var1 = (sumDist/d_n_students);
    var1= (var1/d_max_dist);

    var2 = (totalSesc*0.5);
    var3 = (totalcostCupo /d_n_colegios);
    var4 = (penalty / d_weight_n_students);

    //printf("Previous var1: %d, var2: %d, var3: %d, var4: %d \n", var1,var2,var3,var4);
    
    d_costPrevSolUnitTest[0] = (double)((d_alpha[0] * var1) + (d_alpha[1] * var2) + (d_alpha[2] * var3) + (d_alpha[3] * var4));
    //printf("Previous cost: %.16f \n\n", d_costPrevSolUnitTest[0]);
}