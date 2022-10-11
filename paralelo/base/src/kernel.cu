#include <kernel.cuh>

__constant__ int d_cupoArray[85];
__constant__ double d_alpha[3];

__global__ void newSolution_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_thread,
    const int n_students,
    const int n_colegios,
    const int n_thread,
    const double max_dist,
    const int* __restrict__ d_alumnosSep,
    const int totalVuln,
    const int* __restrict__ d_aluxcol,
    const int* __restrict__ d_aluVulxCol,
    const int* __restrict__ d_currentSolution,
    const double* __restrict__ d_distMat,
    const int* __restrict__ d_shuffle_students,
    const int* __restrict__ d_shuffle_colegios,
    const double* __restrict__ d_currentVars,
    size_t pitch){

    /// Shared Memory
    extern __shared__ double sharedMem[];
    int* aluxcolblock = (int*)sharedMem;
    int* aluVulxColblock = (int*)&aluxcolblock[n_colegios];
    double* solutions =(double*)&aluVulxColblock[n_colegios];
    int* solutions_thread = (int*)&solutions[n_thread];
    /// Inicializa variables en 0
    int aluchange,
            colchange,
            newSchool,
            aluVulCol= 0,
            aluNoVulCol= 0,
            totalAluCol= 0,
            myID = threadIdx.x,
            currentSchool,
            salto= n_thread,
            myCurrentAluxCol = 0,
            myCurrentAluVulxCol = 0,
            myNewAluxCol = 0,
            myNewAluVulxCol = 0;

    double  totalcostCupo= 0.0,
            totalSesc= 0.0,
            var1,
            var2,
            var3,
            sumDist = 0.0;
    /// Inicializa arrays
    aluchange = d_shuffle_students[blockIdx.x];
    colchange = d_shuffle_colegios[threadIdx.x];
    solutions_thread[threadIdx.x] = threadIdx.x;
    currentSchool = d_currentSolution[aluchange];
    newSchool = colchange;
    sumDist= d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];

    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia
    sumDist-=cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual 
    totalcostCupo-=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela nueva
    totalcostCupo-=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2)));
    ////////////////////////////////////////////////////////////////
    /////// Realiza Movimiento
    ////////////////////////////////////////////////////////////////
    //ELimina el estudiante de la escuela actual
    myCurrentAluxCol = d_aluxcol[currentSchool]-1;
    myCurrentAluVulxCol = d_aluVulxCol[currentSchool]-d_alumnosSep[aluchange];
    //Asigna al estudiante a la nueva escuela
    myNewAluxCol = d_aluxcol[newSchool]+1;
    myNewAluVulxCol = d_aluVulxCol[newSchool]+d_alumnosSep[aluchange];

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////
    sumDist+=cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + newSchool]);
    // seg de la escuela actual
    totalAluCol = myCurrentAluxCol;
    aluVulCol = myCurrentAluVulxCol;
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual
    totalcostCupo+=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    
    // seg de la escuela antigua
    totalAluCol = myNewAluxCol;
    aluVulCol = myNewAluVulxCol;
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela antigua
    totalcostCupo+=cu_round_n(((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2))));


    var1 = (sumDist/n_students);
    var1= (var1/max_dist);
    //cout << var1 << "\n";
    var2 = (totalSesc/2.0);
    //cout << var2 << "\n";
    var3 = (totalcostCupo /n_colegios);
    solutions[myID] =  (double)((d_alpha[0]*var1)+(d_alpha[1]*var2)+(d_alpha[2]*var3));

    __syncthreads();
    while(salto){
        if(salto-(myID+1)>myID){
            if(currentSchool!=colchange && d_shuffle_colegios[salto-(myID+1)]!=currentSchool){
                if(solutions[myID]>solutions[salto-(myID+1)]){
                    solutions[myID]=solutions[salto-(myID+1)];
                    solutions_thread[myID]=solutions_thread[salto-(myID+1)];
                }
            }
            else{
                if(currentSchool==colchange){
                    solutions[myID]=solutions[salto-(myID+1)];
                    solutions_thread[myID]=solutions_thread[salto-(myID+1)];
                }
            }
        }
        salto = (salto/2)+(salto&(2-1));
        if(salto==1){
            salto = 0;
        }
        __syncthreads();
    }
    if(myID==0)
    {
        d_array_current_Solution[blockIdx.x] = solutions[myID];
        d_array_current_Solution_thread[blockIdx.x] = solutions_thread[myID];

    }
}

__global__ void reduce_block_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_thread,
    int *d_array_current_Solution_block,
    const int n_block){

    extern __shared__ double sharedMem[];
    double* solutions =(double*)sharedMem;
    int* solutions_block = (int*)&solutions[n_block];
    int* solutions_thread = (int*)&solutions_block[n_block];

    int myID = threadIdx.x;
    int salto= n_block;
    solutions[myID] = d_array_current_Solution[myID];
    solutions_thread[myID] = d_array_current_Solution_thread[myID];
    solutions_block[myID]= myID;
    __syncthreads();
    while(salto){
        if(salto-(myID+1)>myID){
            if(solutions[myID]>solutions[salto-(myID+1)]){
                solutions[myID]=solutions[salto-(myID+1)];
                solutions_thread[myID]=solutions_thread[salto-(myID+1)];
                solutions_block[myID]=solutions_block[salto-(myID+1)];
            }
        }
        salto = (salto/2)+(salto&(2-1));
        if(salto==1){
            salto = 0;
        }
        __syncthreads();
    }
    if(myID==0)
    {
        d_array_current_Solution[myID] = solutions[myID];
        d_array_current_Solution_thread[myID]= solutions_thread[myID];
        d_array_current_Solution_block[myID] = solutions_block[myID];
    }
}

__global__ void calculateSolution(
    double *d_array_current_Solution,
    int *d_array_current_Solution_thread,
    int *d_array_current_Solution_block,
    const int* __restrict__ d_shuffle_students,
    const int* __restrict__ d_shuffle_colegios,
    const int n_students,
    const int n_colegios,
    const int n_thread,
    const double max_dist,
    const int* __restrict__ d_alumnosSep,
    int totalVuln,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    const double* __restrict__ d_distMat,
    size_t pitch,
    double *d_currentVars,
    double *d_costCurrentSolution){

    int aluchange,
    colchange,
    newSchool,
    aluVulCol= 0,
    aluNoVulCol= 0,
    totalAluCol= 0,
    currentSchool;

    double  totalcostCupo= 0.0,
            totalSesc= 0.0,
            var1,
            var2,
            var3,
            sumDist = 0.0;
    /// Inicializa arrays

    aluchange = d_shuffle_students[d_array_current_Solution_block[0]];
    colchange = d_shuffle_colegios[d_array_current_Solution_thread[0]];
    currentSchool = d_currentSolution[aluchange];
    newSchool = colchange;
    
    
    sumDist= d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];


    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia
    sumDist-=cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);

    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual 

    totalcostCupo-=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    
    totalSesc-=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela nueva

    totalcostCupo-=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2)));
    
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
    sumDist+=cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + newSchool]);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual

    totalcostCupo+=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));

    // seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool];

    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela antigua

    totalcostCupo+=cu_round_n(((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2))));

    d_currentVars[0] = sumDist;
    d_currentVars[1] = totalSesc;
    d_currentVars[2] = totalcostCupo;

    var1 = (sumDist/n_students);
    var1= (var1/max_dist);
    //cout << var1 << "\n";
    var2 = (totalSesc/2.0);
    //cout << var2 << "\n";
    var3 = (totalcostCupo /n_colegios);
    d_costCurrentSolution[0] =  (double)((d_alpha[0]*var1)+(d_alpha[1]*var2)+(d_alpha[2]*var3));
    d_array_current_Solution[0] = d_costCurrentSolution[0];
    if(d_array_current_Solution[0] != d_costCurrentSolution[0]){
        printf("ERRORRRRRRRRRRRR no son iguales!!!!!!!!!!!!!!!!!!\n");
        printf("%lf\n",d_array_current_Solution[0]);
        printf("%lf\n",d_costCurrentSolution[0]);

    }
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

__device__ double cu_round_n(double x)
{
    double digits = pow(10.0, 16);
    return trunc(x * digits) / digits;
}
