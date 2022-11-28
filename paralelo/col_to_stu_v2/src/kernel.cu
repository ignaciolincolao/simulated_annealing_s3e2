#include <kernel.cuh>

__constant__ int d_cupoArray[85];
__constant__ double d_alpha[3];



__global__ void newSolution_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col,
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
    double* solutions =(double*)sharedMem;
    int* solutions_col = (int*)&solutions[blockDim.x/32+1];
    int* solutions_alu =  (int*)&solutions_col[blockDim.x/32+1];
    /// Inicializa variables en 0
    int aluchange,
            colchange,
            newSchool,
            aluVulCol= 0,
            aluNoVulCol= 0,
            totalAluCol= 0,
            myID = threadIdx.x,
            currentSchool,
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
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    aluchange = d_shuffle_students[tid%n_students]; 
    //aluchange = d_shuffle_students[threadIdx.x]; 
    colchange = d_shuffle_colegios[blockIdx.x%n_colegios];
    currentSchool = d_currentSolution[aluchange];
    //printf("%d|%d|%d|%d\n",colchange,currentSchool,aluchange,tid%n_students);

    double cost_solution = 9999.9;
    int col_solution = colchange;
    int alu_solution = aluchange;


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
    cost_solution =  (double)((d_alpha[0]*var1)+(d_alpha[1]*var2)+(d_alpha[2]*var3));
    //printf("%.16lf %d %d\n",solutions[myID], colchange,aluchange);
    __syncthreads();

    int warpID = threadIdx.x / 32;
    int lane = threadIdx.x % 32;
    #define FULL_MASK 0xffffffff
    // Encuentra minimo a nivel de warp
    for(int salto=32/2; salto>0; salto>>=1){ // salto>>=1 es igual a salto/2 
        double neighbour_solution = __shfl_down_sync(FULL_MASK,cost_solution,salto);
        int co = __shfl_down_sync(FULL_MASK,col_solution,salto);
        int al = __shfl_down_sync(FULL_MASK,alu_solution,salto);
        if(neighbour_solution < cost_solution){
            cost_solution = neighbour_solution;
            col_solution = co;
            alu_solution = al;
        }
        __syncthreads();
    }
    if(lane==0){
        solutions[warpID] = cost_solution;
        solutions_col[warpID] = col_solution;
        solutions_alu[warpID] = alu_solution;
    }
    
    __syncthreads();
    // Encuentra el minimo a nivel de bloque
    if(warpID == 0){
        cost_solution = (myID < blockDim.x/32)?solutions[lane]:9999.9;
        col_solution = (myID < blockDim.x/32)?solutions_col[lane]:0;
        alu_solution = (myID < blockDim.x/32)?solutions_alu[lane]:0;
        //printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        for(int salto=32/2; salto >0; salto>>=1){
            double neighbour_solution = __shfl_down_sync(FULL_MASK,cost_solution,salto);
            int co = __shfl_down_sync(FULL_MASK,col_solution,salto);
            int al = __shfl_down_sync(FULL_MASK,alu_solution,salto);
            if(neighbour_solution < cost_solution){
                cost_solution = neighbour_solution;
                col_solution = co;
                alu_solution = al;
            }
        }
        //__syncthreads();
        if(lane==0){
            d_array_current_Solution[blockIdx.x] = cost_solution;
            d_array_current_Solution_alu[blockIdx.x] = alu_solution;
            d_array_current_Solution_col[blockIdx.x] = col_solution;
            //printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        }
       
    }

    
}

__global__ void reduce_block_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col){

    extern __shared__ double sharedMem[];
    double* solutions =(double*)sharedMem;
    int* solutions_col = (int*)&solutions[blockDim.x/32+1];
    int* solutions_alu =  (int*)&solutions_col[blockDim.x/32+1];
    int myID = threadIdx.x;
    int end = blockDim.x-1;


    double cost_solution = d_array_current_Solution[myID];
    int col_solution = d_array_current_Solution_col[myID];
    int alu_solution = d_array_current_Solution_alu[myID];
    int warpID = threadIdx.x / 32;
    int lane = threadIdx.x % 32;
    #define FULL_MASK 0xffffffff
    
    if(myID==0){
        if(d_array_current_Solution[end] < cost_solution){
            cost_solution = d_array_current_Solution[end];
            col_solution = d_array_current_Solution_col[end];
            alu_solution = d_array_current_Solution_alu[end];
        }
    }




    // Encuentra minimo a nivel de warp
    //printf("%.16lf %d %d\n", cost_solution,col_solution,alu_solution);
    for(int salto=32/2; salto>0; salto>>=1){ // salto>>=1 es igual a salto/2 
        double neighbour_solution = __shfl_down_sync(FULL_MASK,cost_solution,salto);
        int co = __shfl_down_sync(FULL_MASK,col_solution,salto);
        int al = __shfl_down_sync(FULL_MASK,alu_solution,salto);
        if(neighbour_solution < cost_solution){
            cost_solution = neighbour_solution;
            col_solution = co;
            alu_solution = al;
        }
    }
    if(lane==0){
        solutions[warpID] = cost_solution;
        solutions_col[warpID] = col_solution;
        solutions_alu[warpID] = alu_solution;
    }
    


    __syncthreads();
    // Encuentra el minimo a nivel de bloque
    if(warpID == 0){
        cost_solution = (myID < blockDim.x/32)?solutions[lane]:9999.9;
        col_solution = (myID < blockDim.x/32)?solutions_col[lane]:0;
        alu_solution = (myID < blockDim.x/32)?solutions_alu[lane]:0;
        //printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        for(int salto=32/2; salto >0; salto>>=1){
            double neighbour_solution = __shfl_down_sync(FULL_MASK,cost_solution,salto);
            int co = __shfl_down_sync(FULL_MASK,col_solution,salto);
            int al = __shfl_down_sync(FULL_MASK,alu_solution,salto);
            if(neighbour_solution < cost_solution){
                cost_solution = neighbour_solution;
                col_solution = co;
                alu_solution = al;
            }
        }
        //__syncthreads();
        if(lane==0){
            d_array_current_Solution[blockIdx.x] = cost_solution;
            d_array_current_Solution_alu[blockIdx.x] = alu_solution;
            d_array_current_Solution_col[blockIdx.x] = col_solution;
            //printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        }
       
    }
    /*
    __syncthreads();
    while(salto){
        if(salto-(myID+1)>myID){
            if(solutions[salto-(myID+1)]<solutions[myID]){
                solutions[myID]=solutions[salto-(myID+1)];
                solutions_alu[myID]=solutions_alu[salto-(myID+1)];
                solutions_col[myID]=solutions_col[salto-(myID+1)];
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
        //printf("\n %.10lf\n ",solutions[myID]);
        d_array_current_Solution[0] = solutions[0];
        d_array_current_Solution_alu[0] = solutions_alu[0];
        d_array_current_Solution_col[0] = solutions_col[0];
        //printf("%d \t %.20lf | %d %d \n",blockIdx.x,solutions[0],solutions_alu[0],solutions_col[0]);
    }
    */
}

__global__ void calculateSolution(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col,
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

    aluchange = d_array_current_Solution_alu[0];
    colchange = d_array_current_Solution_col[0];
    currentSchool = d_currentSolution[aluchange];
    //printf("%d \t %.20lf | %d %d %d \n",blockIdx.x,d_array_current_Solution[0],d_array_current_Solution_alu[0],d_array_current_Solution_col[0],currentSchool);
    newSchool = colchange;

    
    sumDist= d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    //printf("%lf |%lf |%lf |%lf |%d |%d \n",sumDist,totalSesc,totalcostCupo,d_array_current_Solution[0],aluchange,colchange);

    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia
    sumDist-=cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);
    //printf("%lf \n",sumDist);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual 
    //printf("%lf \n",totalSesc);
    totalcostCupo-=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    //printf("%lf \n",totalcostCupo);
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    
    totalSesc-=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela nueva

    totalcostCupo-=cu_round_n((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2)));
    //printf("a%d \n",newSchool);
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
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
    //printf("%lf \n",totalcostCupo);
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    // seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool];

    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=cu_round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela antigua

    totalcostCupo+=cu_round_n(((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2))));
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    d_currentVars[0] = sumDist;
    d_currentVars[1] = totalSesc;
    d_currentVars[2] = totalcostCupo;
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    var1 = (sumDist/n_students);
    var1= (var1/max_dist);
    //cout << var1 << "\n";
    var2 = (totalSesc/2.0);
    //cout << var2 << "\n";
    var3 = (totalcostCupo /n_colegios);
    d_costCurrentSolution[0] =  (double)((d_alpha[0]*var1)+(d_alpha[1]*var2)+(d_alpha[2]*var3));
    d_array_current_Solution[0] = d_costCurrentSolution[0];
    if(d_array_current_Solution[0]!=d_costCurrentSolution[0]){
        printf("ERRORRRRRRRRRRRR no son iguales!!!!!!!!!!!!!!!!!!\n");
        printf("%.10lf\n",d_array_current_Solution[0]);
        printf("%.10lf\n",d_costCurrentSolution[0]);

    }
    //d_array_current_Solution[0] = d_costCurrentSolution[0];
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

inline __device__ double cu_round_n(double x)
{
    double digits = pow(10.0, 16);
    return trunc(x * digits) / digits;
}
