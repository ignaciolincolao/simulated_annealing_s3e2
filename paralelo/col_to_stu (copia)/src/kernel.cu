#include <kernel.cuh>

__constant__ int d_cupoArray[85];
__constant__ double d_alpha[3];

__global__ void newSolution_kernel(
    double *d_array_current_Solution,
    double *d_array_current_Solution_cup,
    double *d_alpha_current_Solution_seg,
    double *d_alpha_current_Solution_dis,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col,
    const int n_students,
    const int n_colegios,
    const int space_solution,
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

    /// Inicializa variables en 0
    int aluchange,
            colchange,
            newSchool,
            aluVulCol= 0,
            aluNoVulCol= 0,
            totalAluCol= 0,
            myID = blockIdx.x * blockDim.x + threadIdx.x,
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
    aluchange = d_shuffle_students[myID%n_students]; 
    colchange = d_shuffle_colegios[blockIdx.x%n_colegios];
    currentSchool = d_currentSolution[aluchange];
    newSchool = colchange;
    sumDist= d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia
    
    sumDist-=(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual 
    totalcostCupo-=((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela nueva
    totalcostCupo-=((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2)));
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
    sumDist+=(d_distMat[aluchange * pitch / sizeof(double) + newSchool]);
    // seg de la escuela actual
    totalAluCol = myCurrentAluxCol;
    aluVulCol = myCurrentAluVulxCol;
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual
    totalcostCupo+=((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    
    // seg de la escuela antigua
    totalAluCol = myNewAluxCol;
    aluVulCol = myNewAluVulxCol;
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela antigua
    totalcostCupo+=(((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2))));


    var1 = (sumDist/n_students);
    var1= (var1/max_dist);

    var2 = (totalSesc/2.0);
    //cout << var2 << "\n";
    var3 = (totalcostCupo /n_colegios);
    d_array_current_Solution[myID]= (d_alpha[0]*var1+(d_alpha[1]*var2+(d_alpha[2]*var3+0)));
    d_array_current_Solution_col[myID] = colchange;
    d_array_current_Solution_alu[myID] = aluchange;
    d_array_current_Solution_cup[myID] = totalcostCupo;
    d_alpha_current_Solution_seg[myID] = totalSesc;
    d_alpha_current_Solution_dis[myID] = sumDist;

    
}

__global__ void reduce_block_kernel(
    double *d_array_current_Solution,
    double *d_array_current_Solution_cup,
    double *d_alpha_current_Solution_seg,
    double *d_alpha_current_Solution_dis,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col,
    int space_solution,
    int* d_alumnosSep,
    int* d_aluxcol,
    int* d_aluVulxCol,
    int* d_currentSolution,
    double *d_currentVars,
    double *d_costCurrentSolution
    ){

    int myID = blockIdx.x * blockDim.x + threadIdx.x;
    //printf("%d %d %d\n",myID,space_solution,blockDim.x);
    int salto= space_solution;
    while(salto){
        if(salto-(myID+1)>myID){
            if(d_array_current_Solution[salto-(myID+1)]<d_array_current_Solution[myID]){
                d_array_current_Solution[myID]=d_array_current_Solution[salto-(myID+1)];
                d_array_current_Solution_cup[myID]=d_array_current_Solution_cup[salto-(myID+1)];
                d_alpha_current_Solution_seg[myID]=d_alpha_current_Solution_seg[salto-(myID+1)];
                d_alpha_current_Solution_dis[myID]=d_alpha_current_Solution_dis[salto-(myID+1)];
                d_array_current_Solution_alu[myID]=d_array_current_Solution_alu[salto-(myID+1)];
                d_array_current_Solution_col[myID]=d_array_current_Solution_col[salto-(myID+1)];
            }
        }
        salto = (salto/2)+(salto&(2-1));
        if(salto==1){
            salto = 0;
        }
        __syncthreads();
    }
    if(myID == 0){
        d_costCurrentSolution[0] = d_array_current_Solution[myID];


        d_currentVars[0] = d_alpha_current_Solution_dis[myID];
        d_currentVars[1] = d_alpha_current_Solution_seg[myID];
        d_currentVars[2] = d_array_current_Solution_cup[myID];

        int aluchange = d_array_current_Solution_alu[0];
        int newSchool = d_array_current_Solution_col[0];
        int currentSchool = d_currentSolution[aluchange];

        d_aluxcol[currentSchool]-=1;
        d_aluVulxCol[currentSchool]-=d_alumnosSep[aluchange];
        
        //Asigna al estudiante a la nueva escuela
        d_currentSolution[aluchange] = newSchool;
        d_aluxcol[newSchool]+=1;
        d_aluVulxCol[newSchool]+=d_alumnosSep[aluchange];
    }
    //printf("%lf\n",d_array_current_Solution[0]);
    
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
    //printf("%d,%d, %d\n",aluchange,currentSchool,colchange);
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
    sumDist-=(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);
    //printf("%lf \n",sumDist);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual 
    //printf("%lf \n",totalSesc);
    totalcostCupo-=((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    //printf("%lf \n",totalcostCupo);
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    
    totalSesc-=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela nueva

    totalcostCupo-=((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2)));
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
    sumDist+=(d_distMat[aluchange * pitch / sizeof(double) + newSchool]);
    
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));
    // costcupo escuela actual

    totalcostCupo+=((double)totalAluCol*fabs(((double)d_cupoArray[currentSchool]-totalAluCol)/pow(((double)d_cupoArray[currentSchool]/2),2)));
    //printf("%lf \n",totalcostCupo);
    //printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    // seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool];

    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela antigua

    totalcostCupo+=(((double)totalAluCol*fabs(((double)d_cupoArray[newSchool]-totalAluCol)/pow(((double)d_cupoArray[newSchool]/2),2))));
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
    //printf("%lf\n",d_costCurrentSolution[0]);
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
