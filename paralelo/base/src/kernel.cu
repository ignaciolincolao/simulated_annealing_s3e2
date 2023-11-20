#include <kernel.cuh>

__constant__ double d_alpha[4];
__constant__ int d_n_students;
__constant__ int d_n_colegios;
__constant__ double d_max_dist;
__constant__ int d_totalVuln;

__global__ void newSolution_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col,
    const int *__restrict__ d_cupoArray,
    const int *__restrict__ d_alumnosSep,
    const int *__restrict__ d_aluxcol,
    const int *__restrict__ d_aluVulxCol,
    const int *__restrict__ d_currentSolution,
    const double *__restrict__ d_distMat,
    const int *__restrict__ d_shuffle_students,
    const int *__restrict__ d_shuffle_colegios,
    const double *__restrict__ d_currentVars,
    size_t pitch,
    size_t penalty) {

    /// Shared Memory
    extern __shared__ double sharedMem[];
    double *solutions = (double *)sharedMem;
    int *solutions_col = (int *)&solutions[(blockDim.x >> 5) + 1]; // blockDim.x>>5  --> blockDim.x/32
    int *solutions_alu = (int *)&solutions_col[(blockDim.x >> 5) + 1];
    /// Inicializa variables en 0
    int aluchange,
        newSchool,
        aluVulCol = 0,
        aluNoVulCol = 0,
        totalAluCol = 0,
        myID = threadIdx.x,
        currentSchool;

    double totalcostCupo = 0.0,
           totalSesc = 0.0,
           sumDist = 0.0;

    /// Inicializa arrays
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    aluchange = d_shuffle_students[tid % d_n_students];
    // aluchange = d_shuffle_students[threadIdx.x];
    newSchool = d_shuffle_colegios[blockIdx.x % d_n_colegios];
    currentSchool = d_currentSolution[aluchange];
    // printf("%d|%d|%d|%d\n",newSchool,currentSchool,aluchange,tid%d_n_students);

    double cost_solution;
    int col_solution = newSchool;
    int alu_solution = aluchange;

    sumDist = d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia

    sumDist -= cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    // cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc -= cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));
    // costcupo escuela actual

    totalcostCupo -= cu_round_n((double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2));

    // seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    // cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc -= cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));

    // costcupo escuela nueva
    totalcostCupo -= cu_round_n((double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2));

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////
    sumDist += cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + newSchool]);
    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool] - 1;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluVulCol -= d_alumnosSep[aluchange];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));
    // costcupo escuela actual
    totalcostCupo += cu_round_n((double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2));

    // seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool] + 1;
    aluVulCol = d_aluVulxCol[newSchool];
    aluVulCol += d_alumnosSep[aluchange];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));

    // costcupo escuela antigua
    totalcostCupo += cu_round_n(((double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2)));

    cost_solution = d_alpha[0] * (sumDist / (d_n_students * d_max_dist));
    cost_solution += d_alpha[1] * (totalSesc * 0.5);
    cost_solution += d_alpha[2] * (totalcostCupo / d_n_colegios);
    cost_solution += d_alpha[3] * penalty;
    // printf("%.16lf %d %d\n",solutions[myID], colchange,aluchange);
    __syncthreads();

    // Encuentra el minimo
    int warpID = threadIdx.x / 32;
    int lane = threadIdx.x % 32;
#define FULL_MASK 0xffffffff
    // Encuentra minimo a nivel de warp
    for (int salto = 32 / 2; salto > 0; salto >>= 1) { // salto>>=1 es igual a salto/2
        /*
        double neightbour_solutions[2] = {cost_solution, __shfl_down_sync(FULL_MASK,cost_solution,salto)};
        int cols_solutions[2] = {col_solution, __shfl_down_sync(FULL_MASK,col_solution,salto)};
        int alus_solutions[2] = {alu_solution, __shfl_down_sync(FULL_MASK,alu_solution,salto)};
        int pos = (neightbour_solutions[1] < neightbour_solutions[0]);
        cost_solution = neightbour_solutions[pos];
        col_solution = cols_solutions[pos];
        alu_solution = alus_solutions[pos];
        */

        double neighbour_solution = __shfl_down_sync(FULL_MASK, cost_solution, salto);
        int co = __shfl_down_sync(FULL_MASK, col_solution, salto);
        int al = __shfl_down_sync(FULL_MASK, alu_solution, salto);
        if (neighbour_solution < cost_solution) {
            cost_solution = neighbour_solution;
            col_solution = co;
            alu_solution = al;
        }
    }
    if (lane == 0) {
        solutions[warpID] = cost_solution;
        solutions_col[warpID] = col_solution;
        solutions_alu[warpID] = alu_solution;
    }

    __syncthreads();
    // Encuentra el minimo a nivel de bloque
    if (warpID == 0) {
        /*
        int select= myID < blockDim.x/32;
        double cost_solution_select[2] = {9999.9,solutions[lane]};

        cost_solution = cost_solution_select[select];
        col_solution = solutions_col[lane];
        alu_solution = solutions_alu[lane];
        */

        cost_solution = (myID < blockDim.x / 32) ? solutions[lane] : 9999;
        col_solution = (myID < blockDim.x / 32) ? solutions_col[lane] : -1;
        alu_solution = (myID < blockDim.x / 32) ? solutions_alu[lane] : -1;
        // printf("laneID= %d blockIdx.x= %d | %.16lf %d %d\n",lane, blockIdx.x, cost_solution, alu_solution, col_solution);
        for (int salto = 32 / 2; salto > 0; salto >>= 1) {
            /*
            double neightbour_solutions[2] = {cost_solution, __shfl_down_sync(FULL_MASK,cost_solution,salto)};
            int cols_solutions[2] = {col_solution, __shfl_down_sync(FULL_MASK,col_solution,salto)};
            int alus_solutions[2] = {alu_solution, __shfl_down_sync(FULL_MASK,alu_solution,salto)};
            int pos = (neightbour_solutions[1] < neightbour_solutions[0]);
            cost_solution = neightbour_solutions[pos];
            col_solution = cols_solutions[pos];
            alu_solution = alus_solutions[pos];
            */
            double neighbour_solution = __shfl_down_sync(FULL_MASK, cost_solution, salto);
            int co = __shfl_down_sync(FULL_MASK, col_solution, salto);
            int al = __shfl_down_sync(FULL_MASK, alu_solution, salto);
            if (neighbour_solution < cost_solution) {
                cost_solution = neighbour_solution;
                col_solution = co;
                alu_solution = al;
            }
        }
        __syncthreads();
        if (lane == 0) {
            d_array_current_Solution[blockIdx.x] = cost_solution;
            d_array_current_Solution_alu[blockIdx.x] = alu_solution;
            d_array_current_Solution_col[blockIdx.x] = col_solution;
            // printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        }
    }
}

__global__ void reduce_block_kernel(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col) {

    extern __shared__ double sharedMem[];
    double *solutions = (double *)sharedMem;
    int *solutions_col = (int *)&solutions[blockDim.x / 32 + 1];
    int *solutions_alu = (int *)&solutions_col[blockDim.x / 32 + 1];
    int myID = threadIdx.x;
    int end = blockDim.x - 1;

    double cost_solution = d_array_current_Solution[myID];
    int col_solution = d_array_current_Solution_col[myID];
    int alu_solution = d_array_current_Solution_alu[myID];
    int warpID = threadIdx.x / 32;
    int lane = threadIdx.x % 32;
#define FULL_MASK 0xffffffff

    if (myID == 0) {
        if (d_array_current_Solution[end] < cost_solution) {
            cost_solution = d_array_current_Solution[end];
            col_solution = d_array_current_Solution_col[end];
            alu_solution = d_array_current_Solution_alu[end];
        }
    }

    // Encuentra minimo a nivel de warp
    // printf("%.16lf %d %d\n", cost_solution,col_solution,alu_solution);
    for (int salto = 32 / 2; salto > 0; salto >>= 1) { // salto>>=1 es igual a salto/2
        /*
        double neightbour_solutions[2] = {cost_solution, __shfl_down_sync(FULL_MASK,cost_solution,salto)};
        int cols_solutions[2] = {col_solution, __shfl_down_sync(FULL_MASK,col_solution,salto)};
        int alus_solutions[2] = {alu_solution, __shfl_down_sync(FULL_MASK,alu_solution,salto)};
        int pos = (neightbour_solutions[1] < neightbour_solutions[0]);
        cost_solution = neightbour_solutions[pos];
        col_solution = cols_solutions[pos];
        alu_solution = alus_solutions[pos];
        */

        double neighbour_solution = __shfl_down_sync(FULL_MASK, cost_solution, salto);
        int co = __shfl_down_sync(FULL_MASK, col_solution, salto);
        int al = __shfl_down_sync(FULL_MASK, alu_solution, salto);
        if (neighbour_solution < cost_solution) {
            cost_solution = neighbour_solution;
            col_solution = co;
            alu_solution = al;
        }
    }
    if (lane == 0) {
        solutions[warpID] = cost_solution;
        solutions_col[warpID] = col_solution;
        solutions_alu[warpID] = alu_solution;
    }

    __syncthreads();
    // Encuentra el minimo a nivel de bloque
    if (warpID == 0) {
        /*
        int select= myID < blockDim.x/32;
        double cost_solution_select[2] = {9999.9,solutions[lane]};

        cost_solution = cost_solution_select[select];
        col_solution = solutions_col[lane];
        alu_solution = solutions_alu[lane];
        */
        cost_solution = (myID < blockDim.x / 32) ? solutions[lane] : 9999;
        col_solution = (myID < blockDim.x / 32) ? solutions_col[lane] : -1;
        alu_solution = (myID < blockDim.x / 32) ? solutions_alu[lane] : -1;
        // printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        for (int salto = 32 / 2; salto > 0; salto >>= 1) {
            /*
            double neightbour_solutions[2] = {cost_solution, __shfl_down_sync(FULL_MASK,cost_solution,salto)};
            int cols_solutions[2] = {col_solution, __shfl_down_sync(FULL_MASK,col_solution,salto)};
            int alus_solutions[2] = {alu_solution, __shfl_down_sync(FULL_MASK,alu_solution,salto)};
            int pos = (neightbour_solutions[1] < neightbour_solutions[0]);
            cost_solution = neightbour_solutions[pos];
            col_solution = cols_solutions[pos];
            alu_solution = alus_solutions[pos];
            */
            double neighbour_solution = __shfl_down_sync(FULL_MASK, cost_solution, salto);
            int co = __shfl_down_sync(FULL_MASK, col_solution, salto);
            int al = __shfl_down_sync(FULL_MASK, alu_solution, salto);
            if (neighbour_solution < cost_solution) {
                cost_solution = neighbour_solution;
                col_solution = co;
                alu_solution = al;
            }
        }
        //__syncthreads();
        if (lane == 0) {
            d_array_current_Solution[blockIdx.x] = cost_solution;
            d_array_current_Solution_alu[blockIdx.x] = alu_solution;
            d_array_current_Solution_col[blockIdx.x] = col_solution;
            // printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        }
    }
}

__global__ void reduce_block_max(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col) {

    extern __shared__ double sharedMem[];
    double *solutions = (double *)sharedMem;
    int *solutions_col = (int *)&solutions[blockDim.x / 32 + 1];
    int *solutions_alu = (int *)&solutions_col[blockDim.x / 32 + 1];
    int myID = threadIdx.x;
    int end = blockDim.x - 1;

    double cost_solution = d_array_current_Solution[myID];
    int col_solution = d_array_current_Solution_col[myID];
    int alu_solution = d_array_current_Solution_alu[myID];
    int warpID = threadIdx.x / 32;
    int lane = threadIdx.x % 32;
#define FULL_MASK 0xffffffff

    if (myID == 0) {
        if (d_array_current_Solution[end] < cost_solution) {
            cost_solution = d_array_current_Solution[end];
            col_solution = d_array_current_Solution_col[end];
            alu_solution = d_array_current_Solution_alu[end];
        }
    }

    // Encuentra minimo a nivel de warp
    // printf("%.16lf %d %d\n", cost_solution,col_solution,alu_solution);
    for (int salto = 32 / 2; salto > 0; salto >>= 1) { // salto>>=1 es igual a salto/2
        double neightbour_solutions[2] = {cost_solution, __shfl_down_sync(FULL_MASK, cost_solution, salto)};
        int cols_solutions[2] = {col_solution, __shfl_down_sync(FULL_MASK, col_solution, salto)};
        int alus_solutions[2] = {alu_solution, __shfl_down_sync(FULL_MASK, alu_solution, salto)};
        int pos = (neightbour_solutions[1] < neightbour_solutions[0]);
        cost_solution = neightbour_solutions[pos];
        col_solution = cols_solutions[pos];
        alu_solution = alus_solutions[pos];
        /*
        double neighbour_solution = __shfl_down_sync(FULL_MASK,cost_solution,salto);
        int co = __shfl_down_sync(FULL_MASK,col_solution,salto);
        int al = __shfl_down_sync(FULL_MASK,alu_solution,salto);
        if(neighbour_solution < cost_solution){
            cost_solution = neighbour_solution;
            col_solution = co;
            alu_solution = al;
        }
        */
    }
    if (lane == 0) {
        solutions[warpID] = cost_solution;
        solutions_col[warpID] = col_solution;
        solutions_alu[warpID] = alu_solution;
    }

    __syncthreads();
    // Encuentra el minimo a nivel de bloque
    if (warpID == 0) {
        cost_solution = (myID < blockDim.x / 32) ? solutions[lane] : 0.000000000000;
        col_solution = (myID < blockDim.x / 32) ? solutions_col[lane] : 0;
        alu_solution = (myID < blockDim.x / 32) ? solutions_alu[lane] : 0;
        // printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        for (int salto = 32 / 2; salto > 0; salto >>= 1) {
            double neightbour_solutions[2] = {cost_solution, __shfl_down_sync(FULL_MASK, cost_solution, salto)};
            int cols_solutions[2] = {col_solution, __shfl_down_sync(FULL_MASK, col_solution, salto)};
            int alus_solutions[2] = {alu_solution, __shfl_down_sync(FULL_MASK, alu_solution, salto)};
            int pos = (neightbour_solutions[1] < neightbour_solutions[0]);
            cost_solution = neightbour_solutions[pos];
            col_solution = cols_solutions[pos];
            alu_solution = alus_solutions[pos];
            /*
            double neighbour_solution = __shfl_down_sync(FULL_MASK,cost_solution,salto);
            int co = __shfl_down_sync(FULL_MASK,col_solution,salto);
            int al = __shfl_down_sync(FULL_MASK,alu_solution,salto);
            if(neighbour_solution < cost_solution){
                cost_solution = neighbour_solution;
                col_solution = co;
                alu_solution = al;
            }
            */
        }
        //__syncthreads();
        if (lane == 0) {
            d_array_current_Solution[blockIdx.x] = cost_solution;
            d_array_current_Solution_alu[blockIdx.x] = alu_solution;
            d_array_current_Solution_col[blockIdx.x] = col_solution;
            // printf("laneID= %d %.16lf %d %d\n",lane, cost_solution, alu_solution, col_solution);
        }
    }
}

__global__ void calculateSolution(
    double *d_array_current_Solution,
    int *d_array_current_Solution_alu,
    int *d_array_current_Solution_col,
    const int *__restrict__ d_cupoArray,
    const int *__restrict__ d_alumnosSep,
    int *d_aluxcol,
    int *d_aluVulxCol,
    int *d_currentSolution,
    const double *__restrict__ d_distMat,
    size_t pitch,
    size_t penalty,
    double *d_currentVars,
    double *d_costCurrentSolution) {

    int aluchange,
        colchange,
        newSchool,
        aluVulCol = 0,
        aluNoVulCol = 0,
        totalAluCol = 0,
        currentSchool;

    double totalcostCupo = 0.0,
           totalSesc = 0.0,
           var1,
           var2,
           var3,
           sumDist = 0.0;

    size_t var4;
    /// Inicializa arrays

    aluchange = d_array_current_Solution_alu[0];
    colchange = d_array_current_Solution_col[0];
    currentSchool = d_currentSolution[aluchange];
    // printf("%d \t %.20lf | %d %d %d \n",blockIdx.x,d_array_current_Solution[0],d_array_current_Solution_alu[0],d_array_current_Solution_col[0],currentSchool);
    newSchool = colchange;

    sumDist = d_currentVars[0];
    totalSesc = d_currentVars[1];
    totalcostCupo = d_currentVars[2];
    // printf("%lf |%lf |%lf |%lf |%d |%d \n",sumDist,totalSesc,totalcostCupo,d_array_current_Solution[0],aluchange,colchange);

    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia
    sumDist -= cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + currentSchool]);
    // printf("%lf \n",sumDist);
    //  seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];

    // cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc -= cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));
    // costcupo escuela actual
    // printf("%lf \n",totalSesc);
    totalcostCupo -= cu_round_n((double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2));
    // printf("%lf \n",totalcostCupo);
    // printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    //  seg de la escuela nueva
    totalAluCol = d_aluxcol[newSchool];
    // cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol = totalAluCol - aluVulCol;

    totalSesc -= cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));

    // costcupo escuela nueva

    totalcostCupo -= cu_round_n((double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2));
    // printf("a%d \n",newSchool);
    // printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    ////////////////////////////////////////////////////////////////
    /////// Realiza Movimiento
    ////////////////////////////////////////////////////////////////
    // ELimina el estudiante de la escuela actual
    d_aluxcol[currentSchool] -= 1;
    d_aluVulxCol[currentSchool] -= d_alumnosSep[aluchange];
    // Asigna al estudiante a la nueva escuela
    d_currentSolution[aluchange] = newSchool;
    d_aluxcol[newSchool] += 1;
    d_aluVulxCol[newSchool] += d_alumnosSep[aluchange];

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////
    sumDist += cu_round_n(d_distMat[aluchange * pitch / sizeof(double) + newSchool]);

    // seg de la escuela actual
    totalAluCol = d_aluxcol[currentSchool];
    aluVulCol = d_aluVulxCol[currentSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));
    // costcupo escuela actual

    totalcostCupo += cu_round_n((double)totalAluCol * fabs((double)d_cupoArray[currentSchool] - totalAluCol) / pow(((double)d_cupoArray[currentSchool] * 0.5), 2));
    // printf("%lf \n",totalcostCupo);
    // printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    //  seg de la escuela antigua
    totalAluCol = d_aluxcol[newSchool];

    aluVulCol = d_aluVulxCol[newSchool];
    aluNoVulCol = totalAluCol - aluVulCol;
    totalSesc += cu_round_n(fabs((aluVulCol / (double)d_totalVuln) - (aluNoVulCol / (double)(d_n_students - d_totalVuln))));

    // costcupo escuela antigua

    totalcostCupo += cu_round_n(((double)totalAluCol * fabs((double)d_cupoArray[newSchool] - totalAluCol) / pow(((double)d_cupoArray[newSchool] * 0.5), 2)));
    // printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    d_currentVars[0] = sumDist;
    d_currentVars[1] = totalSesc;
    d_currentVars[2] = totalcostCupo;
    // printf("%lf %lf %lf %d %d %d\n",sumDist,totalSesc,totalcostCupo,aluchange,colchange,currentSchool);
    var1 = (sumDist / d_n_students);
    var1 = (var1 / d_max_dist);
    // cout << var1 << "\n";
    var2 = (totalSesc * 0.5);
    // cout << var2 << "\n";
    var3 = (totalcostCupo / d_n_colegios);
    // std::cout << var4 << "\n";
    var4 = penalty;

    d_costCurrentSolution[0] = (double)((d_alpha[0] * var1) + (d_alpha[1] * var2) + (d_alpha[2] * var3) + (d_alpha[3] * var4));
    d_array_current_Solution[0] = d_costCurrentSolution[0];
    if (d_array_current_Solution[0] != d_costCurrentSolution[0]) {
        printf("ERRORRRRRRRRRRRR no son iguales!!!!!!!!!!!!!!!!!!\n");
        printf("%.10lf\n", d_array_current_Solution[0]);
        printf("%.10lf\n", d_costCurrentSolution[0]);
        printf("%.10lf\n", var1);
        printf("%.10lf\n", var2);
        printf("%.10lf\n", var3);
        printf("%.10lf\n", d_alpha[0]);
        printf("%.10lf\n", d_alpha[1]);
        printf("%.10lf\n", d_alpha[2]);
    }
    // d_array_current_Solution[0] = d_costCurrentSolution[0];
}

__global__ void copyMemSolution(
    int *solution,
    int *new_solution,
    int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += stride) {
        solution[i] = new_solution[i];
    }
}
__global__ void copyMemCol(
    int *col,
    int *new_col,
    int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < N; i += stride) {
        col[i] = new_col[i];
    }
}
__global__ void copyVars(
    double *var,
    double *new_var) {

    var[threadIdx.x] = new_var[threadIdx.x];
}

__global__ void copyCost(
    double *costCurrentSolution,
    double *new_costCurrentSolution) {

    costCurrentSolution[0] = new_costCurrentSolution[0];
}

inline __device__ double cu_round_n(double x) {
    // pow(10.0, 16);
    double digits = 10000000000000000;
    return trunc(x * digits) / digits;
}

__device__ size_t calcPenalty(double *currentSolution) {
}
