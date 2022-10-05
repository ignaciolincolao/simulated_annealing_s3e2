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
        i = 0,
        x = 0,
        aluVulCol= 0,
        aluNoVulCol= 0,
        totalAluCol= 0,
        myID = threadIdx.x,
        school_alu_change,
        salto= n_thread;

double  totalcostCupo= 0.0,
        totalSesc= 0.0,
        var1,
        var2,
        var3,
        result= 0.0;
/// Inicializa arrays
aluchange = d_shuffle_students[blockIdx.x];
colchange = d_shuffle_colegios[threadIdx.x];
solutions_thread[threadIdx.x] = threadIdx.x;

/// Recopila la informacion que existe en memoria global
/// a shared memory
school_alu_change = d_currentSolution[aluchange];
for (i = threadIdx.x; i< n_colegios; i=i+n_thread){
    aluxcolblock[i] = d_aluxcol[i];
    aluVulxColblock[i] = d_aluVulxCol[i];
    if(i == school_alu_change){
        aluxcolblock[school_alu_change]-=1;
        aluVulxColblock[school_alu_change]-=d_alumnosSep[aluchange];
    }
}

/// Calcula la distancia total
for (x = 0 ; x < n_students ; x++) {
    if (x != aluchange) {
        result += d_distMat[x * pitch / sizeof(double) + d_currentSolution[x]];
    }
    else {
        result += d_distMat[x * pitch / sizeof(double) + colchange];
    }
}
__syncthreads();
/// Calcula el costo cupo y la cantidad de segregación total
for(int n=0; n<n_colegios; n++){
    totalAluCol = aluxcolblock[n];
    aluVulCol = aluVulxColblock[n];
    if(n == colchange){
        totalAluCol+=1;
        aluVulCol+=d_alumnosSep[aluchange];
    }
    aluNoVulCol =totalAluCol - aluVulCol;
    // Calcula el costo cupo
    totalcostCupo+=totalAluCol*fabs((d_cupoArray[n]-totalAluCol)/pow(((double)d_cupoArray[n]/2),2));
    // Calcula el total sesc
    totalSesc+=((double)1/2)*fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln)));
}

var1 = d_alpha[0]*((result/(double)n_students)/(double)max_dist);
var2 = d_alpha[1]*totalSesc;
var3 = d_alpha[2]*(totalcostCupo/n_colegios);
solutions[myID] = var1+var2+var3;

__syncthreads();
while(salto){
    if(salto-(myID+1)>myID){
        if(school_alu_change!=colchange && d_shuffle_colegios[salto-(myID+1)]!=school_alu_change){
            if(solutions[myID]>solutions[salto-(myID+1)]){
                solutions[myID]=solutions[salto-(myID+1)];
                solutions_thread[myID]=solutions_thread[salto-(myID+1)];
            }
        }
        else{
            if(school_alu_change==colchange){
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

