#ifndef STRUCT_CUH
#define STRUCT_CUH



struct DataResult {
    double costSolution; 
    int col;
    int stu;

    __host__ __device__ DataResult() : costSolution(0), col(-1), stu(-1) {}
   __host__ __device__ DataResult(double val, int c, int s) : costSolution(val), col(c), stu(s) {}
};



__host__ __device__ inline bool operator<(const DataResult& lhs, const DataResult& rhs) {
    return lhs.costSolution < rhs.costSolution;
}



#endif