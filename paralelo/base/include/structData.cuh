#ifndef STRUCT_CUH
#define STRUCT_CUH


#ifdef __CUDACC__
  #define HD __host__ __device__
#else
  #define HD
#endif

struct DataResult {
    double costSolution;
    int col;
    int stu;

    HD DataResult() : costSolution(0.0), col(-1), stu(-1) {}
    HD DataResult(double val, int c, int s) : costSolution(val), col(c), stu(s) {}

    // Necesario para thrust::sort en device y usable también en host
    HD bool operator<(const DataResult& other) const {
        return costSolution < other.costSolution;
    }
};


struct GPU_move {
    double currentVars[4];   
    int aluxcol[6]; //1: colegio origen 2: cantidad nueva colegio origen 3: cantidad nueva vulnerables colegio origen
                    //4: colegio destino 5: cantidad nueva colegio destino 6: cantidad nueva vulnerables colegio origen

};
#undef HD

#endif