nsys profile --stats=true -o parelelo-report-5 paralelo 
nsys profile --stats=true --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true -o parelelo-reduce-1 paralelo
nsys profile --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true --stats=true -o parelelo-reduce-3 paralelo