// Habilitar el registro de la CPU
sudo sysctl -w kernel.perf_event_paranoid=2
// Hacer las pruebas
nsys profile --stats=true -o parelelo-report-5 paralelo 
nsys profile --stats=true --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true -o parelelo-reduce-1 paralelo
sudo nsys profile --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true --stats=true -o parelelo-reduce-3 paralelo

// Con metricas
sudo nsys profile --gpu-metrics-device=0 paralelo
sudo nsys profile --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true --stats=true --gpu-metrics-device=0 -o paralelo-metrics paralelo



flag

cmake -DENABLE_SAVE_DATA=OFF .. 
cmake -DENABLE_SAVE_DATA=ON ..