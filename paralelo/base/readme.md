


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
cmake -DBUILD_SHARED_LIB=ON ..


cmake -DBUILD_SHARED_LIB=OFF -DENABLE_SAVE_DATA=ON -DENABLE_OPEN_RECORD_INFO=OFF -DENABLE_OPEN_RECORD_GRAPHICS=OFF -DENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION=OFF -DENABLE_OPEN_RECORD_REGISTER=ON -DENABLE_OPEN_RECORD_MOVE_SOLUTION=OFF -DCMAKE_BUILD_TYPE=Release ..


Generar txt con las asignaciones en base a MRUN y RBD
-DENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION_RBD=ON

Generar txt con los nuevos puntajes simce
-DENABLE_OPEN_RECORD_SIMCE_UPDATE=ON

PD: en windows no es necesario poner los ..

// Instalación de libreria

1-  entrar al build
2- dentro escribir:
    cmake ..
    sudo cmake -P cmake_install.cmake



Libreria necesaria es:
nlohmnann_json/3.11.3
sudo apt-get install nlohmann-json3-dev

para windows se puede utilizar vcpkg y decirle a cmake que lea las librerias de vcpkg con la flag:
-DCMAKE_TOOLCHAIN_FILE=(ruta a vcpkg)/vcpkg/scripts/buildsystems/vcpkg.cmake

nsys profile --stats=true --cuda-memory-usage=true -o parelelo-reduce-2 paralelo

cmake -DBUILD_SHARED_LIB=OFF -DENABLE_SAVE_DATA=ON -DENABLE_OPEN_RECORD_INFO=OFF -DENABLE_OPEN_RECORD_GRAPHICS=OFF -DENABLE_OPEN_RECORD_GRAPHICS_BEST_SOLUTION=OFF -DENABLE_OPEN_RECORD_REGISTER=ON -DENABLE_OPEN_ITERATION_REGISTER=OFF -DENABLE_OPEN_RECORD_MOVE_SOLUTION=OFF -DENABLE_GPU_RECORD_TIME=OFF -DCMAKE_BUILD_TYPE=Release ..