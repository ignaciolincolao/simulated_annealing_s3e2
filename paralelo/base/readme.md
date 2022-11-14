// Habilitar el registro de la CPU para Nsight 
sudo sysctl -w kernel.perf_event_paranoid=2
//Habilitar los permisos para Nsight Compute
#sudo nano /etc/modprobe.d/cuda.conf
// Dentro del archivo ingresar
options nvidia "NVreg_RestrictProfilingToAdminUsers=0"
#guardar el archivo
#actualizar initframe
sudo update-initramfs -u

// Hacer las pruebas
nsys profile --stats=true -o parelelo-report-5 paralelo 
nsys profile --stats=true --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true -o parelelo-reduce-1 paralelo
sudo nsys profile --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true --stats=true -o parelelo-reduce-3 paralelo

// Con metricas
sudo nsys profile --gpu-metrics-device=0 paralelo
sudo nsys profile --cuda-memory-usage=true --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true --stats=true --gpu-metrics-device=0 -o paralelo-metrics paralelo