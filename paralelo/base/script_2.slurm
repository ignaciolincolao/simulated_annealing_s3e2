#!/bin/bash
#SBATCH --job-name=codigo-prueba
#SBATCH --output=salida.txt
#SBATCH --ntasks=1


cd ./build
./grid_executable 1602.26 0.999 2.78021 9.89461
