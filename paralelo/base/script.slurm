#!/bin/bash
#SBATCH --job-name=codigo-prueba
#SBATCH --output=salida.txt
#SBATCH --ntasks=1


cd ./build
./grid_executable 100000 0.98 1 2
