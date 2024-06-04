#!/bin/bash
#SBATCH --job-name=codigo-prueba
#SBATCH --output=salida.txt
#SBATCH --ntasks=1


cd ./build
./grid_executable
