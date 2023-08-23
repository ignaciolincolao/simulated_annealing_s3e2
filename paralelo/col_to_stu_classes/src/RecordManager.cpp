#include <RecordManager.hpp>

#include <iomanip>

RecordManager::RecordManager(SimulatedParams &saParams_, RecordParams &params_)
    : saParams(saParams_), rMgrParams(params_)
{

    path_names[0] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info.txt";
    path_names[1] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info-register.txt";
    path_names[2] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info-graphics.txt";
}

RecordManager::~RecordManager()
{
    if (info.is_open())
        info.close();
    if (infoRegister.is_open())
        infoRegister.close();
    if (infoGraphics.is_open())
        infoGraphics.close();
}

void RecordManager::open_file(const std::size_t n_file, std::ofstream &file)
{
    try
    {
        file.open(path_names[n_file], std::ios::out | std::ios::app);
    }
    catch (const std::ofstream::failure &error)
    {
        std::cerr << "[Error]: No se pudo abrir el archivo: " << path_names[n_file] << std::endl;
    }
}

/*
 * Genera los archivos que contienen información de los estados de estudiantes y escuelas durante
 * la ejecución del algoritmo
 */
void RecordManager::openRecordInfo()
{
    open_file(0, info);
}

void RecordManager::openRecordRegister()
{
    open_file(1, infoRegister);
}

void RecordManager::openRecordGraphics()
{
    open_file(2, infoGraphics);
}

void RecordManager::closeRecordInfo()
{
    info.close();
}

void RecordManager::closeRecordRegister()
{
    infoRegister.close();
}

void RecordManager::closeRecordGraphics()
{
    infoGraphics.close();
}

/*
 * Escribe la información en los archivos previamente generados
 */
void RecordManager::SaveInfoInit(double costBestSolution,
                                 double meanDist,
                                 double S,
                                 double costCupo)
{
    info << "--------------- Primeros datos -------------\n";
    info << "Primer costo de solución: " << costBestSolution << "\n";
    info << "Primer distancia: " << meanDist << "\n";
    info << "Primer Segregación: " << S << "\n";
    info << "Primer CostoCupo: " << costCupo << "\n\n";
}

void RecordManager::SaveInfoFinish(
    double costPreviousSolution,
    double costBestSolution,
    double costCurrentSolution,
    double time_taken,
    double meanDist,
    double S,
    double costCupo)
{
    info << "--------------- Resultado Final ----------------"
         << "\n";
    info << "Numero de Ciclos: " << saParams.count << "\n";
    info << "Costo de la solución previa: " << costPreviousSolution << "\n";
    info << "Costo de la mejor solución: " << costBestSolution << "\n";
    info << "Costo de la solución actual: " << costCurrentSolution << "\n";
    info << "Tiempo de ejecución de SA: " << time_taken << "\n";
    info << "distancia: " << meanDist << "\n";
    info << "Segregación: " << S << "\n";
    info << "CostoCupo: " << costCupo << "\n";
    info << "--------------- Finalizo con exito ----------------"
         << "\n";
}

void RecordManager::SaveInfoRegisterInit()
{
}

void RecordManager::SaveInfoRegisterFinish()
{
}

void RecordManager::SaveGraphicsInit(double meanDist, double S, double costCupo, double costCurrentSolution)
{
    infoGraphics << std::setprecision(13);
    infoGraphics << saParams.count << ","
                 << meanDist / saParams.max_dist << ","             // Distancia promedio recorrida por los estudiantes normalizada
                 << meanDist << ","                                 // Distancia promedio recorrida por los estudiantes
                 << S << ","                                        // Indice de duncan
                 << costCupo << ","                                 // Costo cupo de las escuelas
                 << costCurrentSolution << ","                      // Solución actual
                 << saParams.temp << std::setprecision(13) << "\n"; // Temperatura actual
}

void RecordManager::SaveGraphicsFinish()
{
}
// void RecordManager::recordEnd()
// {

//   for (x = 0; x < n_students; x++)
//   {
//       info_graficos_bestSolution << bestSolution[x] << ",";
//   }

//   for (x = 0; x < vector_count.size(); x++)
//   {
//       info_graficos << vector_count.at(x) << ","
//                     << vector_meanDist.at(x) / max_dist << "," // Distancia promedio recorrida por los estudiantes normalizada
//                     << vector_meanDist.at(x) << ","
//                     << vector_segregation.at(x) << ","
//                     << vector_costoCupo.at(x) << ","
//                     << vector_costCurrentSolution.at(x) << ","
//                     << std::fixed << vector_temp.at(x) << std::setprecision(13) << "\n";
//   }

//   ///////////////////////////////////////////////////
//   /// Almacenamiento de datos
//   ///////////////////////////////////////////////////
//   cout.precision(dbl::max_digits10);
//   cout << "--------------- Resultado Final ----------------"
//        << "\n";
//   cout << "Numero de Ciclos " << count << "\n";
//   cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
//   cout << "Costo de la mejor solución: " << costBestSolution << "\n";
//   cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
//   cout << "Tiempo de ejecución de SA: " << time_taken << "\n";
//   cout << "distancia: " << meanDist(bestSolution, distMat) << "\n";
//   cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
//   cout << "CostoCupo: " << costCupo(bestSolution, cupoArray) << "\n";

//   cout << "Cal costo " << calCosto(bestSolution, distMat, ptr_alpha, alumnosSep, totalVuln, cupoArray) << endl;
//   cout << "Costo de: " << costBestSolution << "\n";

//   // cout << fixed << setprecision(70) << endl;
//   cout << sumDist(bestSolution, distMat) << "\n";
//   cout << bestVars[0] << endl;
//   cout << sumS(bestSolution, alumnosSep, totalVuln) << "\n";
//   cout << bestVars[1] << endl;
//   cout << sumCostCupo(bestSolution, cupoArray) << "\n";
//   cout << bestVars[2] << endl;
//   cout << "Tiempo de ejecución de SA get_result: " << vector_time1 << "\n";

//   cout << "--------------- Finalizo con exito ----------------"
//        << "\n";

//   info << "--------------- Resultado Final ----------------"
//        << "\n";
//   info << "Numero de Ciclos " << count << "\n";
//   info << "Costo de la solución previa: " << costPreviousSolution << "\n";
//   info << "Costo de la mejor solución: " << costBestSolution << "\n";
//   info << "Costo de la solución actual: " << costCurrentSolution << "\n";
//   info << "Tiempo de ejecución de SA: " << time_taken << "\n";
//   info << "distancia: " << meanDist(bestSolution, distMat) << "\n";
//   info << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
//   info << "CostoCupo: " << costCupo(bestSolution, cupoArray) << "\n";
//   info << "--------------- Finalizo con exito ----------------"
//        << "\n";

//   info_test << fixed << time_taken << setprecision(9) << ","
//             << costBestSolution << ","
//             << meanDist(bestSolution, distMat) / max_dist
//             << "," << meanDist(bestSolution, distMat)
//             << "," << S(bestSolution, alumnosSep, totalVuln)
//             << "," << costCupo(bestSolution, cupoArray)
//             << "," << count
//             << "," << fixed << temp_init << setprecision(13)
//             << "," << fixed << cooling.getTemp() << setprecision(13)
//             << "," << min_temp
//             << "," << seed
//             << "," << alpha1
//             << "," << alpha2
//             << "," << alpha3
//             << "," << alpha[0]
//             << "," << alpha[1]
//             << "," << alpha[2]
//             << "," << coolingRate
//             << "," << k_reheating_init
//             << "," << e_const
//             << "," << n_reheating
//             << "," << len1_init
//             << "," << len2_init
//             << "," << len3_init
//             << "," << len4_init
//             << "," << len1
//             << "," << len2
//             << "," << len3
//             << "," << len4
//             << "," << Th
//             << "," << n_block
//             << "," << n_thread
//             << "," << name_exp << "\n";

//   info_graficos_bestSolution.close();
//   cout << ".";
//   info_graficos.close();
//   cout << ".";
//   info_test.close();
//   info.close();
//   cout << ".\n";
//   cout << " Archivos Guardado"
//        << "\n";
// }
