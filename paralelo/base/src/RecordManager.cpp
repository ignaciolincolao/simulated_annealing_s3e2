#include <RecordManager.hpp>

#include <iomanip>

RecordManager::RecordManager(SimulatedParams &saParams_, RecordParams &params_)
    : saParams(saParams_), rMgrParams(params_)
{

    path_names[0] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info.txt";
    path_names[1] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info-register.txt";
    path_names[2] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info-graphics.txt";
    path_names[3] = rMgrParams.ruta_save + rMgrParams.prefijo_save + "-info-graphicsBestSolution.txt";
}

RecordManager::~RecordManager()
{
    if (info.is_open())
        info.close();
    if (infoRegister.is_open())
        infoRegister.close();
    if (infoGraphics.is_open())
        infoGraphics.close();
    if (infoGraphicsBestSolution.is_open())
        infoGraphicsBestSolution.close();
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

void RecordManager::openRecordGraphicsBestSolution()
{
    open_file(3, infoGraphicsBestSolution);
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

void RecordManager::closeRecordGraphicsBestSolution()
{
    infoGraphicsBestSolution.close();
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

void RecordManager::SaveGraphicsBestSolution(int *solution)
{
    for (std::size_t i{}; i < saParams.n_students; i++)
        infoGraphicsBestSolution << solution[i] << ",";
    infoGraphicsBestSolution << "\n";
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
    for (std::size_t x = 0; x < vector_count.size(); x++)
    {
        infoGraphics << vector_count.at(x) << ","
                     << vector_meanDist.at(x) / saParams.max_dist << "," // Distancia promedio recorrida por los estudiantes normalizada
                     << vector_meanDist.at(x) << ","
                     << vector_segregation.at(x) << ","
                     << vector_costoCupo.at(x) << ","
                     << vector_costCurrentSolution.at(x) << ","
                     << std::fixed << vector_temp.at(x) << std::setprecision(13) << "\n";
    }
}

void RecordManager::SaveInfoRegister(
    double time_taken,
    double costBestSolution,
    double meanDist,
    double S,
    double costCupo,
    double coolingRate,
    double k_reheating_init,
    double e_const,
    int n_reheating,
    int len1_init,
    int len2_init,
    double len3_init,
    double len4_init,
    int len1,
    int len2,
    double len3,
    double len4,
    double Th,
    int n_block,
    int n_thread)
{
    infoRegister << std::fixed << time_taken << std::setprecision(9) << ","
                 << costBestSolution << ","
                 << meanDist / saParams.max_dist
                 << "," << meanDist
                 << "," << S
                 << "," << costCupo
                 << "," << saParams.count
                 << "," << std::fixed << saParams.temp_init << std::setprecision(13)
                 << "," << std::fixed << saParams.temp << std::setprecision(13)
                 << "," << saParams.min_temp
                 << "," << saParams.seed
                 << "," << saParams.alpha1
                 << "," << saParams.alpha2
                 << "," << saParams.alpha3
                 << "," << saParams.alpha[0]
                 << "," << saParams.alpha[1]
                 << "," << saParams.alpha[2]
                 << "," << coolingRate
                 << "," << k_reheating_init
                 << "," << e_const
                 << "," << n_reheating
                 << "," << len1_init
                 << "," << len2_init
                 << "," << len3_init
                 << "," << len4_init
                 << "," << len1
                 << "," << len2
                 << "," << len3
                 << "," << len4
                 << "," << Th
                 << "," << n_block
                 << "," << n_thread
                 << "," << rMgrParams.name_exp << "\n";
}