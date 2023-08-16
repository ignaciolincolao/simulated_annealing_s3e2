#include <RecordManager.hpp>

RecordManager::RecordManager(const std::string& ruta_save, const std::string& prefijo_save){
    std::string infotxt = ruta_save + prefijo_save + "-info.txt";
    info.open(nameInfo);

    std::string nameinfo_test = ruta_save + prefijo_save + "-info-test.txt";
    infoRegister.open(nameInfoRegister);

    std::string name_info_graficos = ruta_save + prefijo_save + "-info-graficos.txt";
    infoGraphics.open(nameInfoGraphics);
}

RecordManage::~RecordManager(){
    if (info.is_open()) {
        info.close();
    }
    if (infoRegister.is_open()) {
        infoRegister.close();
    }
    if (infoGraphics.is_open()) {
        infoGraphics.close();
    }
}

void RecordManage::recordEnd(){

    for(x = 0; x < n_students; x++) {
        info_graficos_bestSolution << bestSolution[x] << ",";
    }
    

    for(x=0; x<vector_count.size(); x++){
        info_graficos << vector_count.at(x) << "," 
                    << vector_meanDist.at(x)/max_dist << "," // Distancia promedio recorrida por los estudiantes normalizada
                    << vector_meanDist.at(x) << "," 
                    << vector_segregation.at(x) << "," 
                    << vector_costoCupo.at(x) << "," 
                    << vector_costCurrentSolution.at(x) << "," 
                    << fixed << vector_temp.at(x) << setprecision(13) << "\n";
    }




    ///////////////////////////////////////////////////
    /// Almacenamiento de datos
    ///////////////////////////////////////////////////
    cout.precision(dbl::max_digits10);
    cout << "--------------- Resultado Final ----------------" << "\n";
    cout << "Numero de Ciclos " << count << "\n";
    cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
    cout << "Costo de la mejor solución: " << costBestSolution << "\n";
    cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
    cout << "Tiempo de ejecución de SA: " << time_taken << "\n";
    cout << "distancia: " << meanDist(bestSolution,distMat) << "\n";
    cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << "CostoCupo: " << costCupo(bestSolution,cupoArray) << "\n";

    cout << "Cal costo " << calCosto(bestSolution, distMat, ptr_alpha, alumnosSep, totalVuln, cupoArray) << endl;
    cout << "Costo de: " << costBestSolution << "\n";


    //cout << fixed << setprecision(70) << endl;
    cout << sumDist(bestSolution,distMat) << "\n";
    cout << bestVars[0] << endl;
    cout << sumS(bestSolution, alumnosSep, totalVuln) << "\n";
    cout << bestVars[1] << endl;
    cout << sumCostCupo(bestSolution,cupoArray) << "\n";
    cout << bestVars[2] << endl;
    cout << "Tiempo de ejecución de SA get_result: " << vector_time1 << "\n";

    cout << "--------------- Finalizo con exito ----------------" << "\n";


    info << "--------------- Resultado Final ----------------" << "\n";
    info << "Numero de Ciclos " << count << "\n";
    info << "Costo de la solución previa: " << costPreviousSolution << "\n";
    info << "Costo de la mejor solución: " << costBestSolution << "\n";
    info << "Costo de la solución actual: " << costCurrentSolution << "\n";
    info << "Tiempo de ejecución de SA: " << time_taken << "\n";
    info << "distancia: " << meanDist(bestSolution,distMat) << "\n";
    info << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    info << "CostoCupo: " << costCupo(bestSolution,cupoArray) << "\n";
    info << "--------------- Finalizo con exito ----------------" << "\n";


    info_test << fixed << time_taken << setprecision(9) << "," 
            << costBestSolution << "," 
            << meanDist(bestSolution,distMat)/max_dist 
            << "," << meanDist(bestSolution,distMat) 
            << "," << S(bestSolution, alumnosSep, totalVuln) 
            << "," << costCupo(bestSolution,cupoArray) 
            << "," << count 
            << "," << fixed << temp_init << setprecision(13) 
            << "," << fixed << cooling.getTemp() << setprecision(13) 
            << "," << min_temp 
            << "," << seed
            << "," << alpha1 
            << "," << alpha2 
            << "," << alpha3 
            << "," << alpha[0]
            << "," << alpha[1]
            << "," << alpha[2]
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
            << ","<< name_exp << "\n";

    info_graficos_bestSolution.close();
    cout << ".";
    info_graficos.close();
    cout << ".";
    info_test.close();
    info.close();
    cout << ".\n";
    cout << " Archivos Guardado" << "\n";
}

void RecordManager::RecordInfo(){
}
void RecordManager::RecordRegister(){
}
/*
* Genera los archivos que contienen información de los estados de estudiantes y escuelas durante
* la ejecución del algoritmo
*/
void RecordManager::RecordGraphics(){
    ///////////////////////////////////////////////////
    /// Posicion estudiantes
    ///////////////////////////////////////////////////

    ofstream info_graficos_bestSolution;
    string name_info_graficos_bestSolution = ruta_save + prefijo_save +"-info-graficos_bestSolution.txt"; // concatenar
    info_graficos_bestSolution.open(name_info_graficos_bestSolution);
    for(x = 0; x < n_students; x++){
        info_graficos_bestSolution << currentSolution[x] << ",";
    }
    info_graficos_bestSolution << "\n";
}

void RecordManager::SaveInfo(){};
void RecordManager::SaveInfoRegister(){};
void RecordManager::SaveGraphics(){};
