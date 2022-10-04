
#include <sas.hpp>
#include <ExplorationCriterion.hpp>
#include <AcceptanceCriterion.hpp>
#include <TemperatureLength.hpp>
#include <ReheatingMethods.hpp>
#include <CoolingScheme.hpp>

#include <limits>

#define DECIMAL 16
typedef std::numeric_limits< long double > dbl;

///////////////////////////////////////////////////
/// Variables globales.
///////////////////////////////////////////////////




long double sasFunc() {
    int x=0,z=0;
    int totalVuln=0;
    cout.precision(dbl::max_digits10);
    //cout << fixed << setprecision(70) << endl;
    //srand(time(NULL));
    ///////////////////////////////////////////////////
    /// Genera archivo de almacenamiento de datos
    ///////////////////////////////////////////////////

    /*
    * Prepara archivos para guardar los datos
    */
    
    ofstream info;
    string infotxt = ruta_save + prefijo_save +"-info.txt"; 
    info.open(infotxt);
    /*
    * 
    */
    ofstream info_test;
    string nameinfo_test = ruta_save + prefijo_save+"-info-test.txt"; 
    info_test.open(nameinfo_test);
    /*
    * Genera los archivos que contienen información de los estados de estudiantes y escuelas durante
    * la ejecución del algoritmo
    */
    ofstream info_graficos;
    string name_info_graficos = ruta_save + prefijo_save +"-info-graficos.txt";
    info_graficos.open(name_info_graficos);

    ///////////////////////////////////////////////////
    /// Datos colegios
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura Info_colegio
    ///////////////////////////////////////////////////
    Info_colegio *ptr_colegios;
    vector<Info_colegio> colegios;
    getDataSchool(colegios);
    ptr_colegios = colegios.data();
    n_colegios = colegios.size();

    ///////////////////////////////////////////////////
    /// Datos Alumnos
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura info_student
    ///////////////////////////////////////////////////
    Info_alu *ptr_students;
    vector<Info_alu> students;
    getDataStudents(students,totalVuln);
    ptr_students = students.data();
    n_students = students.size();

    ///////////////////////////////////////////////////
    /// Inicializa Variables y arreglos
    ///////////////////////////////////////////////////

    int aluVulxCol[n_colegios], aluxcol[n_colegios];
    int previousAluxCol[n_colegios];
    int previousAluVulxCol[n_colegios];
    int bestAluxCol[n_colegios];
    int bestAluVulxCol[n_colegios];

    int *previousSolution= nullptr;
    int *bestSolution= nullptr;
    int *currentSolution=nullptr;
    long double **distMat=nullptr;
    int *cupoArray=nullptr;
    int *alumnosSep=nullptr;

    
    long double  costBestSolution,
        costPreviousSolution,
        costCurrentSolution,
        *ptr_alpha = &alpha[0];
    
    int count=0;

    previousSolution = (int *)malloc(sizeof(int)*n_students);
    bestSolution=(int *)malloc(sizeof(int)*n_students);
    currentSolution=(int *)malloc(sizeof(int)*n_students);
    distMat=(long double **)malloc(sizeof(long double)*n_students);
    cupoArray=(int *)malloc(sizeof(int)*n_colegios);
    alumnosSep = (int *)malloc( sizeof(int)*n_students);
    for(x=0; x < n_students; x++) {
        distMat[ x ]=(long double *)malloc(sizeof(long double)*n_colegios);
    }

    ///////////////////////////////////////////////////
    /// Se asignan las escuelas un arreglo que y estudiantes a la escuela
    /// las escuelas tendran como identificación el indice
    /// y currentSolution tiene como indice al estudiante y el valor del indice a la escuela que asignada
    ///////////////////////////////////////////////////
    initializeArray(aluxcol, 
                previousAluxCol, 
                bestAluxCol, 
                aluVulxCol, 
                previousAluVulxCol, 
                bestAluVulxCol, 
                alumnosSep,
                students,
                colegios);
    assignSchoolToArray(previousSolution, bestSolution, currentSolution, ptr_colegios, ptr_students, cupoArray);
    calcDist(ptr_colegios, ptr_students, distMat);
    max_dist = getMaxDistance(distMat);
    normalizedAlpha(alpha);

    ///////////////////////////////////////////////////
    /// Registro de datos
    ///////////////////////////////////////////////////
    costBestSolution=calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
    costPreviousSolution=costBestSolution;
    costCurrentSolution=costBestSolution;
    
    
    cout << "--------------- Primeros datos -------------" << "\n";
    cout << "Primer costo de solución: " << costBestSolution << "\n";
    cout << "Primer distancia: " << meanDist(currentSolution,distMat) << "\n";
    cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    cout << "Primer CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";

    info      << "--------------- Primeros datos -------------" << "\n";
    info      << "Primer costo de solución: " << costBestSolution << "\n";
    info      << "Primer distancia: " << meanDist(currentSolution,distMat) << "\n";
    info      << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    info      << "Primer CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";


    ///////////////////////////////////////////////////
    /// Generación de archivos que almacenan información de los graficos
    ///////////////////////////////////////////////////

    info_graficos << setprecision(13);
    info_graficos << count << "," 
                << meanDist(currentSolution,distMat)/max_dist << "," // Distancia promedio recorrida por los estudiantes normalizada
                << meanDist(currentSolution,distMat) << "," // Distancia promedio recorrida por los estudiantes
                << S(currentSolution, alumnosSep, totalVuln) << "," // Indice de duncan
                << costCupo(currentSolution,cupoArray) << "," // Costo cupo de las escuelas
                << costCurrentSolution << "," // Solución actual
                << temp << setprecision(13) << "\n"; // Temperatura actual

    count++;
    ///////////////////////////////////////////////////
    /// Genera arreglos que contendran valores del 0 hasta n_students y n_colegios
    ///////////////////////////////////////////////////
    int *shuffle_student = (int *)malloc(sizeof(int)*n_students);
    int *shuffle_colegios = (int *)malloc(sizeof(int)*n_colegios);
    for (int i = 0; i < n_students; i++) {
        shuffle_student[i] = i;
    }
    for (int i=0; i < n_colegios; i++){
        shuffle_colegios[i]=i;
    }
    ///////////////////////////////////////////////////
    /// Posicion estudiantes
    ///////////////////////////////////////////////////

    ofstream info_graficos_bestSolution;
    string name_info_graficos_bestSolution = ruta_save + prefijo_save +"-info-graficos_bestSolution.txt"; // concatenar
    info_graficos_bestSolution.open(name_info_graficos_bestSolution);
    for(x=0;x<n_students;x++){
        info_graficos_bestSolution << currentSolution[x] << ",";
    }
    info_graficos_bestSolution << "\n";

    ///////////////////////////////////////////////////
    /// Genera distribuciones para seleccionar un estudiante y una escuela al azar
    ///////////////////////////////////////////////////

    dist = uniform_int_distribution<int>(0, n_students-1);
    dist2 = uniform_int_distribution<int>(0, n_colegios-1);

    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////
    auto start = chrono::high_resolution_clock::now();
    ///////////////////////////////////////////////////
    /// Comienza a ejecutarse el algoritmo de SA
    ///////////////////////////////////////////////////


    vector<long double> vector_costCurrentSolution;
    vector<long double> vector_meanDist;
    vector<long double> vector_segregation;
    vector<long double> vector_costoCupo;
    vector<long double> vector_temp;
    vector<int> vector_count;

    
    int count_rechaso=0;
    int reheating = 0;
    int c_accepta = 0;
    int c_cooling_temperature = 0;
    int valmaxheating=n_colegios;
    int count_reheating = 0;
    long double bestTemp = 0;
    long double k_reheating_init = k_reheating;
    long double temp_init = temp;
    int count_trials = 0;
    float len1_init = len1;
    float len2_init = len2;
    long double len3_init = len3;
    long double len4_init = len4;

    ////////////////////////////////////////////////////////////////////////
    // VARIABLES DE PRUEBA
    ////////////////////////////////////////////////////////////////////////

    long double costCurrentSolutionV2 = costCurrentSolution;
    long double *currentVars = (long double *)malloc(3 * sizeof(long double));
    long double *previousVars = (long double *)malloc(3 * sizeof(long double));
    long double *bestVars = (long double *)malloc(3 * sizeof(long double));

    currentVars[0] = sumDist(currentSolution,distMat);
    currentVars[1] = sumS(currentSolution, alumnosSep, totalVuln);
    currentVars[2] = sumCostCupo(currentSolution,cupoArray);
    previousVars[0] = currentVars[0];
    previousVars[1] = currentVars[1];
    previousVars[2] = currentVars[2];


    auto start_compare = chrono::high_resolution_clock::now();
    auto end_compare = chrono::high_resolution_clock::now();
    long double time_taken_v1 = chrono::duration_cast<chrono::nanoseconds>(end_compare - start_compare).count();
    long double time_taken_v2;
    long double vector_time1 =0;
    long double vector_time2=0;
    
    while(temp > min_temp){

        memcpy(currentSolution,previousSolution,sizeof(int)*n_students);
        memcpy(aluxcol,previousAluxCol,sizeof(int)*n_colegios);
        memcpy(aluVulxCol,previousAluVulxCol,sizeof(int)*n_colegios);
        memcpy(currentVars,previousVars,sizeof(long double)*3);
        ////////////////////////////////

        ///////////////////////////////////////////////////
        ///  Selecciona aleatoria mente a los alumnos
        ///////////////////////////////////////////////////
        
        start_compare = chrono::high_resolution_clock::now();
        costCurrentSolution = solutionNE1_v2(n_students,n_colegios,totalVuln,aluxcol,aluVulxCol,cupoArray,distMat,currentSolution,costCurrentSolution,ptr_alpha,shuffle_student,shuffle_colegios,alumnosSep,currentVars, max_dist);
        end_compare = chrono::high_resolution_clock::now();
        time_taken_v1 = chrono::duration_cast<chrono::nanoseconds>(end_compare - start_compare).count();
        time_taken_v1 *= 1e-9;
        vector_time1+=time_taken_v1;
        
        
        // Verifica si el costo actual es mejor que la mejor solución
        // en el caso que el costo actual es menor a la mejor solución, acepta la solución y los
        // guarda en el estado como mejor solución




        if(costCurrentSolution < costBestSolution){

            memcpy(bestSolution,currentSolution,sizeof(int)*n_students);
            memcpy(previousSolution,currentSolution,sizeof(int)*n_students);
            memcpy(previousAluxCol,aluxcol,sizeof(int)*n_colegios);
            memcpy(previousAluVulxCol,aluVulxCol,sizeof(int)*n_colegios);
            memcpy(previousVars,currentVars,sizeof(long double)*3);
            memcpy(bestVars,currentVars,sizeof(long double)*3);


            costBestSolution=costCurrentSolution;
            costPreviousSolution=costCurrentSolution;

            
            vector_costCurrentSolution.push_back(costCurrentSolution);
            vector_meanDist.push_back(meanDist(currentSolution,distMat));
            vector_segregation.push_back(S(currentSolution, alumnosSep, totalVuln));
            vector_costoCupo.push_back(costCupo(currentSolution,cupoArray));
            vector_temp.push_back(temp);
            vector_count.push_back(count);

            c_accepta++;
            count_rechaso=0;
        }
        // En el caso que el la solución actual sea mas alta intenta aceptar una peor solución en base
        // a la función acepta
        else{
            if(metropolisAC1(costPreviousSolution,costCurrentSolution)==1){
                memcpy(previousSolution,currentSolution,sizeof(int)*n_students);
                memcpy(previousAluxCol,aluxcol,sizeof(int)*n_colegios);
                memcpy(previousAluVulxCol,aluVulxCol,sizeof(int)*n_colegios);
                memcpy(previousVars,currentVars,sizeof(long double)*3);

                costPreviousSolution=costCurrentSolution;

                count_rechaso=0;
                c_accepta++;
            }
            else{
                count_rechaso++;
            }
        }

        if(temperatureTL7(temp, c_cooling_temperature, c_accepta, len1, len2, n_colegios, coolingRate,count)){
        //if(temperatureTL8(temp, c_cooling_temperature, count_trials, len1, len2, coolingRate)){
        //if(temperatureTL9(temp, c_cooling_temperature, count_trials, len3, len4, coolingRate)){
        //if(temperatureTL11(temp, c_cooling_temperature, count_trials, len3, len4, coolingRate)){
            coolingCS2(temp,coolingRate);
        }

        //reheatingTR11(temp, k_reheating, n_reheating, count_rechaso);
        //reheatingTR12(temp, k_reheating, n_reheating, count);
        //reheatingTR13(temp, k_reheating, n_reheating, c_cooling_temperature);
        //reheatingTR14(temp, k_reheating, k_reheating_init, n_reheating, count_rechaso, e_const);
        count_trials++;
        count++;

    }

    ///////////////////////////////////////////////////
    /// Obtiene el tiempo de ejecución
    ///////////////////////////////////////////////////
    auto end = chrono::high_resolution_clock::now();
    long double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    for(x=0;x<n_students;x++){
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

    cout << "Cal costo" << calCosto(bestSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray) << endl;
    cout << "Costo de:" << costBestSolution << "\n";


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
            << "," << fixed << temp << setprecision(13) 
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


    return (costBestSolution);

}

///////////////////////////////////////////////////
///////////////////////////////////////////////////


///////////////////////////////////////////////////
/// Calcula el costo
///////////////////////////////////////////////////
long double calCosto(int currentSolution[], long double **distMat, const long double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]){
    long double var1 = meanDist(currentSolution,distMat)/max_dist;
    //cout << "distancia: " << var1 << "\n";
    long double var2 = S(currentSolution, alumnosSep, totalVuln);
    //cout << "Segregación: " << var2 << "\n";
    long double var3 = costCupo(currentSolution,cupoArray);
    //cout << "CostoCupo: " << var3 << "\n";
    return (long double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
}

///////////////////////////////////////////////////
/// Distancia promedio que recorreran los estudiantes
///////////////////////////////////////////////////
long double meanDist(const int currentSolution[], long double  **distMat){
    long double sumDist=0.0;
    for(int i=0;i<n_students;i++){
        sumDist+=round_n(distMat[i][currentSolution[i]]); // distMat[estudiante][escuela]
    }
    //cout << "meanDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << n_student << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist/n_students;
}

long double sumDist(const int currentSolution[], long double  **distMat){
    long double sumDist=0.0;
    for(int i=0;i<n_students;i++){
        sumDist+=round_n(distMat[i][currentSolution[i]]); // distMat[estudiante][escuela]
    }
    //cout << "sumDist: " << sumDist << endl;
    //cout << "Numero de estudiantes: " << n_student << "  |  Suma de distancias:" << sumDist << "\n";
    return sumDist;
}


///////////////////////////////////////////////////
/// Calcula segregación por duncan
///////////////////////////////////////////////////

long double S(const int currentSolution[],const int alumnosSep[], int totalVuln){
    long double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=round_n(fabs((aluVulCol/(long double)totalVuln)-(aluNoVulCol/(long double)(n_students-totalVuln))));
        }
    }
    return totalSesc/2.0;
}

long double sumS(const int currentSolution[],const int alumnosSep[], int totalVuln){
    long double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=round_n(fabs((aluVulCol/(long double)totalVuln)-(aluNoVulCol/(long double)(n_students-totalVuln))));
        }
    }
    return totalSesc;
}


///////////////////////////////////////////////////
/// Calcula el costo de tener los estudiantes en las escuelas
///////////////////////////////////////////////////

long double costCupo(int currentSolution[],int cupoArray[]){
    long double totalcostCupo = 0.0;
    int totalAluCol = 0;
    long double a = 0.0;
    for(int j=0;j<n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        totalcostCupo+=round_n((long double)totalAluCol*fabs(((long double)cupoArray[j]-totalAluCol)/pow(((long double)cupoArray[j]/2),2)));
    }
    return totalcostCupo/n_colegios;
}



long double sumCostCupo(int currentSolution[],int cupoArray[]){
    long double totalcostCupo = 0.0;
    int totalAluCol = 0;
    for(int j=0;j<n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        totalcostCupo+= round_n((long double)totalAluCol*fabs(((long double)cupoArray[j]-totalAluCol)/pow(((long double)cupoArray[j]/2),2)));
    }
    return totalcostCupo;
}
///////////////////////////////////////////////////
/// Genera una nueva solución en donde asigna a un estudiante a una escuela
/// aleatoriamente
///////////////////////////////////////////////////

void newSolution(int currentSolution[],const int previousSolution[]){
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> dist(0, n_students);
    random_device rd2;
    mt19937 mt2(rd2());
    uniform_int_distribution<int> dist2(0, n_colegios);
    int selectStudent=dist(mt);
    int selectSchool = dist2(mt2);
    for(int x=0; x<n_students; x++){
        if(x == selectStudent) {
            currentSolution[x] = selectSchool;
        }
        else {
            currentSolution[x] = previousSolution[x];
        }
    }

}




///////////////////////////////////////////////////
/// Asigna a las soluciones la escuela actual Solo se utiliza al inicio
///////////////////////////////////////////////////
void assignSchoolToArray(int previousSolution[], int bestSolution[], int currentSolution[], Info_colegio *ptr_colegios, Info_alu *ptr_students, int cupoArray[]){
    Info_alu *ptr_aux = ptr_students;
    for(int x=0;x < n_colegios;x++){
        for(int y=0; y < n_students; y++){
            if(ptr_colegios->rbd == ptr_students->rbd){
                previousSolution[y] = x;
                bestSolution[y] = x;
                currentSolution[y] = x;
            }
            ptr_students++;

        }
        /*
         * cupoArray sera un arreglo que por indice es la escuela y su valor sera el cupo que posee esa escuela
         * se asume que las escuelas pueden tener sobre cupo.
         */

        cupoArray[x] = ptr_colegios->num_alu+ ((int)((ptr_colegios->num_alu*10)/100));

        ptr_students = ptr_aux;
        ptr_colegios++;
    }
}
///////////////////////////////////////////////////
/// Crea una matriz de distancia donde x es el estudiante, y es la escuela
///////////////////////////////////////////////////
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, long double **distMat){
    Info_colegio *ptr_aux = ptr_colegios;
    for(int x=0;x < n_students ;x++){
        for(int y=0; y < n_colegios; y++){
            distMat[x][y] = sqrt( pow((ptr_students->latitude - ptr_colegios->latitude),2)+pow((ptr_students->longitude - ptr_colegios->longitude),2))/1000;
            ptr_colegios++;

        }
        ptr_colegios = ptr_aux;
        ptr_students++;
    }
}

///////////////////////////////////////////////////
/// newSolution_v2, tiene como entrada la información de los estado actual de la solución, y alcula de inmediato la
/// distancia promedio, el costocupo y segregación total.
///////////////////////////////////////////////////
long double newSolution_v2(int n_students,int n_colegios,int totalVuln,int aluxcol[],int aluVulxCol[],int cupoArray[],long double **distMat, int currentSolution[], const long double ptr_alpha[]){
    long double sumDist=0;
    long double mean=0.0;
    long double totalcostCupo = 0.0;
    long double totalSesc = 0.0;
    int aluVulCol, aluNoVulCol,totalAluCol;
    for(int i=0;i<n_students;i++){
        sumDist+=distMat[i][currentSolution[i]]; // distMat[estudiante][escuela]
    }
    for(int n=0; n<n_colegios; n++) {
        totalAluCol = aluxcol[n];
        aluVulCol = aluVulxCol[n];
        aluNoVulCol =totalAluCol - aluVulCol;
        // Calcula el costo cupo
        totalcostCupo+=totalAluCol*fabs((cupoArray[n]-totalAluCol)/pow(((long double)cupoArray[n]/2),2));
        // Calcula el total sesc
        totalSesc+=((long double)1/2)*fabs((aluVulCol/(long double)totalVuln)-(aluNoVulCol/(long double)(n_students-totalVuln)));
    }
    totalcostCupo = totalcostCupo/n_colegios;
    long double var1 = (sumDist/(double)n_students)/max_dist;
    //cout << var1 << "\n";
    long double var2 = totalSesc;
    //cout << var2 << "\n";
    long double var3 = totalcostCupo;
    //cout << var3 << "\n";
    return (long double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
}

void shuffle(int *values, const int max_change, uniform_int_distribution<int> distri) {
    int randvalue1,randvalue2,tem_value;
    for (int i = 0; i<max_change; i++) {
        randvalue1 = distri(mt);
        randvalue2 = i;
        tem_value = values[randvalue1];
        values[randvalue1] = values[randvalue2];
        values[randvalue2] = tem_value;
    }
}

void getDataSchool(vector<Info_colegio> &colegios){
    string line_colegios;
    ifstream info_school("colegios_utm.txt"); // concatenar
    int cx = 0;
    while (getline(info_school, line_colegios)) {
        stringstream linestream(line_colegios);
        string data;
        colegios.push_back(Info_colegio());
        getline(linestream, data, ',');
        colegios[cx].rbd = stoi(data);
        getline(linestream, data, ',');
        colegios[cx].latitude = stod(data);
        getline(linestream, data, ',');
        colegios[cx].longitude = stod(data);
        getline(linestream, data, ',');
        colegios[cx].num_alu = stoi(data);
        getline(linestream, data, ',');
        colegios[cx].prioritario = stoi(data);
        cx++;
    }
    info_school.close();
}

void getDataStudents(vector<Info_alu> &students, int &totalVuln)
{
    string line_student;
    ifstream info_student("alumnos_utm.txt"); // concatenar
    int cx = 0;
    while (getline(info_student, line_student)) {
        stringstream linestream(line_student);
        string data;
        students.push_back(Info_alu());
        getline(linestream, data, ',');
        students[cx].rbd = stoi(data);
        getline(linestream, data, ',');
        students[cx].latitude = stod(data);
        getline(linestream, data, ',');
        students[cx].longitude = stod(data);
        getline(linestream, data, ',');
        students[cx].sep = stoi(data);
        if (students[cx].sep == 1) {
            totalVuln++;
        }
        cx++;

    }
    info_student.close();
}

////////////////////////////////////////////////
////// Obtiene la maxima distancia que un estudiante podria llegar a recorrer
///////////////////////////////////////////////////
long double getMaxDistance(long double **distMat){
    long double max = 0;
    for(int i=0;i<n_students;i++){
        for(int x=0;x<n_colegios;x++){
            if(distMat[i][x]>max){
                max = distMat[i][x];
            }
        }
    }
    return max;
}

///////////////////////////////////////////////////
/// Calcula el valor de los alpha
///////////////////////////////////////////////////
void normalizedAlpha(long double alpha[3])
{
    long double sumaAlpha = 0.0;
    for(int x=0; x<3; x++){
        sumaAlpha +=alpha[x];
    }
    for(int x=0; x<3; x++){
        alpha[x]= alpha[x]/(long double)sumaAlpha;
    }
}


///////////////////////////////////////////////////
/// Asigna Información de las escuelas a best, previus y current soluciones
///////////////////////////////////////////////////
void initializeArray(int *aluxcol, int *previousAluxCol, int *bestAluxCol, int *aluVulxCol, int *previousAluVulxCol, int *bestAluVulxCol, int *alumnosSep, vector<Info_alu> &students,vector<Info_colegio> &colegios)
{
    for(int x = 0; x < n_colegios; x++){
        aluxcol[x] = colegios[x].num_alu;
        previousAluxCol[x] = colegios[x].num_alu;
        bestAluxCol[x] = colegios[x].num_alu;
        aluVulxCol[x] = colegios[x].prioritario;
        previousAluVulxCol[x] = colegios[x].prioritario;
        bestAluVulxCol[x] = colegios[x].prioritario;

    }
    ///////////////////////////////////////////////////
    /// Se crear un arreglo donde el el valor es la posición del estudiante sep
    ///////////////////////////////////////////////////
    for(int x=0; x < n_students; x++) {
        alumnosSep[x] = students[x].sep;
    }
}


long double round_n(long double x)
{
    long double digits = pow(10.0, DECIMAL);
    return trunc(x * digits) / digits;
}


