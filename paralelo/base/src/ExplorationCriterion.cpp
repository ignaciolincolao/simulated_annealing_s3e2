#include <ExplorationCriterion.hpp>


double solutionNE1(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep)
    {
    int aluchange,colchange;
    shuffle(shuffle_student,1,dist);
    shuffle(shuffle_colegios,1,dist2);
    aluchange=shuffle_student[0];
    

    if(shuffle_colegios[0] != currentSolution[aluchange]){ 
        colchange = shuffle_colegios[0];
    }
    else{
        colchange = shuffle_colegios[1];
    }
    //ELimina el estudiante de la escuela actual
    aluxcol[currentSolution[aluchange]]-=1;
    aluVulxCol[currentSolution[aluchange]]-=alumnosSep[aluchange];
    //Asigna al estudiante a la nueva escuela
    currentSolution[aluchange] = colchange;
    aluxcol[colchange]+=1;
    aluVulxCol[colchange]+=alumnosSep[aluchange];


        // Obtiene el costo Actual

    return newSolution_v2(n_students,n_colegios,totalVuln,aluxcol,aluVulxCol,cupoArray,distMat,currentSolution,ptr_alpha);
}

double solutionNE1_v2(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep,
        double *currentVars,
        double max_dist)
    {
    int aluchange,colchange, currentSchool, newSchool,oldSchool;
    shuffle(shuffle_student,1,dist);
    shuffle(shuffle_colegios,1,dist2);
    aluchange=shuffle_student[0];
    

    if(shuffle_colegios[0] != currentSolution[aluchange]){ 
        colchange = shuffle_colegios[0];
    }
    else{
        colchange = shuffle_colegios[1];
    }

    double sumDist= 0.0;
    double totalSesc =0.0;
    double totalcostCupo =0.0;

    sumDist= currentVars[0];
    totalSesc = currentVars[1];
    totalcostCupo = currentVars[2];

    int aluVulCol, aluNoVulCol,totalAluCol;

    currentSchool = currentSolution[aluchange];
    newSchool = colchange;
    
    ////////////////////////////////////////////////////////////////
    /////// Descuenta antes de mover
    ////////////////////////////////////////////////////////////////
    // Distancia

    sumDist-=round_n(distMat[aluchange][currentSchool]);

    // seg de la escuela actual
    totalAluCol = aluxcol[currentSchool];
    //cout << "Alumnos actual escuela "<< totalAluCol << " " << endl;
    aluVulCol = aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela actual
    
    totalcostCupo-=round_n((double)totalAluCol*fabs(((double)cupoArray[currentSchool]-totalAluCol)/pow(((double)cupoArray[currentSchool]/2),2)));


    // seg de la escuela nueva
    totalAluCol = aluxcol[newSchool];
    //cout << "Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc-=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela nueva
    totalcostCupo-=round_n((double)totalAluCol*fabs(((double)cupoArray[newSchool]-totalAluCol)/pow(((double)cupoArray[newSchool]/2),2)));

    //ELimina el estudiante de la escuela actual
    aluxcol[currentSchool]-=1;
    aluVulxCol[currentSchool]-=alumnosSep[aluchange];
    //Asigna al estudiante a la nueva escuela
    currentSolution[aluchange] = newSchool;
    aluxcol[newSchool]+=1;
    aluVulxCol[newSchool]+=alumnosSep[aluchange];

    ////////////////////////////////////////////////////////////////
    ////// Calculó despues de mover
    //////////////////////////////////////////////////////////////
    sumDist+=round_n(distMat[aluchange][newSchool]);
    // seg de la escuela actual
    totalAluCol = aluxcol[currentSchool];
    //cout << "movio Alumnos nueva escuela "<< totalAluCol << " " << endl;
    aluVulCol = aluVulxCol[currentSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela actual
    totalcostCupo+=round_n((double)totalAluCol*fabs(((double)cupoArray[currentSchool]-totalAluCol)/pow(((double)cupoArray[currentSchool]/2),2)));

    // seg de la escuela antigua
    totalAluCol = aluxcol[newSchool];
    //cout << "mobio Alumnos antigua escuela "<< totalAluCol << " " << endl;
    aluVulCol = aluVulxCol[newSchool];
    aluNoVulCol =totalAluCol - aluVulCol;
    totalSesc+=round_n(fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln))));

    // costcupo escuela antigua
    totalcostCupo+=round_n(((double)totalAluCol*fabs(((double)cupoArray[newSchool]-totalAluCol)/pow(((double)cupoArray[newSchool]/2),2))));
    
    currentVars[0] = sumDist;
    currentVars[1] = totalSesc;
    currentVars[2] = totalcostCupo;

    //cout << fixed << setprecision(100) << endl;
    double var1 = (currentVars[0]/n_students);
    var1= (var1/max_dist);
    //cout << var1 << "\n";
    double var2 = (currentVars[1]/2.0);
    //cout << var2 << "\n";
    double var3 = (currentVars[2] /n_colegios);

    return ((double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3)));
}

double solutionNE3(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep)
    {

    double **mSolution = (double**)malloc(n_block * sizeof(double));
    int **colchange = (int**)malloc(n_block * sizeof(int));
    for (int x=0; x<n_block; x++){
        colchange[x] = (int*)malloc(n_thread * sizeof(int));
        mSolution[x] = (double*)malloc(n_thread * sizeof(double));
    }

    shuffle(shuffle_student,n_block,dist);
    shuffle(shuffle_colegios,n_thread,dist2);

    int *c_aluxcol=(int *)malloc(sizeof(int)*n_colegios);
    int *c_aluVulxCol=(int *)malloc(sizeof(int)*n_colegios);
    int *c_currentSolution =(int *)malloc(sizeof(int)*n_students);

    double minSolution = 0;

    int cAlu = 0;
    int cCol = 0;
    int *aluchange = (int*)malloc(n_block * sizeof(int));
    int bestAluchange;
    int bestColchange;

    for(int i =0; i < n_block; i++)
    {
        cCol = 0;
        aluchange[i] = shuffle_student[cAlu];
        for(int j=0; j < n_thread; j++)
        {
            // Si el colegio que selecciono es distinto al colegio que ya pertenece el estudiante
            if(shuffle_colegios[cCol] != currentSolution[shuffle_student[cAlu]]){ 
                colchange[i][j] = shuffle_colegios[cCol];
            }
            else{
                cCol++;
                colchange[i][j] = shuffle_colegios[cCol];
            }
            memcpy(c_aluxcol,aluxcol,sizeof(int)*n_colegios);
            memcpy(c_aluVulxCol,aluVulxCol,sizeof(int)*n_colegios);
            memcpy(c_currentSolution,currentSolution,sizeof(int)*n_students);
        
            //ELimina el estudiante de la escuela actual
            c_aluxcol[c_currentSolution[aluchange[i]]]-=1; ///
            c_aluVulxCol[c_currentSolution[aluchange[i]]]-=alumnosSep[aluchange[i]]; ///
            c_aluxcol[colchange[i][j]]+=1; ///
            c_aluVulxCol[colchange[i][j]]+=alumnosSep[aluchange[i]]; ///
            c_currentSolution[aluchange[i]] = colchange[i][j]; ///
            mSolution[i][j] = newSolution_v2(n_students,n_colegios,totalVuln,c_aluxcol,c_aluVulxCol,cupoArray,distMat,c_currentSolution,ptr_alpha);
            cCol++;
        }
        cAlu++;
    }
    minSolution= mSolution[0][0];
    bestAluchange = aluchange[0];
    bestColchange = colchange[0][0];

    for(int i =0; i < n_block; i++)
    {
        for(int j=0; j < n_thread; j++)
        {
        
            if(mSolution[i][j] <= minSolution){
                minSolution = mSolution[i][j];
                bestAluchange=aluchange[i];
                bestColchange = colchange[i][j];
            }
        }
    }


    //ELimina el estudiante de la escuela actual

    aluxcol[currentSolution[bestAluchange]]-=1; ///
    aluVulxCol[currentSolution[bestAluchange]]-=alumnosSep[bestAluchange]; ///
    aluxcol[bestColchange]+=1; ///
    aluVulxCol[bestColchange]+=alumnosSep[bestAluchange]; ///
    currentSolution[bestAluchange] = bestColchange; ///


    // For que busca al menor
    for(int i =0; i < n_block; i++)
    {
        free(colchange[i]);
        free(mSolution[i]);
    }
    free(c_aluxcol);
    free(c_aluVulxCol);
    free(c_currentSolution);
    free(aluchange);

    
    return minSolution;

        // Obtiene el costo Actual

    
}



double solutionNE4(int n_students,
        int n_colegios,
        int totalVuln,
        int *aluxcol,
        int *aluVulxCol,
        int *cupoArray,
        double **distMat, 
        int *currentSolution,
        double costCurrentSolution,
        const double *ptr_alpha,
        int *shuffle_student,
        int *shuffle_colegios,
        int *alumnosSep)
    {
    bool foundBetter = false;
    double** mSolution = (double**)malloc(n_block * sizeof(double));
    int** colchange = (int**)malloc(n_block * sizeof(int));
    for (int x=0; x<n_block; x++){
        colchange[x] = (int*)malloc(n_thread * sizeof(int));
        mSolution[x] = (double*)malloc(n_thread * sizeof(double));
    }

    shuffle(shuffle_student,n_block,dist);
    shuffle(shuffle_colegios,n_thread,dist2);

    int *c_aluxcol=(int *)malloc(sizeof(int)*n_colegios);
    int *c_aluVulxCol=(int *)malloc(sizeof(int)*n_colegios);
    int *c_currentSolution =(int *)malloc(sizeof(int)*n_students);

    double minSolution = 0;

    int cAlu = 0;
    int cCol = 0;
    int *aluchange = (int*)malloc(n_block * sizeof(int));
    int bestAluchange;
    int bestColchange;

    for(int i =0; i < n_block; i++)
    {
        cCol = 0;
        aluchange[i] = shuffle_student[cAlu];
        for(int j=0; j < n_thread; j++)
        {
            // Si el colegio que selecciono es distinto al colegio que ya pertenece el estudiante
            if(shuffle_colegios[cCol] != currentSolution[shuffle_student[cAlu]]){ 
                colchange[i][j] = shuffle_colegios[cCol];
            }
            else{
                cCol++;
                colchange[i][j] = shuffle_colegios[cCol];
            }
            memcpy(c_aluxcol,aluxcol,sizeof(int)*n_colegios);
            memcpy(c_aluVulxCol,aluVulxCol,sizeof(int)*n_colegios);
            memcpy(c_currentSolution,currentSolution,sizeof(int)*n_students);
        
            //ELimina el estudiante de la escuela actual
            c_aluxcol[c_currentSolution[aluchange[i]]]-=1; ///
            c_aluVulxCol[c_currentSolution[aluchange[i]]]-=alumnosSep[aluchange[i]]; ///
            c_aluxcol[colchange[i][j]]+=1; ///
            c_aluVulxCol[colchange[i][j]]+=alumnosSep[aluchange[i]]; ///
            c_currentSolution[aluchange[i]] = colchange[i][j]; ///
            mSolution[i][j] = newSolution_v2(n_students,n_colegios,totalVuln,c_aluxcol,c_aluVulxCol,cupoArray,distMat,c_currentSolution,ptr_alpha);
            cCol++;
            if(mSolution[i][j] < costCurrentSolution){
                foundBetter = true;
                minSolution = mSolution[i][j];
                bestAluchange = aluchange[i];
                bestColchange = colchange[i][j];
                break;
            }
        }
        cAlu++;
        if(foundBetter){
            break;
        }
    }
    if(!foundBetter){
        minSolution= mSolution[0][0];
        bestAluchange = aluchange[0];
        bestColchange = colchange[0][0];

        for(int i =0; i < n_block; i++)
        {
            for(int j=0; j < n_thread; j++)
            {
            
                if(mSolution[i][j] <= minSolution){
                    minSolution = mSolution[i][j];
                    bestAluchange=aluchange[i];
                    bestColchange = colchange[i][j];
                }
            }
        }

    }
    aluxcol[currentSolution[bestAluchange]]-=1; ///
    aluVulxCol[currentSolution[bestAluchange]]-=alumnosSep[bestAluchange]; ///
    aluxcol[bestColchange]+=1; ///
    aluVulxCol[bestColchange]+=alumnosSep[bestAluchange]; ///
    currentSolution[bestAluchange] = bestColchange; ///
    
    // For que busca al menor
    for(int i =0; i < n_block; i++)
    {
        free(colchange[i]);
        free(mSolution[i]);
    }
    free(c_aluxcol);
    free(c_aluVulxCol);
    free(c_currentSolution);
    free(aluchange);

    
    return minSolution;

        // Obtiene el costo Actual

    
}


