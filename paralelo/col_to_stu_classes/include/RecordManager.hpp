#ifndef RECORD_MANAGER_H
#define RECORD_MANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <utils/SAParameters.hpp>

struct RecordParams
{
    std::string prefijo_save;
    std::string ruta_save;
};

class RecordManager
{
private:
    std::ofstream info;
    std::ofstream infoRegister;
    std::ofstream infoGraphics;
    std::array<std::string, 3> path_names;
    RecordParams &rMgrParams;
    SimulatedParams &saParams;

private:
    void open_file(std::size_t n_file, std::ofstream &file);

public:
    std::vector<double> vector_costCurrentSolution;
    std::vector<double> vector_meanDist;
    std::vector<double> vector_segregation;
    std::vector<double> vector_costoCupo;
    std::vector<double> vector_temp;
    std::vector<int> vector_count;
    std::vector<double> vector_historyCostSolution;
    std::vector<double> vector_historyTemp;
    std::vector<double> vector_historymeanDist;
    std::vector<double> vector_historymeanDistNorm;
    std::vector<double> vector_historySegregation;
    std::vector<double> vector_historycostoCupo;
    std::vector<bool> vector_historyAcceptSolution;
    std::vector<int> vector_historyAsign;
    std::vector<std::tuple<int, int>> vector_historyMove;

public:
    RecordManager(SimulatedParams &saParams_, RecordParams &params_);

    RecordParams &getRmgrParams() { return rMgrParams; }

    void openRecordInfo();
    void openRecordRegister();
    void openRecordGraphics();

    void closeRecordInfo();
    void closeRecordRegister();
    void closeRecordGraphics();

    void SaveInfoInit(double costBestSolution,
                      double meanDist,
                      double S,
                      double costCupo);

    void SaveInfoFinish(double costPreviousSolution,
                        double costBestSolution,
                        double costCurrentSolution,
                        double time_taken,
                        double meanDist,
                        double S,
                        double costCupo);

    void SaveInfoRegisterInit();

    void SaveInfoRegisterFinish();

    void SaveGraphicsInit(double meanDist,
                          double S,
                          double costCupo,
                          double costCurrentSolution);

    void SaveGraphicsFinish();

    ~RecordManager();
};

#endif