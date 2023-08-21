#ifndef RECORD_MANAGER_H
#define RECORD_MANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>

// struct RecordParams
// {
//     std::string name_info;
//     std::string name_info_register;
//     std::string name_info_graphics;
// };

class RecordManager
{
private:
    std::ofstream info;
    std::ofstream infoRegister;
    std::ofstream infoGraphics;
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
    std::array<std::string, 3> path_names;
    void open_file(std::size_t);

public:
    RecordManager(const std::string &ruta_save_, const std::string &prefijo_save_);
    void initRecordInfo();
    void initRecordRegister();
    void initRecordGraphics();
    void SaveInfo();
    void SaveInfoRegister();
    void SaveGraphics();
    ~RecordManager();
};

#endif