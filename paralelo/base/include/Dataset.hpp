#ifndef DATASET_H
#define DATASET_H

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

///////////////////////////////////////////////////
/// Estructura de datos de los colegios.
///////////////////////////////////////////////////
struct Info_colegio {
    double latitude = 0.0;
    double longitude = 0.0;
    int num_alu = 0;
    int rbd = 0;
    int prioritario = 0;
};
///////////////////////////////////////////////////
/// Estructura de alumnos
///////////////////////////////////////////////////
struct Info_alu {
    int rbd = 0;
    int sep = 0;
    double latitude = 0.0;
    double longitude = 0.0;
    std::array<uint8_t, 5> choices;
};

class Dataset {

private:
    std::size_t totalVuln;
    std::vector<Info_colegio> colegios;
    std::vector<Info_alu> students;
    void getDataSchool(std::string fileName_school);
    void getDataStudents(std::string fileName_students);
    void getDataParents(std::string fileName_parents);

public:
    Dataset(std::string fileName_school, std::string fileName_students, std::string fileName_parents);
    Info_alu *get_students();
    Info_colegio *get_colleges();
    std::size_t getn_students();
    std::size_t getn_colleges();
    std::size_t get_totalVuln();

    // void toCSV(std::string fileName);
};

#endif
