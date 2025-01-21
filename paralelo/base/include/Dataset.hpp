#ifndef DATASET_H
#define DATASET_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <array>

using std::string;
using std::stof;
using std::stoi;
using std::stod;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::random_device;
using std::stringstream;
using std::getline;
using std::fixed;
using std::array;

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
    array<uint8_t, 5> choices; //el oscar tiene 5 preferencias
};

class Dataset {
    public:
        Dataset(std::string fileName_school, std::string fileName_students, std::string fileName_parents);
        std::vector<Info_colegio> colegios;
        std::vector<Info_alu> students;
        Info_colegio *ptr_colegios;
        Info_alu *ptr_students;
        int totalVuln;
        int n_colegios;
        int n_students;
        void getDataSchool(std::string fileName_school, std::vector<Info_colegio> &colegios);
        void getDataStudents(std::string fileName_students, std::vector<Info_alu> &students, int &totalVuln);
        void getDataParents(std::string fileName_parents, std::vector<Info_alu>& students);
        void printParentChoices() {
            for (size_t i = 0; i < students.size(); ++i) {
                std::cout << "Student " << i + 1 << " choices: ";
                for (size_t j = 0; j < students[i].choices.size(); ++j) {
                    std::cout << static_cast<int>(students[i].choices[j]) << " ";
                }
                std::cout << std::endl;
            }
        }

        //void toCSV(std::string fileName);
};


#endif
