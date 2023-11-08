#include <Dataset.hpp>

Dataset::Dataset(std::string fileName_school, std::string fileName_students) {
    totalVuln = 0;
    ///////////////////////////////////////////////////
    /// Datos colegios
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura Info_colegio
    ///////////////////////////////////////////////////
    getDataSchool(fileName_school, colegios);
    ptr_colegios = colegios.data();
    n_colegios = colegios.size();
    ///////////////////////////////////////////////////
    /// Datos Alumnos
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura info_student
    ///////////////////////////////////////////////////
    getDataStudents(fileName_students, students, totalVuln);
    ptr_students = students.data();
    n_students = students.size();
}

void Dataset::getDataSchool(std::string fileName_school, std::vector<Info_colegio> &colegios) {
    string line_colegios;
    ifstream info_school(fileName_school); // concatenar
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

void Dataset::getDataStudents(std::string fileName_students, std::vector<Info_alu> &students, int &totalVuln) {
    string line_student;
    ifstream info_student(fileName_students); // concatenar
    int cx = 0;
    while (getline(info_student, line_student)) {
        stringstream linestream(line_student);
        string data;
        students.push_back(Info_alu());
        getline(linestream, data, ',');
        students[cx].rbd = std::stoi(data);
        getline(linestream, data, ',');
        students[cx].latitude = std::stod(data);
        getline(linestream, data, ',');
        students[cx].longitude = std::stod(data);
        getline(linestream, data, ',');
        students[cx].sep = std::stoi(data);
        if (students[cx].sep == 1)
            totalVuln++;
        cx++;
    }
    info_student.close();
}

/*
void Dataset::toCSV(std::string fileName) {
    try {
        std::ofstream pw(fileName);
        if (pw.is_open()) {
            for (int i = 0; i < X.size(); i++) {
                pw << X[i] << "," << Y[i] << "\n";
            }
            pw.close();
        }
    } catch (const std::exception &e) {
        std::cerr << e.what();
    }
}
*/
