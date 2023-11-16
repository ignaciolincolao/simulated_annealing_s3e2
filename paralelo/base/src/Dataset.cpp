#include <Dataset.hpp>

Dataset::Dataset(std::string fileName_school, std::string fileName_students, std::string fileName_parents) {
    totalVuln = 0;
    ///////////////////////////////////////////////////
    /// Datos colegios
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura Info_colegio
    ///////////////////////////////////////////////////
    getDataSchool(fileName_school);
    ///////////////////////////////////////////////////
    /// Datos Alumnos
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura info_student
    ///////////////////////////////////////////////////
    getDataStudents(fileName_students);
    getDataParents(fileName_parents);
}

Info_colegio *Dataset::get_colleges() { return colegios.data(); }
Info_alu *Dataset::get_students() { return students.data(); }
std::size_t Dataset::getn_colleges() { return colegios.size(); }
std::size_t Dataset::getn_students() { return students.size(); }
std::size_t Dataset::get_totalVuln() { return totalVuln; }

void Dataset::getDataSchool(std::string fileName_school) {
    std::string line_colegios;
    std::ifstream info_school(fileName_school); // concatenar
    int cx = 0;
    while (getline(info_school, line_colegios)) {
        std::stringstream linestream(line_colegios);
        std::string data;
        colegios.push_back(Info_colegio());
        getline(linestream, data, ',');
        colegios[cx].rbd = std::stoi(data);
        getline(linestream, data, ',');
        colegios[cx].latitude = std::stod(data);
        getline(linestream, data, ',');
        colegios[cx].longitude = std::stod(data);
        getline(linestream, data, ',');
        colegios[cx].num_alu = std::stoi(data);
        getline(linestream, data, ',');
        colegios[cx].prioritario = std::stoi(data);
        cx++;
    }
    info_school.close();
}

void Dataset::getDataStudents(std::string fileName_students) {
    std::string line_student;
    std::ifstream info_student(fileName_students); // concatenar
    int cx = 0;
    while (getline(info_student, line_student)) {
        std::stringstream linestream(line_student);
        std::string data;
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

void Dataset::getDataParents(std::string fileName_parents) {
    std::string line_parent;
    std::ifstream info_parent(fileName_parents);
    std::size_t cx = 0;

    while (getline(info_parent, line_parent)) {
        std::stringstream line_stream(line_parent);
        std::string data;
        getline(line_stream, data, ',');
        students[cx].choices[0] = std::stoi(data);
        getline(line_stream, data, ',');
        students[cx].choices[1] = std::stoi(data);
        getline(line_stream, data, ',');
        students[cx].choices[2] = std::stoi(data);
        getline(line_stream, data, ',');
        students[cx].choices[3] = std::stoi(data);
        getline(line_stream, data, ',');
        students[cx].choices[4] = std::stoi(data);
        cx++;
    }
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
