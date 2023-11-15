#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <sstream>
#include <json.hpp>
#include <vector>
using json = nlohmann::json;

double string_to_double(std::string& s) {
    std::stringstream ss;
        double d;
        ss << s;
        ss >> d;
        return d;
}

void write_csv_matrix(std::ofstream &stream, Eigen::Ref <Eigen::MatrixXd> &m) {
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            if (j != m.cols() - 1)
                stream << m(i, j) << ",";
            else
                stream << m(i, j) << std::endl;
        }
    }
}

void read_csv_matrix(std::ifstream &stream, Eigen::MatrixXd &m, const int &str_nums, const int &col_nums) {
        std::string line;
        std::string part;
        std::string sep = ",";
        m.resize(str_nums, col_nums);
        int i = 0;
        while(getline(stream, line)) {
                int finish = 0;
                for (int j = 0; j < col_nums - 1; ++j) {
                        part = line.substr(0, line.find(sep));
                        finish = line.find(sep) + 1;
                        line.erase(0, finish);
                        m(i, j) = string_to_double(part);
                }
                m(i, col_nums - 1) = string_to_double(line);
                ++i;
        }
}

void csv_to_array_of_strings(std::ifstream &stream, std::vector<std::vector<std::string>> &arr, const int &str_nums, const int &col_nums) {
    std::string line;
    std::string part;
    std::string sep = ",";

    arr.resize(str_nums);

    int i = 0;
    while(getline(stream, line)) {
        int finish = 0;
        for (int j = 0; j < col_nums - 1; ++j) {
            part = line.substr(0, line.find(sep));
            finish = line.find(sep) + 1;
            line.erase(0, finish);
            arr[i].push_back(part);
        }
        arr[i].push_back(line);
        ++i;
    }
}

void set_hr(const double &HR, const std::string &config_path) {
    json config;
    std::ifstream config_stream;
    config_stream.open(config_path);
    config_stream >> config;
    config_stream.close();
    config["MetaVertices"]["Heart_AdValves_healthy"]["Heart_rate"] = HR;
    std::ofstream outfile;
    outfile.open(config_path);
    outfile << config;
    outfile.close();
}

