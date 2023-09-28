#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <fstream>
#include <Eigen/Dense>

double string_to_double(std::string& s);
void write_csv_matrix(std::ofstream &stream, Eigen::Ref <Eigen::MatrixXd> &m);
void read_csv_matrix(std::ifstream &stream, Eigen::MatrixXd &m, const int &str_nums, const int &col_nums);

void set_hr(const double &HR, const std::string &config_path);

#endif // MATRIX_UTILS_H
