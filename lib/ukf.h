#ifndef UKF_H
#define UKF_H

#include <Eigen/Dense>
#include "task.h"
#include <array>
#include <string>
#include "rescaler.h"

class Ukf
{
public:
    Ukf(const double& P_sys, const double& P_dis, const double& SV, const std::string blood_config, const std::string heart_config, const std::string base_path, bool restore_from_backup = false);
    void execute_pipeline();
    void set_up_task(Task &task, const Eigen::Vector<double, 8> &Teta);
    static void set_up_task(Task &task, const Eigen::Vector<double, 8> &Teta, const std::string blood_config, const std::string heart_config, const std::string base_path, Rescaler &rescaler);
    static void run_task(Task &task, Eigen::Vector3d& X);

private:
    std::string base_path;
    std::string log_path;
    std::ofstream logger;
    std::ofstream backup_stream_U;
    std::ifstream backup_istream_U;
    std::ofstream backup_stream_X_n;
    std::ifstream backup_istream_X_n;
    std::ofstream backup_stream_Teta_n;
    std::ifstream backup_istream_Teta_n;
    std::ofstream backup_stream_L_x;
    std::ifstream backup_istream_L_x;
    std::ofstream backup_stream_L_teta;
    std::ifstream backup_istream_L_teta;
    std::string path_to_backup;
    bool restore_from_backup;
    std::string blood_config;
    std::string heart_config;
    const static int p = 8; // number of params
    const static int k = 3; // number of target variables
    const double sigma_param_sqr = 1e-2;
    const double sigma_obs_sqr = 1.0;
    const double alpha = 1.0 / (p + 1); // weight of sigma-point
    Eigen::DiagonalMatrix<double, k> P_z;
    Eigen::DiagonalMatrix<double, p> P_teta;
    Eigen::Matrix<double, p, p + 1> I; // matrix of sigma points
    Eigen::DiagonalMatrix<double, p + 1> D; // matrix of sigma points weights
    Eigen::Matrix<double, k, p + 1> Z; // matrix of measurements in sigma-points
    Eigen::Vector3d X_n;
    Eigen::Vector<double, p> Teta_n;
    Eigen::Vector<double, p> Teta_initial;
    Eigen::Matrix<double, k, p> L_x;
    Eigen::Matrix<double, p, p> L_teta;
    Eigen::Matrix<double, k, p> L_gamma;
    Eigen::Matrix<double, p, p> U;
    std::array<Eigen::Vector3d, p + 1> arr_X;
    std::array<Eigen::Vector3d, p + 1> arr_X_next;
    std::array<Eigen::Vector<double, p>, p + 1> arr_Teta;
    std::array<Eigen::Vector<double, p>, p + 1> arr_Teta_next;
    Eigen::Vector3d X_next_min;
    Eigen::Vector<double, p> Teta_next_min;
    Eigen::Matrix<double, k, p + 1> Gamma;
    Eigen::Matrix<double, k, p + 1> H; // X_n+1_minus in matrix form
    Eigen::Matrix<double, p, p + 1> H_teta; // Teta_n+1_minus in matrix form
    Rescaler rescaler;

    void define_initial_conditions();
    void define_matrix_sigma_points();
    void prediction();
    void update_covariances();
    void correction();
    void correct_Teta(Eigen::Vector<double, p>& Teta);
    double get_error(Eigen::Vector3d& X_res);
    void collect_backup_data();
};

#endif // UKF_H
