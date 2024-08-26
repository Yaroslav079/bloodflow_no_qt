#include "ukf.h"
#include <cmath>
#include <random>
#include "task.h"
#include "matrix_utils.h"
#include "rescaler.h"
#include <string>
#include <ctime>
#include <filesystem>

std::string getEpochTimeString() {
    std::time_t t = std::time(nullptr);   // Get current time
    std::stringstream ss;
    ss << std::time(nullptr);
    return ss.str();
}

Ukf::Ukf(const double& P_sys, const double& P_dis, const double& SV, std::string blood_config, std::string heart_config, const std::string base_path, bool restore_from_backup) {
    P_z.diagonal() << sigma_obs_sqr, sigma_obs_sqr, sigma_obs_sqr;
    for (int i = 0; i < p + 1; ++i) {
        Z(0, i) = P_sys;
        Z(1, i) = P_dis;
        Z(2, i) = SV;
    }

    this -> restore_from_backup = restore_from_backup;
    this -> blood_config.assign(blood_config);
    this -> heart_config.assign(heart_config);

    rescaler.init(blood_config);

    log_path.append(base_path).append("/data/log/");

    if (this -> restore_from_backup == false) {
        path_to_backup.append(base_path).append("/data/back/").append(getEpochTimeString());
        std::filesystem::create_directory(path_to_backup);
        path_to_backup.append("/");
    }
    else {
        std::string backup_folder_name;
        std::cout << "Name of foler with backup files : ";
        std::cin >> backup_folder_name;
        path_to_backup.append(base_path).append("/data/back/").append(backup_folder_name);
        path_to_backup.append("/");
    }

    this -> base_path = base_path;
};

void Ukf::set_up_task(Task &task, Eigen::Vector<double, 12> &Teta) {
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;

    std::cout << "params to run : \n";
    std::cout << Teta << std::endl;
    task.set_path_to_brachial_data(this -> base_path);
    rescaler.create_wk_distribution(Teta(6), Teta(5), Teta(7), Teta(10), Teta(11));
    task.set_up(this -> blood_config, this -> heart_config);
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_R4(Teta(0));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_R1(Teta(1));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_PveinPressure(Teta(2));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_Kp(Teta(3));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_Kf(Teta(4));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_I4(Teta(8));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_I1(Teta(9));
}

// unused now, for final test only
void Ukf::set_up_task(Task &task, const Eigen::Vector<double, 12> &Teta, const std::string blood_config, const std::string heart_config, const std::string base_path, Rescaler &rescaler) {
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;

    std::cout << "params to run : \n";
    std::cout << Teta << std::endl;
    task.set_path_to_brachial_data(base_path);
    rescaler.create_wk_distribution(Teta(6), Teta(5), Teta(7), Teta(10), Teta(11));
    task.set_up(blood_config, heart_config);
    
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_R4(Teta(0));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_R1(Teta(1));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_PveinPressure(Teta(2));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_Kp(Teta(3));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_Kf(Teta(4));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_I4(Teta(8));
    dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> set_I1(Teta(9));
}

void Ukf::run_task(Task &task, Eigen::Vector3d& X) {
    try {
        task.run_full();
    }
    catch(...) {
        std::cout << "Error was catched in Ukf::run_task()" << std::endl;
        throw("Calculation error");
    }
    task.dump_data();

    X(0) = task.P_sys;
    X(1) = task.P_dis;
    X(2) = task.get_strokeVolume();
}

void Ukf::define_matrix_sigma_points() {
    const double alpha_sigma = p / (p + 1.0);
    logger << alpha_sigma << std::endl;
    I(0, 0) = - 1.0 / pow(2 * alpha_sigma, 0.5);
    I(0, 1) = 1.0 / pow(2 * alpha_sigma, 0.5);
    for (int i = 2; i < p + 1; ++i) {
        I(0, i) = 0.0;
    }
    for (int d = 1; d < p; ++d) {
        for (int i = 0; i < p + 1; ++i) {
            if (i <= d) {
                I(d, i) = 1.0 / pow(alpha_sigma * (d + 1) * (d + 2), 0.5);
            }
            else if (i == d + 1) {
                I(d, i) = - (d + 1) / pow(alpha_sigma * (d + 1) * (d + 2), 0.5);
            }
            else {
                I(d, i) = 0.0;
            }
        }
    }
    I = I * pow(p, 0.5);

    Eigen::RowVector<double, p + 1> d_vector;
    for (int i = 0; i < p + 1; ++i)
        d_vector[i] = alpha;
    D = d_vector.asDiagonal();

    logger << "I" << std::endl;
    logger << I << std::endl;

    logger << "D" << std::endl;
    logger << D.toDenseMatrix() << std::endl;
};

void Ukf::define_initial_conditions() {
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;

    /*
    old-fashioned first definitions based on inner model params

    Task task;
    task.set_path_to_brachial_data(this -> base_path);
    task.set_up(this -> blood_config, this -> heart_config);
    run_task(task, X_n);

    Teta_n(0) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_R4();
    Teta_n(1) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_R1();
    Teta_n(2) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_PveinPressure();
    Teta_n(5) = rescaler.find_total_compliance();
    Teta_n(6) = rescaler.find_total_resistance();
    Teta_n(3) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_Kp();
    Teta_n(4) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_Kf();
    */

    // new way of defining based on closer base
    Eigen::MatrixXd Teta_dynamic;
    std::ifstream param_file(base_path + "/data/back/base.csv");
    read_csv_matrix(param_file, Teta_dynamic, 12, 1);
    Teta_n = Teta_dynamic;
    param_file.close();
    Task task;
    set_up_task(task, Teta_n);
    run_task(task, X_n);

    Teta_initial = Teta_n;

    Eigen::RowVector<double, p> p_teta_vector;
    for (int i = 0; i < p; ++i)
        p_teta_vector[i] = sigma_param_sqr * pow(Teta_initial(i), 2.0);
    P_teta = p_teta_vector.asDiagonal();

    if (restore_from_backup) {
        backup_istream_Teta_n.open(path_to_backup + "Teta_n.csv");
        Eigen::MatrixXd Teta_n_dyn;
        read_csv_matrix(backup_istream_Teta_n, Teta_n_dyn, 11, 1);
        Teta_n = Teta_n_dyn;
        backup_istream_Teta_n.close();

        backup_istream_X_n.open(path_to_backup + "X_n.csv");
        Eigen::MatrixXd X_n_dyn;
        read_csv_matrix(backup_istream_X_n, X_n_dyn, 3, 1);
        X_n = X_n_dyn;
        backup_istream_X_n.close();

        logger << "restored \n";
        logger << X_n << std::endl;
        logger << Teta_n << std::endl;
    }

    logger << "P_teta" << std::endl;
    logger << P_teta.toDenseMatrix() << std::endl;
    logger << "P_z" << std::endl;
    logger << P_z.toDenseMatrix() << std::endl;

    logger << X_n << std::endl;
    logger << "norm : " << get_error(X_n) << std::endl;
    logger << Teta_n << std::endl;

    if (restore_from_backup) {
        logger << "restored\n";
        backup_istream_L_x.open(path_to_backup + "L_x.csv");
        Eigen::MatrixXd L_x_dyn;
        read_csv_matrix(backup_istream_L_x, L_x_dyn, 3, 12);
        L_x = L_x_dyn;
        backup_istream_L_x.close();

        backup_istream_L_teta.open(path_to_backup + "L_teta.csv");
        Eigen::MatrixXd L_teta_dyn;
        read_csv_matrix(backup_istream_L_teta, L_teta_dyn, 12, 12);
        L_teta = L_teta_dyn;
        backup_istream_L_teta.close();

        backup_istream_U.open(path_to_backup + "U.csv");
        Eigen::MatrixXd U_dyn;
        read_csv_matrix(backup_istream_U, U_dyn, 12, 12);
        U = U_dyn;
        backup_istream_U.close();
    }
    else {
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                if (i != j)
                    L_teta(i, j) = 0.0;
                else
                    L_teta(i, j) = 1.0;
            }
        }
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < p; ++j) {
                L_x(i, j) = 0.0;
            }
        }
        U = P_teta.inverse();
    }

    logger << "L_x : \n" << L_x << std::endl;
    logger << "L_teta : \n" << L_teta << std::endl;
    logger << "U : \n" << U << std::endl;
};

void Ukf::prediction() {
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;

    Eigen::Matrix<double, p, p> Cn( U.inverse().llt().matrixL() ); // Cholesky decomposition
    logger << X_n << std::endl;
    logger << Teta_n << std::endl;
    for (int i = 0; i < p + 1; ++i) {
        arr_X[i] = X_n + L_x * Cn.transpose() * I.col(i);
        arr_Teta[i] = Teta_n + L_teta * Cn.transpose() * I.col(i);
        correct_Teta(arr_Teta[i]);
    }

    for (int i = 0; i < p + 1; ++i) {
        Task task;
        std::cout << i << std::endl;
        // set_up_task(task, arr_Teta[i]);
        // logger << i << " :\n"<< arr_Teta[i] << std::endl;
        // run_task(task, arr_X_next[i]);
        try {
            set_up_task(task, arr_Teta[i]);
            logger << i << " :\n"<< arr_Teta[i] << std::endl;
            run_task(task, arr_X_next[i]);
        }
        catch(...) {
            std::cout << "Error was catched in Ukf::prediction()" << std::endl;
            throw("Calculation error");
        }
        logger << i << "_res :\n"<< arr_X_next[i] << std::endl;
        arr_Teta_next[i] = arr_Teta[i];
    }

    X_next_min.setZero();
    Teta_next_min.setZero();

    for (int i = 0; i < p + 1; ++i) {
        X_next_min += alpha * arr_X_next[i];
        Teta_next_min += alpha * arr_Teta_next[i];
    }

    logger << "X_next_min : " << X_next_min << std::endl;
    logger << "Teta_next_min : " << Teta_next_min << std::endl;

    for (int i = 0; i < p + 1; ++i) {
        for (int j = 0; j < k; ++j) {
            H(j, i) = arr_X_next[i](j);
        }
    }

    for (int i = 0; i < p + 1; ++i) {
        for (int j = 0; j < p; ++j) {
            H_teta(j, i) = arr_Teta_next[i](j);
        }
    }

    logger << "H : " << H << std::endl;
    logger << "H_teta : " << H_teta << std::endl;

    Gamma = Z - H;

    logger << "Gamma : " << Gamma << std::endl;
};

void Ukf::update_covariances(){
    L_x = H * D * I.transpose();
    L_teta = H_teta * D * I.transpose();
    L_gamma = Gamma * D * I.transpose();
    U = I * D * I.transpose()  + L_gamma.transpose() * P_z.inverse()  * L_gamma;

    logger << "L_x : \n" << L_x << std::endl;
    logger << "L_teta : \n" << L_teta << std::endl;
    logger << "L_gamma : \n" << L_gamma << std::endl;
    logger << "U : \n" << U << std::endl;
};

void Ukf::correction(){
    Eigen::Vector3d Gamma_sum;
    Gamma_sum.setZero();
    for (int i = 0; i < p + 1; ++i) {
        Gamma_sum += alpha * Gamma.col(i);
    }
    X_n = X_next_min - L_x * U.inverse() * L_gamma.transpose() * P_z.inverse() * Gamma_sum;
    Teta_n = Teta_next_min - L_teta * U.inverse() * L_gamma.transpose() * P_z.inverse() * Gamma_sum;

    correct_Teta(Teta_n);

    logger << "X_n : " << X_n << std::endl;
    logger << "Teta_n : " << Teta_n << std::endl;
};

void Ukf::execute_pipeline(){
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;

    std::string final_log_path = log_path + getEpochTimeString() + ".txt";
    std::cout << "path to logs in this case : " << final_log_path << std::endl;
    logger.open(final_log_path);
    logger << "to find : ";
    logger << std::to_string(Z(0, 0)) + "," + std::to_string(Z(1, 0)) + "," + std::to_string(Z(2, 0)) << std::endl;
    logger << "for " + blood_config << std::endl;
    logger << "for " + heart_config << std::endl << std::endl;
    Eigen::Vector<double, p> Teta_before;
    Eigen::Vector<double, p> Teta_after;
    Eigen::Vector3d X_res;
    define_matrix_sigma_points();
    define_initial_conditions();
    double norm = get_error(X_n);
    int iter_num = 0;
    while (norm > 0.05) {
        Teta_before = Teta_n;
        try {
            prediction();
        }
        catch (...) {
            std::cout << "Error was catched in Ukf::execute_pipeline()" << std::endl;
            logger << "Error occured during prediction" << std::endl;
            logger << "Fixing results, go for next case" << std::endl;
            break;
        }
        update_covariances();
        correction();
        collect_backup_data();
        Teta_after = Teta_n;
        /*
        Task task;
        set_up_task(task, Teta_after);
        run_task(task, X_res);
        norm = get_error(X_res);
        */
        norm = get_error(X_n);
        logger << "iter_num : " << ++iter_num << std::endl;
        logger << "norm : " << norm << std::endl;
        if (norm < 0.05) {
            Task task;
            set_up_task(task, Teta_after);
            run_task(task, X_res);
            norm = get_error(X_res);
            logger << "norm_res : " << norm << std::endl;
            logger << "X_res : \n" << X_res << std::endl;
        }
        // logger << "X_res : \n" << X_res << std::endl;
        logger << "-------------------------------------------------------------------" << std::endl;
    }
    logger << "parameters were selected for (" << Z(0, 0) << "," << Z(1, 0) << "," << Z(2, 0)  << ")" << std::endl;
    logger << "Teta :\n" << Teta_n << std::endl;
    logger << "max_residual : " << norm << std::endl;
    logger << "iter_num : " << iter_num << std::endl;
    logger.close();
};

void Ukf::correct_Teta(Eigen::Vector<double, p>& Teta){
  for (int i = 0; i < p; ++i) {
      if (Teta(i) < 1e-15) {
          Teta(i) = Teta_initial(i) * pow(2.0, Teta(i));
      }
  }
};

void Ukf::collect_backup_data() {
    backup_stream_U.open(path_to_backup + "U.csv");
    Eigen::Ref <Eigen::MatrixXd> U_dyn(U);
    write_csv_matrix(backup_stream_U, U_dyn);
    backup_stream_U.close();

    backup_stream_X_n.open(path_to_backup + "X_n.csv");
    Eigen::Ref <Eigen::MatrixXd> X_n_dyn(X_n);
    write_csv_matrix(backup_stream_X_n, X_n_dyn);
    backup_stream_X_n.close();

    backup_stream_Teta_n.open(path_to_backup + "Teta_n.csv");
    Eigen::Ref <Eigen::MatrixXd> Teta_n_dyn(Teta_n);
    write_csv_matrix(backup_stream_Teta_n, Teta_n_dyn);
    backup_stream_Teta_n.close();

    backup_stream_L_x.open(path_to_backup + "L_x.csv");
    Eigen::Ref <Eigen::MatrixXd> L_x_dyn(L_x);
    write_csv_matrix(backup_stream_L_x, L_x_dyn);
    backup_stream_L_x.close();

    backup_stream_L_teta.open(path_to_backup + "L_teta.csv");
    Eigen::Ref <Eigen::MatrixXd> L_teta_dyn(L_teta);
    write_csv_matrix(backup_stream_L_teta, L_teta_dyn);
    backup_stream_L_teta.close();

    std::cout << "backup done; path_to_backup : " << path_to_backup << std::endl;
    logger << "backup done; path_to_backup : " << path_to_backup << std::endl;
}

double Ukf::get_error(Eigen::Vector3d& X_res) {
    double Psys = Z(0, 0);
    double Pdis = Z(1, 0);
    double SV = Z(2, 0);
    double eps = 1.0;
    // eps = fabs(X_res(0) - Psys) / Psys + fabs(X_res(1) - Pdis) / Pdis + fabs(X_res(2) - SV) / SV;

    if ((fabs(X_res(0) - Psys) / Psys) > (fabs(X_res(1) - Pdis) / Pdis))
        eps = fabs(X_res(0) - Psys) / Psys;
    else
        eps = fabs(X_res(1) - Pdis) / Pdis;
    if (eps < fabs(X_res(2) - SV) / SV)
        eps = fabs(X_res(2) - SV) / SV;

    return eps;
};

