#include "bayes_runner.h"
#include <filesystem>
#include <iostream>
#include <chrono>
#include <thread>
#include <fstream>
#include "matrix_utils.h"
#include "ukf.h"
#include "task.h"
#include <string>

BayesRunner::BayesRunner(std::string base_path, const std::string blood_config, const std::string heart_config) {
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;

    this -> base_path.assign(base_path);
    this -> blood_config.assign(blood_config);
    this -> heart_config.assign(heart_config);
    this -> read_path.assign(base_path + "/data/share_folder/next_point.csv");
    this -> write_path.assign(base_path + "/data/share_folder/next_target.csv");
    rescaler.init(this -> blood_config);

    if (!std::filesystem::exists(read_path)) {
        std::ofstream read_file(read_path);
        read_file.close();
    }
    last_modification = std::filesystem::last_write_time(read_path).time_since_epoch().count();

    // the way of first run like in ukf now

    Eigen::MatrixXd Teta_dynamic;
    std::ifstream param_file(base_path + "/data/back/base.csv");
    read_csv_matrix(param_file, Teta_dynamic, 12, 1);
    Teta = Teta_dynamic;
    run_task();

    /*
    Task task;
    task.set_path_to_brachial_data(base_path);
    task.set_up(this -> blood_config, this -> heart_config);
    Ukf::run_task(task, X);

    Teta(0) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_R4();
    Teta(1) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_R1();
    Teta(2) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_PveinPressure();
    Teta(3) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_Kp();
    Teta(4) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_Kf();    
    Teta(5) = rescaler.find_total_compliance();
    Teta(6) = rescaler.find_total_resistance();
    Teta(7) = rescaler.get_p_out();
    Teta(8) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_I4();
    Teta(9) = dynamic_cast<Heart_AdValves*>(task.true_0d_heart) -> get_I1();
    Teta(10) = 1.0;
    Teta(11) = 1.0;
    */

    write_initial_point();
}

void BayesRunner::write_initial_point() {
    std::ofstream ostream;
    ostream.open(base_path + "/data/share_folder/params_0.csv");
    Eigen::Ref <Eigen::MatrixXd> Teta_ref(Teta);
    write_csv_matrix(ostream, Teta_ref);
    ostream.close();


    ostream.open(base_path + "/data/share_folder/targets_0.csv");
    Eigen::Ref <Eigen::MatrixXd> X_ref(X);
    write_csv_matrix(ostream, X_ref);
    ostream.close();
}

void BayesRunner::wait_point() {
    const std::filesystem::path path = read_path;
    bool not_changed = true;
    const std::chrono::seconds pause(10);
    time_t current_modification;

    while (not_changed) {
        current_modification = std::filesystem::last_write_time(path).time_since_epoch().count();
        std::cout << current_modification << std::endl;
        std::cout << last_modification << std::endl;
        if (current_modification > last_modification) {
            last_modification = current_modification;
            not_changed = false;
        }
        else {
            std::cout << "sleep\n";
            std::this_thread::sleep_for(pause);
        }
    }
}

void BayesRunner::read_point() {
    std::ifstream istream;
    Eigen::MatrixXd Teta_dyn;

    istream.open(read_path);
    read_csv_matrix(istream, Teta_dyn, 12, 1);
    Teta = Teta_dyn;
}

void BayesRunner::run_task() {
    Task task;
    Ukf::set_up_task(task, Teta, blood_config, heart_config, base_path, rescaler);
    Ukf::run_task(task, X);
    std::cout << "res : \n";
    std::cout << X << std::endl;
}

void BayesRunner::write_results() {
    std::ofstream ostream;
    ostream.open(write_path);
    Eigen::Ref <Eigen::MatrixXd> X_ref(X);
    write_csv_matrix(ostream, X_ref);
}

void BayesRunner::write_error() {
    std::ofstream ostream;
    ostream.open(write_path);
    ostream << "err\n";
  }

void BayesRunner::execute_pipeline() {
    while (true) {
        wait_point();
        read_point();
        try {
            run_task();
            write_results();
        }
        catch (char const *err) {
            std::cout << "from catch block\n";
            std::cout << err << std::endl;
            write_error();
        }
    }
}
