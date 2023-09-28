#ifndef DETAILED_RUN_H
#define DETAILED_RUN_H

#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include "rescaler.h"
#include "task.h"
#include "ukf.h"

class DetailedRun{
public:
    DetailedRun(const std::string base_path, const std::string &blood_path, const std::string &heart_path, Eigen::MatrixXd &Teta, const int &age, const double &HR) {
        this -> Teta = Teta;
        this -> HR = HR;
        this -> base_path = base_path;
        this -> blood_path = blood_path;
        this -> heart_path = heart_path;
        rescaler.init(this -> blood_path);
    }

    void run_task() {
        Task task;
        task.set_base_path(base_path);
        Ukf::set_up_task(task, Teta, blood_path, heart_path, base_path, rescaler);

        task.run_full(true);
        task.dump_data();

        X(0) = task.P_sys;
        X(1) = task.P_dis;
        X(2) = task.get_strokeVolume();

    }


private:
    Eigen::Vector3d X;
    Eigen::MatrixXd Teta;
    double HR;
    std::string base_path;
    std::string blood_path;
    std::string heart_path;
    Rescaler rescaler;
};

#endif // DETAILED_RUN_H
