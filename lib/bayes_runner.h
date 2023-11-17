#ifndef BAYES_RUNNER_H
#define BAYES_RUNNER_H

#include <string>
#include <Eigen/Dense>
#include "rescaler.h"

class BayesRunner
{
public:
    BayesRunner(std::string base_path, const std::string blood_config, const std::string heart_config);
    void execute_pipeline();

private:
    std::string base_path;
    std::string blood_config;
    std::string heart_config;
    std::string read_path;
    std::string write_path;
    time_t last_modification;
    Eigen::Vector<double, 10> Teta;
    Eigen::Vector3d X;

    Rescaler rescaler;

    void wait_point();
    void read_point();
    void run_task();
    void write_results();
    void write_error();
    void write_initial_point();
};

#endif // BAYES_RUNNER_H
