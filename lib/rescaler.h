#ifndef RESCALER_H
#define RESCALER_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <json.hpp>
using json = nlohmann::json;

class Rescaler {
private:
    const double p = 2.27;
    const double resistence_rel = 0.25; // R1 / R2 = 0.25 => R1 + R2 = R2 * 1.25

    std::map<std::string, double> diameter_info;
    double diameter_sum; // (d_i / d_1)^p
    std::string id_1;
    double d_1;

    std::string config_path;
    json config;

    double previous_p0_coeff = 1.0;
    double previous_pwv_coeff = 1.0;

    void write_json(const std::string &path, json json) {
        std::ofstream outfile;
        outfile.open(path);
        outfile << json;
        outfile.close();
    };

    std::pair<double, double> split_resistance(const double &R) {
        double R_2 = R / (1.0 + resistence_rel);
        double R_1 = R_2 * resistence_rel;
        return {R_1, R_2};
    }

    void set_p_out(const double &p_out) {
        for (auto mv: config["MetaVertices"].items()) {
            if (mv.value()["Type"].template get<std::string>() == "Windkessel_vertex") {
                mv.value()["P_out"] = p_out;
            }
        }
    }

    void set_p0(const double &coeff_p0) {
        double p0_pr;
        for (auto mv: config["Edges"].items()) {
            if (mv.value()["Type"].template get<std::string>() == "Classic") {
                p0_pr = mv.value()["p0"].template get<double>();
                mv.value()["p0"] = coeff_p0 * p0_pr / previous_p0_coeff;
            }
        }
        previous_p0_coeff = coeff_p0;
    }

    void set_pwv(const double &coeff_pwv) {
        double c_pr;
        for (auto mv: config["Edges"].items()) {
            if (mv.value()["Type"].template get<std::string>() == "Classic") {
                c_pr = mv.value()["C"].template get<double>();
                mv.value()["C"] = coeff_pwv * c_pr / previous_pwv_coeff;
            }
        }
        previous_pwv_coeff = coeff_pwv;
    }

public:
    Rescaler() {}

    void init(const std::string &blood_config_path) {
        config_path = blood_config_path;
        std::ifstream config_stream;
        config_stream.open(config_path);
        config_stream >> config;
        config_stream.close();

        for (auto mv: config["MetaVertices"].items()) {
            if (mv.value()["Type"].template get<std::string>() == "Windkessel_vertex") {
                std::string id = mv.key();
                double diameter = mv.value()["diameter"].template get<double>();
                diameter_info[id] = diameter;
            }
        }

        id_1 = diameter_info.begin() -> first;
        d_1 =  diameter_info.begin() -> second;

        for (const auto &item : diameter_info)
            diameter_sum += pow(item.second / d_1, p);
    }

    void create_wk_resistance_distribution(const double &total_res) {
        double R_1 = total_res * diameter_sum;
        std::pair<double, double> R_1_pair = split_resistance(R_1);
        config["MetaVertices"][id_1]["R1"] = R_1_pair.first;
        config["MetaVertices"][id_1]["R2"] = R_1_pair.second;
        for (const auto &item : diameter_info) {
            double R = R_1 * pow(d_1 / item.second, p);
            std::pair<double, double> R_pair = split_resistance(R);
            config["MetaVertices"][item.first]["R1"] = R_pair.first;
            config["MetaVertices"][item.first]["R2"] = R_pair.second;
        }
    }

    void create_wk_compliance_distribution(const double &total_comp) {
        double C_1 = total_comp / diameter_sum;
        config["MetaVertices"][id_1]["C"] = C_1;
        for (const auto &item : diameter_info) {
            double C = C_1 * pow(item.second / d_1, p);
            config["MetaVertices"][item.first]["C"] = C;
        }
    }

    void create_wk_distribution(const double &total_res, const double &total_comp, const double &p_out, const double &coeff_p0, const double &coeff_pwv) {
        create_wk_resistance_distribution(total_res);
        create_wk_compliance_distribution(total_comp);
        set_p_out(p_out);
        set_p0(coeff_p0);
        set_pwv(coeff_pwv);
        write_json(config_path, config);
    }

    double find_total_compliance() {
        double total_comp = 0.0;
        for (const auto &item : diameter_info) {
            total_comp += (double)(config["MetaVertices"][item.first]["C"]);
        }
        return total_comp;
    }

    double find_total_resistance() {
        double reverse_total_res = 0.0;
        for (const auto &item : diameter_info) {
            double R = (double)(config["MetaVertices"][item.first]["R1"]) + (double)(config["MetaVertices"][item.first]["R2"]);
            reverse_total_res += 1.0 / R;
        }
        return 1.0 / reverse_total_res;
    }

};

#endif
