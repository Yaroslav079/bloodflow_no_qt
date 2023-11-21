#include "task.h"
#include <iostream>
#include <fstream>
#include <omp.h>
#include "csv_reader.h"
#include <filesystem>

#include <json.hpp>
using json = nlohmann::json;

Task::~Task()
{
    for (auto &e: edges)
        delete e;
    for (auto &v: vertices)
        delete v;
    for (auto &v: s_vertices)
        delete v;
}

bool Task::is_valid()
{//TODO
    return virtual_time < time_max ;//&& dt > 0;
}

void Task::run_full(bool detailed_dump)
{//TODO
    while (is_valid()) {
        run_step(0, detailed_dump);
    }

}

json config;
json heart_config;
void Task::set_up(const std::string & main_config_path, const std::string & heart_additional_config_path)
{
    Eigen::setNbThreads(1);
    true_0d_heart = nullptr;
    heart = nullptr;

    std::ifstream configFile;
    configFile.open(main_config_path);
    if (!configFile.is_open())
        throw(GraphConfigError("Cannot open main config file"));

    configFile >> config;
    configFile.close();

    std::ifstream configFile_heart;
    configFile_heart.open(heart_additional_config_path);
    if (!configFile_heart.is_open())
        throw(GraphConfigError("Cannot open heart config file"));
    configFile_heart >> heart_config;
    configFile_heart.close();

    virtual_time = config["Time0"].get<double>();
    computation_time = config["TimeComp0"].get<double>();
    Courant_number = config["CourantNumber"].get<double>();
    double dt = config["InitialTimeStep"].get<double>();
    max_abs_eigenvalue_div_dx = Courant_number/dt;


    density = config["BloodDensity"].get<double>();
    viscosity = config["BloodViscosity"].get<double>();
    zeta = config["Zeta"].get<double>();
    dx_step = config["dx"].get<double>();
    gamma = config["NumSchemeGamma"].get<double>();
    terminalResistanceScaler = config["terminalResistanceScaler"].get<double>();


    auto adan56 = config.find("ADAN56_Global");
    double E_Young = 0, P_exterior = 0, P_diastolic = 0, ca = 0, cb = 0, cc = 0, cd = 0;
    if (adan56 != config.end()) {
        E_Young = adan56.value()["E_Young"].get<double>();
        P_exterior = adan56.value()["P_exterior"].get<double>();
        P_diastolic = adan56.value()["P_diastolic"].get<double>();
        ca = adan56.value()["ca"].get<double>();
        cb = adan56.value()["cb"].get<double>();
        cc = adan56.value()["cc"].get<double>();
        cd = adan56.value()["cd"].get<double>();
    }

    for (auto simpleVertex : config["SimpleVertices"].items()) {
        enum {x_coord = 0, y_coord = 1};
        Simple_vertex * sv = new Simple_vertex(simpleVertex.key(), simpleVertex.value()[x_coord].get<double>(), simpleVertex.value()[y_coord].get<double>());
        s_vert_map[sv->get_id()] = sv;
        s_vertices.push_back(sv);
    }

    for (auto edge: config["Edges"].items()) {
        /* for now name of edge is a key and duplicated keys simply ignored by json parser
        if (edge_map.find(edge.key()) != edge_map.end()) {
            throw(GraphConfigError(edge.key() + ": name of this edge is not unique"));
        }
        */

        Simple_vertex * a, * b;
        try {
            a = s_vert_map.at(edge.value()["v1"].get<std::string>());
            b = s_vert_map.at(edge.value()["v2"].get<std::string>());
        } catch (const std::out_of_range &) {
            throw(GraphConfigError(edge.key() + ": cannot find its vertices among SimpleVertices"));
        }

        const std::string edge_type = edge.value()["Type"].get<std::string>();

        double p0;
        if (edge.value().contains("p0")) {
            p0 = edge.value()["p0"].get<double>();
        } else {
            p0 = 0;
            // std::cout << edge.key() << ": no p0" << std::endl;
        }



        Edge * e;
        if (edge_type == "Classic") {
            e = new EdgeClassicTubeLaw(/*name*/ edge.key(), /*first*/ a, /*last*/ b,
                                       /*dx_step*/ dx_step, /*len*/ edge.value()["length"].get<double>(),
                    /*density*/ density, /*viscosity*/ viscosity, /*zeta*/ zeta, /*gamma*/ gamma,
                    /*width*/ edge.value()["width"].get<double>(), /*c*/ edge.value()["C"].get<double>(), p0
                    );
        }
        /*
        else if (edge_type == "ADAN56") {
            if (adan56 == config.end())
                throw(GraphConfigError("No global ADAN56 data specified but there is an edge of type ADAN56"));

            auto adan_edge = new EdgeAdan56TubeLaw(edge.key(), a, b,
                        dx_step, edge.value()["length"].get<double>(),
                    density,  viscosity,  zeta, gamma,
                    P_exterior, P_diastolic, ca, cb, cc, cd, E_Young,
                    edge.value()["RadProx"].get<double>(),
                    edge.value()["RadDistal"].get<double>()
                    );
            //std::cout << adan_edge->get_id() << " Diam: " << adan_edge->mean_diameter() << std::endl;
            e = adan_edge;
        } */
        else {
            throw(GraphConfigError(edge.key() + ": unsupported type of edge"));
        }

        if (edge.value().contains("Study") && edge.value()["Study"].get<int>() != 0) {
            edges_to_study.push_back(e);
        }

        edge_map[e->get_id()] = e;
        edges.push_back(e);
        s_edges.push_back(e);
    }

    for (auto mv: config["MetaVertices"].items()) {
        get_metavertex_from_config(mv);
    }

    for (auto mv: heart_config["MetaVertices"].items()) {
        get_metavertex_from_config(mv);
    }

    double period = true_0d_heart -> get_period();
    for (auto &e: edges)
        e -> set_period(period);
    for (auto &wk: wk_vertices)
        wk -> set_period(period);
/*
    //now check whether we need to use additional heart config file
    if (heart) {
        //already have a heart from main config so we dont need additional config and
        //no SimpleVertices should be left
        if (!s_vert_map.empty())
            throw(GraphConfigError("Some SimpleVertices do not have their personal MetaVertices"));
        if (!heart_additional_config_path.empty())
            throw(GraphConfigError("Main config file already has heart model"));
    } else {
        //no heart yet so try to use additional config file
        //we also have to keep some SimpleVertices for our heart metavertex
        if (!heart_additional_config_path.empty()) {
            std::ifstream configFile;
            configFile.open(heart_additional_config_path);
            if (!configFile.is_open())
                throw(GraphConfigError("Cannot open heart config file"));
            json heart_config;
            configFile >> heart_config;
            configFile.close();
            //here you should get heart and additional metavertices like internal in case of heart without pump
            for (auto mv: heart_config["MetaVertices"].items()) {
                get_metavertex_from_config(mv);
            }
        } else {
            throw(GraphConfigError("Main config file does not have heart model, provide additional config file with heart"));
        }
        //now check did we finally get heart or not
        if (!heart)
            throw(GraphConfigError("No heart model provided in any of config files"));
        //no SimpleVertices should be left at this point
        if (!s_vert_map.empty())
            throw(GraphConfigError("Some SimpleVertices do not have their personal MetaVertices"));
    }
*/

/*
    for (auto &v: s_vertices)
        std::cout << *v << std::endl;
    for (auto &e: edges) {
        std::cout << *e << std::endl;
        e->print_info(std::cout);
    }

    if (heart) {
        for (ParameterBlock & block: heart->parameterBlocks)
            for (Parameter & param: block)
                std::cout << param.first << ": " << *param.second << std::endl;
    }




    std::cout << std::endl << std::endl;
    std::cout << "*************PRETTY PRINT*************" << std::endl;
    std::cout << "1. EDGES" << std::endl;
    for (auto &e: edges)
        e->print_info();

    std::cout << "2. WK ENDS" << std::endl;
    for (auto &wk: wk_vertices)
        wk->print_info();

    std::cout << "3. CORONARY ENDS" << std::endl;
    for (auto &c: terminal_vertices)
        c->print_info();
    std::cout << "*************PRETTY PRINT*************" << std::endl;
*/
}

void Task::run_step(double sim_timer, bool detailed_dump = false)
{
    double t1 = omp_get_wtime();
    double local_max_abs_eigenvalue_div_dx = max_abs_eigenvalue_div_dx;
#pragma omp parallel firstprivate(sim_timer, local_max_abs_eigenvalue_div_dx)
{
    double time_cur = virtual_time;
    double dt;
    while (!(sim_timer < 0)) {

        dt = Courant_number / local_max_abs_eigenvalue_div_dx;
        time_cur += dt;
        sim_timer -= dt;

#pragma omp single
{
        max_abs_eigenvalue_div_dx = 0;
}

#pragma omp for reduction(max: max_abs_eigenvalue_div_dx) schedule(dynamic, 1)
        for (unsigned long i = 0; i < edges.size(); i++) {
            edges[i]->set_dt(dt);
            edges[i]->set_T(time_cur);
            max_abs_eigenvalue_div_dx = std::max(max_abs_eigenvalue_div_dx, edges[i]->solve());
        }
local_max_abs_eigenvalue_div_dx = max_abs_eigenvalue_div_dx;

#pragma omp for schedule(dynamic, 2)
        for (unsigned long i = 0; i < vertices.size(); i++) {
            vertices[i]->set_dt(dt);
            vertices[i]->set_T(time_cur);
            vertices[i]->update_boundary();
        }
    }

#pragma omp single
{
    max_abs_eigenvalue_div_dx = local_max_abs_eigenvalue_div_dx;
    virtual_time = time_cur;
}

}
    double t3 = omp_get_wtime();

    computation_time += t3 - t1;
    for (auto &e: edges) {
        if (e->get_id() == "brachial-L") {
            e->print_info(fout_brachial);
        }
    }
    
    if (detailed_dump) {
        for (auto &e: edges) {
            std::ofstream fout = set_path_to_dump(e -> get_id());
            e->print_info(fout);
            fout.close();
        }
        for (auto &wk: wk_vertices) {
            std::ofstream fout = set_path_to_dump(wk -> get_id());
            wk->print_info(virtual_time, fout);
            fout.close();
        }
        std::ofstream fout = set_path_to_dump("heart");
        true_0d_heart -> print_info(fout);
        fout.close();
    }
}


void Task::find_max_abs_eigenvalue_div_dx()
{
    max_abs_eigenvalue_div_dx = 0;
    for (auto &edge: edges) {
        const double tmp_e = edge->get_max_abs_eigenvalue_div_dx();
        if (max_abs_eigenvalue_div_dx < tmp_e)
            max_abs_eigenvalue_div_dx = tmp_e;
    }
}
/*
void Task::find_dt()
{
    find_max_abs_eigenvalue_div_dx();
    if (max_abs_eigenvalue_div_dx > 1e-10) {
        double dt = Courant_number / max_abs_eigenvalue_div_dx;
    } else {
        double dt = -1;
        std::cout << "scheme stopped" << std::endl;
        throw("Scheme stopped");
    }
}
*/


int Task::get_edge_data_size(const std::string & ename, double & len)
{
    for (unsigned i = 0; i < edges.size(); i++) {
        if (edges[i]->get_id() == ename) {
            len = edges[i]->get_len();
            return edges[i]->get_points_num();
        }
    }
    return 0;
}

double Task::get_simTime() const
{
    return virtual_time;
}
double Task::get_computationTime() const
{
    return computation_time;
}

void Task::get_LV_pv(double &p, double &v) const
{
    if (true_0d_heart) {
        p = true_0d_heart->get_LV_P();
        v = true_0d_heart->get_LV_V();
    }
}

void Task::get_LA_pv(double &p, double &v) const
{
    if (true_0d_heart) {
        p = true_0d_heart->get_LA_P();
        v = true_0d_heart->get_LA_V();
    }
}

/*
void Task::set_pump_LV_pressure_mmHg(const double &p)
{
    //Obsolete
    //if (vertexPump)
     //  ;// vertexPump->set_LV_pressure_mmHg(p);

}
*/
double Task::get_strokeVolume()
{
    if (heart)
        return heart->get_stroke_volume_av();

    else
        return -999;
}

double Task::get_aortic_reg() {
    return heart -> get_aortic_regurgitation_fraction();
}

double Task::get_mitral_reg() {
    return heart -> get_mitral_regurgitation_fraction();
}

double Task::get_aortic_valve()
{
    if (true_0d_heart)
        return true_0d_heart->get_aortic_valve();
    else
        return 0;
}

double Task::get_mitral_valve()
{
    if (true_0d_heart)
        return true_0d_heart->get_mitral_valve();
    else
        return 0;
}


void Task::generate_report_speed()
{
    if (false) {
        std::ofstream speedStatsStatus("speedStatsStatus.dat"), speedStats("speedStats.dat");
        speedStatsStatus << "% RPM  | is_pump_backflow | is_av_always_closed | is_mv_always_opened" << std::endl;
        speedStats << "RPM  sv-aortic-valve  sv-pump  sv-total  av-time-opened  mv-time-opened  pulse-pressure-root \
 max-pressure-root  min-pressure-root  min-pump-flow  lv-work";
        for (auto &e: edges_to_study) {
            speedStats << "  " << e->get_id() + "-mean-flow" << "  " << e->get_id() + "-mean-pressure" << "  " << e->get_id() + "-mean-speed";
        }
        speedStats << std::endl;

        double first_speed = 0, last_speed = 20000, step_speed = 500;
        for (double speed = first_speed; speed < last_speed + 1; speed += step_speed) {
            // pump_0d_heart->set_w(speed);
            bool is_pump_backflow, is_av_always_closed,is_mv_always_opened;
            double sv_av, sv_pump;
            double av_time_opened, mv_time_opened;
            double pulse_pressure;
            double max_pressure, min_pressure, min_pump_flow, lv_work;
            generate_report(std::to_string(int(speed)), is_pump_backflow, is_av_always_closed,
                            is_mv_always_opened, sv_av, sv_pump, av_time_opened, mv_time_opened, pulse_pressure,
                            max_pressure, min_pressure, min_pump_flow, lv_work
                            );

            speedStatsStatus << speed << " | " << is_pump_backflow << " | " << is_av_always_closed << " | "
                       << is_mv_always_opened << std::endl;

            speedStats << speed << " " << sv_av << " " << sv_pump << " " << sv_av + sv_pump
                       << " " << av_time_opened << " " << mv_time_opened << " " << pulse_pressure << " " <<
                          max_pressure << " " << min_pressure << " " <<  min_pump_flow << " " << lv_work;

            for (auto &e: edges_to_study) {
                speedStats << " " << e->get_mean_center_flow() << " " << e->get_mean_pressure()/1333.22 << " " << e->get_mean_center_speed();
            }
            speedStats << std::endl;
        }
        speedStats.close();
        speedStatsStatus.close();
    }
}
void Task::generate_report_test_LV_ESPVR()
{
    if (true_0d_heart) {
        std::ofstream espvrStats("espvrStats.dat");
        espvrStats << "%ESPVR | total_sv | av_time_opened | mv_time_opened | pulse_pressure | max_pressure | min_pressure | lv_work" << std::endl;
        const double espvr = true_0d_heart->get_ESPVR();
        for (int denom = 1; espvr/denom > 0.1; denom *= 2) {
            double new_espvr = espvr/denom;
            true_0d_heart->set_ESPVR(new_espvr);
            bool is_pump_backflow, is_av_always_closed, is_mv_always_opened;
            double sv_av, sv_pump;
            double av_time_opened, mv_time_opened;
            double pulse_pressure;
            double max_pressure, min_pressure, min_pump_flow, lv_work;
            generate_report(std::to_string(denom) + "ESVPR_part", is_pump_backflow, is_av_always_closed,
                            is_mv_always_opened, sv_av, sv_pump, av_time_opened, mv_time_opened, pulse_pressure, max_pressure, min_pressure, min_pump_flow, lv_work);

            espvrStats << new_espvr << " " << sv_av + sv_pump << " " << av_time_opened << " " << mv_time_opened << " " << pulse_pressure
                       << " " << max_pressure << " " << min_pressure << " " << " " << lv_work
                       << std::endl;
        }
        //restore
        true_0d_heart->set_ESPVR(espvr);
        espvrStats.close();
    }
}
void Task::generate_report_test_pulm_vein_pressure()
{
    if (false) {
        double first_speed = 0, last_speed = 20000, step_speed = 1000;
        for (double speed = first_speed; speed < last_speed + 1; speed += step_speed) {
            std::ofstream pveinStats(std::string("pvein") + std::to_string(int(speed)) + ".dat");
            pveinStats << "%Pvein | total_sv | av_time_opened | mv_time_opened | pulse_pressure | max_pressure | min_pressure | lv_work" << std::endl;
            const double pvein_orig = true_0d_heart->get_PveinPressure();
            for (int pvein = 5; pvein < 18; pvein++) {
                // pump_0d_heart->set_w(speed);
                true_0d_heart->set_PveinPressure(pvein);
                bool is_pump_backflow, is_av_always_closed, is_mv_always_opened;
                double sv_av, sv_pump;
                double av_time_opened, mv_time_opened;
                double pulse_pressure;
                double max_pressure, min_pressure, min_pump_flow, lv_work;
                generate_report(std::to_string(pvein) + "Pvein", is_pump_backflow, is_av_always_closed,
                                is_mv_always_opened, sv_av, sv_pump, av_time_opened, mv_time_opened, pulse_pressure, max_pressure, min_pressure, min_pump_flow, lv_work);

                pveinStats << pvein << " " << sv_av + sv_pump << " " << av_time_opened << " " << mv_time_opened << " " << pulse_pressure
                           << " " << max_pressure << " " << min_pressure << " " << lv_work
                           << std::endl;
            }
            true_0d_heart->set_PveinPressure(pvein_orig);
        }
        return;
    }



    if (true_0d_heart) {
        std::ofstream pveinStats("pveinPStats.dat");
        pveinStats << "%Pvein | total_sv | av_time_opened | mv_time_opened | pulse_pressure | max_pressure | min_pressure | lv_work" << std::endl;
        const double pvein_orig = true_0d_heart->get_PveinPressure();
        for (int pvein = 5; pvein < 18; pvein++) {
            true_0d_heart->set_PveinPressure(pvein);
            bool is_pump_backflow, is_av_always_closed, is_mv_always_opened;
            double sv_av, sv_pump;
            double av_time_opened, mv_time_opened;
            double pulse_pressure;
            double max_pressure, min_pressure, min_pump_flow, lv_work;
            generate_report(std::to_string(pvein) + "Pvein", is_pump_backflow, is_av_always_closed,
                            is_mv_always_opened, sv_av, sv_pump, av_time_opened, mv_time_opened, pulse_pressure, max_pressure, min_pressure, min_pump_flow, lv_work);

            pveinStats << pvein << " " << sv_av + sv_pump << " " << av_time_opened << " " << mv_time_opened << " " << pulse_pressure
                       << " " << max_pressure << " " << min_pressure << " " << lv_work
                       << std::endl;
        }
        true_0d_heart->set_PveinPressure(pvein_orig);
        pveinStats.close();
    }
}

void Task::generate_report_test_afterload()
{
    if (false) {
        double first_speed = 0, last_speed = 20000, step_speed = 1000;
        for (double speed = first_speed; speed < last_speed + 1; speed += step_speed) {
            // pump_0d_heart->set_w(speed);
            test_afterload(std::string("afterload") + std::to_string(int(speed)) + ".dat");
        }
    } else if (true_0d_heart) {
        test_afterload("afterload.dat");
    }
}

void Task::test_afterload(const std::string &output)
{
    //assume only one end wk vertex
    Windkessel_vertex * wk = wk_vertices[0];
    std::ofstream afterloadStats(output);
    //afterloadStats << "%R1 | total_sv" << std::endl;
    afterloadStats << "%R2 | total_sv" << std::endl;

    //for (int R1 = 30; R1 < 100; R1 += 10) {
    for (int R2 = 1000; R2 < 2000; R2 += 100) {
        //wk->set_R1(R1);
        wk->set_R2(R2);
        run_step(10);
        afterloadStats << R2 << " " << heart->get_stroke_volume_total() << std::endl;
        //afterloadStats << R1 << " " << pump_0d_heart->get_stroke_volume_total() << std::endl;
    }
    wk->set_R2(1500);
    //wk->set_R1(60);
}

void Task::generate_report_basic()
{
    if (true_0d_heart) {
        std::ofstream heartSV(true_0d_heart->get_id() + "_heart.dat");
        heartSV << "sv-total  av-time-opened  mv-time-opened  pulse-pressure-root  max-pressure-root  min-pressure-root  lv-work";

        for (auto &e: edges_to_study) {
            heartSV << "  " << e->get_id() + "-mean-flow" << "  " << e->get_id() + "-mean-pressure" << "  " << e->get_id() + "-mean-speed";
        }
        heartSV << std::endl;

        bool is_pump_backflow, is_av_always_closed, is_mv_always_opened;
        double sv_av, sv_pump;
        double av_time_opened, mv_time_opened;
        double pulse_pressure;
        double max_pressure, min_pressure, min_pump_flow, lv_work;
        generate_report(true_0d_heart->get_id(), is_pump_backflow, is_av_always_closed,
                        is_mv_always_opened, sv_av, sv_pump, av_time_opened, mv_time_opened, pulse_pressure, max_pressure, min_pressure, min_pump_flow, lv_work);
        heartSV << sv_av << " " << av_time_opened << " " << mv_time_opened << " " << pulse_pressure
                << " " << max_pressure << " " << min_pressure << " " << lv_work;
        for (auto &e: edges_to_study) {
            heartSV << " " << e->get_mean_center_flow() << " " << e->get_mean_pressure()/1333.22 << " " << e->get_mean_center_speed();
        }
        heartSV << std::endl;
        heartSV.close();
    }
}

void Task::generate_report(const std::string &outputPrefix, bool &is_pump_backflow, bool &is_av_always_closed,
                           bool &is_mv_always_opened, double & sv_av, double & sv_pump, double & aortic_valve_opened_time,
                           double & mitral_valve_opened_time, double & pulse_pressure,
                           double & max_pressure, double & min_pressure,
                           double & min_pump_flow, double & lv_work)
{
    assert(true_0d_heart);
    std::ofstream report(outputPrefix + std::string("_report.dat"));
    if (!report.is_open()) {
        throw("Report file is not opened");
    }
    report << "time(s) LV-V LV-P LA-V LA-P flow-av av-angle mi-angle P-aortic-root H-pump Q-pump Q-total";

    for (auto &e: edges_to_study)
        report << " " << e->get_id() + "-center-pressure"; ///@todo should specify what value to dump in config file, but now if you want speed instead of pressure, you should change the string to "-center-speed"
    report << std::endl;

    double skipTime = 6 - fmod(virtual_time, 1);
    double recordTime = 1; //one cycle

 //   double startTime = time_cur + skipTime;
  //  double endTime = startTime + recordTime;
   // while (time_cur < startTime) {
    //    run_step(1);
    //}
    run_step(skipTime);
    double startTime = virtual_time - fmod(virtual_time, 1);
    double endTime = startTime + recordTime;
    double timeStep = 0.01;

    is_pump_backflow = 0;
    is_av_always_closed = 1;
    is_mv_always_opened = 1;
    sv_av = 0;
    sv_pump = 0;
   // int i = 0;
    while (virtual_time < endTime + timeStep) {
      //  if (time_cur - startTime + 1e-15 < i*timeStep) {
       //     run_step(1);
       //     continue;
       // }
       // i++;
        report <<  virtual_time - startTime
                << " " << true_0d_heart->get_LV_V()
                << " " << true_0d_heart->get_LV_P()
                << " " << true_0d_heart->get_LA_V()
                << " " << true_0d_heart->get_LA_P()
                << " " << true_0d_heart->get_flow_av()
                << " " << true_0d_heart->get_aortic_valve()
                << " " << true_0d_heart->get_mitral_valve()
                << " " << true_0d_heart->get_aortic_root_pressure();
                //<< " " << edges[1]->Equation_of_state(edges[1]->get_V_s(s_vertices.back(), 0), s_vertices.back(), 0)/1333.22;
            //smelly
        report << " " << 0
               << " " << 0
               << " " << 0;

        for (auto &e: edges_to_study) {
            report << " " <<  e->get_center_pressure(); ///@todo see another todo earlier in this function, if you want to dump speed instead of pressure, call get_center_speed() method
        }

        report << std::endl;

        if (true_0d_heart->get_aortic_valve() > 1)
            is_av_always_closed = 0;
        if (true_0d_heart->get_mitral_valve() < 1)
            is_mv_always_opened = 0;
        run_step(timeStep);
    }
    sv_av = true_0d_heart->get_stroke_volume_av();
    aortic_valve_opened_time = true_0d_heart->get_time_aortic_valve_opened();
    mitral_valve_opened_time = true_0d_heart->get_time_mitral_valve_opened();
    pulse_pressure = true_0d_heart->get_pulsePressure();

    max_pressure = true_0d_heart->get_max_aortic_pressure();
    min_pressure = true_0d_heart->get_min_aortic_pressure();
    lv_work = true_0d_heart->get_LV_work();

    report.close();
}

void Task::generate_report_edge_stats()
{
    //start testing
    run_step(10);

    //write to the file
    std::ofstream report(std::string("edges_stats_report.dat"));
    report << "edge | max_center_speed | min_center_speed | mean_pressure | mean_center_flow" << std::endl;
    for (auto &e: edges) {
        report << e->get_id() << " " << e->get_max_center_speed()
               << " " << e->get_min_center_speed() << " " << e->get_mean_pressure()/1333.22 << " " << e->get_mean_center_flow() << std::endl;
    }

}

void variationWithPerm(std::vector<int> & v)
{//-1, 0 or 1 only
    assert(!v.empty());
    v[0] += 1;
    for (unsigned i = 0; i < v.size(); i++) {
        if (v[i] == 2) {
            v[i] = -1;
            if (i != v.size() - 1)
                v[i+1] += 1;
        } else {
            break;
        }
    }
}

void Task::sensitivityAnalysis()
{
    if (heart) {
        std::cout << "########### SENSITIVITY ANALYSIS ###########" << std::endl;
        const double factor = 0.5;
        std::cout << "Factor: " << factor << std::endl;
        //find reference output for default parameters
        run_step(5);
        double ref_sv = heart->get_stroke_volume_av();
        std::cout << "Reference: " << "sv: " << ref_sv << std::endl;
        for (ParameterBlock & block: heart->parameterBlocks) {
            std::cout << "ParameterBlock" << std::endl;

            std::vector<int> var(block.size(), -1);

            for (Parameter & param: block)
                std::cout << param.first << " ";
            std::cout << std::endl;

            //save default parameters
            std::vector<double> default_block(block.size());
            for (unsigned i = 0; i < block.size(); i++)
                default_block[i] = *block[i].second;

            //iterate over all possible variations
            for (int i = 0; i < pow(3, var.size()); i++) {
                //change parameters in block
                for (unsigned j = 0; j < block.size(); j++)
                    *block[j].second = default_block[j] * (1 + factor * var[j]);

                run_step(5);

                for (unsigned j = 0; j < block.size(); j++)
                    std::cout << var[j] << " " << *block[j].second << " ";
                std::cout << "sv: " << heart->get_stroke_volume_av() << std::endl;
                variationWithPerm(var);
            }
            //restore default parameters
            for (unsigned i = 0; i < block.size(); i++)
                *block[i].second =  default_block[i];

        }
        std::cout << "########### COMPLETE ###########" << std::endl;
    }
}

void Task::generate_report_wk_flows()
{

/*
    if (true_0d_heart) {
        const double espvr = true_0d_heart->get_ESPVR();
        for (int denom = 1; espvr/denom > 0.1; denom *= 2)


    {

                double new_espvr = espvr/denom;
                true_0d_heart->set_ESPVR(new_espvr);

*/


    //start testing
    //run_step(10);

    //write to the file
    //std::ofstream report(std::to_string(denom) + "ESVPR_WK.dat");
    std::ofstream report(std::string("wk_flows_report.dat"));
    double input_sv;
    //double total_out = 0;
    input_sv = heart->get_stroke_volume_av();
    report << "%Heart SV: " << input_sv << std::endl;

  /*  report << "%wk & SV, ml & percentage" << std::endl;
    for (auto & term_v: terminal_vertices) {
        const double term_sv = - term_v->get_stroke_volume();
        total_out += term_sv;
        report << term_v->get_id() << " & " << term_sv << " & "
               << term_sv / input_sv * 100 <<
               "\\\\" << std::endl;
    }
    report << "%Total outflow: " << total_out << std::endl;
    */


    report << "%vessel & SV, ml & percentage" << std::endl;

    for (auto & e: edges) {
        report << e->get_id() << " & " << e->get_mean_center_flow() << " & "
               << fabs(e->get_mean_center_flow())*100/input_sv << " \\\\" << std::endl;
    }

/*

    }



        //restore
        true_0d_heart->set_ESPVR(espvr);
    }
*/
}


GraphConfigError::GraphConfigError(const std::string &what_arg)
    :std::logic_error(what_arg)
{}

GraphConfigError::GraphConfigError(const char *what_arg)
    :std::logic_error(what_arg)
{}

void Task::dump_data() {
    Csv_Reader reader(path_to_brachial_data.c_str());
    std::map <double, double> p_t = reader.getData();

    double P_d = 1e4;
    double P_s = 0.0;
    for (auto itr = p_t.begin(); itr != p_t.end(); ++itr) {
        if (itr -> second > P_s)
            P_s = itr -> second;
        if (itr -> second < P_d)
            P_d = itr -> second;
    }

    P_sys = P_s;
    P_dis = P_d;

    std::cout << P_sys << "," << P_dis << "," << get_strokeVolume() << "," << get_aortic_reg() * 100.0 << "," << get_mitral_reg() * 100.0 << std::endl;
}
