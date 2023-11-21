#ifndef TASK_H
#define TASK_H

#include <Eigen/Dense>
#include <unordered_map>
#include <exception>
#include "edge.h"
#include "vertex.h"
#include "rcrwindkessel.h"
#include "heart_advanced_valves.h"
#include "heart_aortic_reg.h"
#include <fstream>

/**
 * @brief Error class to throw while reading config file with vascular graph
 */
class GraphConfigError
: public std::logic_error
{
public:
  explicit GraphConfigError(const std::string& what_arg);
  explicit GraphConfigError(const char* what_arg);
};


/**
 * @brief Pop key-value pair from map and return value
 */
template<typename Map, typename Key>
auto map_pop(Map & map, const Key & key)
{
    auto iter = map.find(key);
    if (iter == map.end()) {
        throw std::out_of_range("map_pop: no key was found");
    } else {
        const auto res = iter->second;
        map.erase(iter);
        return res;
    }
}


/**
 * @brief Class for problem initialization and run
 */
struct Task
{
    typedef ::Edge<Eigen::Matrix, double, Eigen::Dynamic> Edge;
    typedef ::EdgeClassicTubeLaw<Eigen::Matrix, double, Eigen::Dynamic>  EdgeClassicTubeLaw;
    // typedef ::EdgeAdan56TubeLaw<Eigen::Matrix, double, Eigen::Dynamic>  EdgeAdan56TubeLaw;

    typedef ::Vertex<Edge, Eigen::Matrix, double, Eigen::Dynamic> Vertex;
    typedef ::Internal_vertex<Edge, Eigen::Matrix, double, Eigen::Dynamic> Internal_vertex;
    //typedef ::Free_flow_vertex<Edge, Eigen::Matrix, double, Eigen::Dynamic> Free_flow_vertex;
    typedef ::TerminalVertex<Edge, Eigen::Matrix, double, Eigen::Dynamic> TerminalVertex;
    typedef ::RCR_Windkessel_vertex<Edge, Eigen::Matrix, double, Eigen::Dynamic> Windkessel_vertex;
    typedef ::TerminalResistanceCoronaryArteries<Edge, Eigen::Matrix, double, Eigen::Dynamic> TerminalResistanceCoronaryArteries;

    typedef ::Heart<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart;
    typedef ::Constant_flow_heart<Edge, Eigen::Matrix, double, Eigen::Dynamic> Constant_flow_heart;
    typedef ::Simple_heart<Edge, Eigen::Matrix, double, Eigen::Dynamic> Simple_heart;

    typedef ::True_0d_heart<Edge, Eigen::Matrix, double, Eigen::Dynamic> True_0d_heart;
    typedef ::Heart_AdValves<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_AdValves;
    typedef ::Heart_Aortic_Reg<Edge, Eigen::Matrix, double, Eigen::Dynamic> Heart_Aortic_Reg;

    double virtual_time; // time inside a model
    double time_max = 8.942; // max time inside a model

    double computation_time; // time spent on computations
    double Courant_number;
    double max_abs_eigenvalue_div_dx;

    std::vector<Edge *> edges;
    std::vector<Simple_edge *> s_edges; // pointers to Edge objects from edges vector but casted to Simple_edge

    std::vector<Simple_vertex *> s_vertices;
    std::vector<Vertex *> vertices;

    std::vector<Windkessel_vertex *> wk_vertices;
    std::vector<TerminalVertex *> terminal_vertices;

    std::unordered_map<std::string, Simple_vertex *> s_vert_map;
    std::unordered_map<std::string, Edge *> edge_map;

    std::vector<Edge *> edges_to_study;

    Heart * heart; ///@todo pointer is stored in vertices vector which is cleaned in the Task destructor but better to use smart pointers
    True_0d_heart * true_0d_heart;
    // std::ofstream fout_dump;

    double P_sys;
    double P_dis;

private:
    double density, viscosity, zeta, dx_step, gamma, terminalResistanceScaler;
    std::string path_to_brachial_data;
    std::ofstream fout_brachial;
public:

    template<typename Mjson>
    void get_metavertex_from_config(const Mjson &mv)
    {
        double w = 0; // default pump speed
        const std::string type = mv.value()["Type"].template get<std::string>();
        if (type == "Heart_AdValves") {
            Edge * e;
            try {
                e = edge_map.at(mv.value()["edge"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its edge among Edges"));
            }
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            if (!e->test_connection(sv))
                throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));

            double heart_period = 60.0 / mv.value()["Heart_rate"].template get<double>();
            true_0d_heart =  new Heart_AdValves(mv.key(), e, sv, density, viscosity,
                    mv.value()["L_av"].template get<double>(),
                    mv.value()["L_pu"].template get<double>(),
                    mv.value()["L_mi"].template get<double>(),
                    mv.value()["B_pu_const"].template get<double>(),
                    mv.value()["B_av_denomin"].template get<double>(),
                    mv.value()["B_mi_denomin"].template get<double>(),
                    mv.value()["pulmVeinsPressure"].template get<double>(),
                    mv.value()["LV_initialVolume"].template get<double>(),
                    mv.value()["LV_V0"].template get<double>(),
                    mv.value()["LA_V0"].template get<double>(),
                    mv.value()["LV_ESPVR"].template get<double>(),
                    mv.value()["LV_EDPVR"].template get<double>(),
                    mv.value()["LA_ESPVR"].template get<double>(),
                    mv.value()["LA_EDPVR"].template get<double>(),
                    mv.value()["valvePressureForceCoeff"].template get<double>(),
                    mv.value()["valveFrictionalForceCoeff"].template get<double>(),
                    mv.value()["LV_inertiaCoeff"].template get<double>(),
                    mv.value()["LV_dynamicResistanceCoeff"].template get<double>(),
                    mv.value()["LA_inertiaCoeff"].template get<double>(),
                    mv.value()["LA_dynamicResistanceCoeff"].template get<double>(),
                    heart_period,
                    heart_period * mv.value()["T_sys"].template get<double>()
                    );
            vertices.push_back(true_0d_heart);
            heart = true_0d_heart;
            time_max = 10.0 * heart_period;
        }
        else if (type == "Heart_Aortic_Reg") {
            Edge * e;
            try {
                e = edge_map.at(mv.value()["edge"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its edge among Edges"));
            }
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            if (!e->test_connection(sv))
                throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));

            double heart_period = 60.0 / mv.value()["Heart_rate"].template get<double>();
            true_0d_heart =  new Heart_Aortic_Reg(mv.key(), e, sv, density, viscosity,
                    mv.value()["L_av"].template get<double>(),
                    mv.value()["L_pu"].template get<double>(),
                    mv.value()["L_mi"].template get<double>(),
                    mv.value()["B_pu_const"].template get<double>(),
                    mv.value()["B_av_denomin"].template get<double>(),
                    mv.value()["B_mi_denomin"].template get<double>(),
                    mv.value()["pulmVeinsPressure"].template get<double>(),
                    mv.value()["LV_initialVolume"].template get<double>(),
                    mv.value()["LV_V0"].template get<double>(),
                    mv.value()["LA_V0"].template get<double>(),
                    mv.value()["LV_ESPVR"].template get<double>(),
                    mv.value()["LV_EDPVR"].template get<double>(),
                    mv.value()["LA_ESPVR"].template get<double>(),
                    mv.value()["LA_EDPVR"].template get<double>(),
                    mv.value()["valvePressureForceCoeff"].template get<double>(),
                    mv.value()["valveFrictionalForceCoeff"].template get<double>(),
                    mv.value()["LV_inertiaCoeff"].template get<double>(),
                    mv.value()["LV_dynamicResistanceCoeff"].template get<double>(),
                    mv.value()["LA_inertiaCoeff"].template get<double>(),
                    mv.value()["LA_dynamicResistanceCoeff"].template get<double>(),
                    heart_period,
                    heart_period * mv.value()["T_sys"].template get<double>()
                    );
            vertices.push_back(true_0d_heart);
            heart = true_0d_heart;
            time_max = 10.0 * heart_period;
        }
        else if (type == "Internal") {
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            std::vector<Edge*> this_vertex_edges;
            for (auto ej: mv.value()["edges"]) {
                Edge * e;
                try {
                    e = edge_map.at(ej.template get<std::string>());
                } catch (const std::out_of_range &) {
                    throw(GraphConfigError(mv.key() + ": cannot find its edges among Edges"));
                }
                if (!e->test_connection(sv))
                    throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));
                this_vertex_edges.push_back(e);
            }
            vertices.push_back(new Internal_vertex(mv.key(), this_vertex_edges, sv));
        } else if (type == "Windkessel_vertex") {
            Edge * e;
            try {
                e = edge_map.at(mv.value()["edge"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its edge among Edges"));
            }
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            if (!e->test_connection(sv))
                throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));

            double Pout = mv.value()["P_out"].template get<double>();
            double R1 = terminalResistanceScaler * mv.value()["R1"].template get<double>();
            double R2 = terminalResistanceScaler * mv.value()["R2"].template get<double>();
            double C = mv.value()["C"].template get<double>();
            auto * wk = new Windkessel_vertex(mv.key(), e, sv, Pout, R1, R2, C);
            wk_vertices.push_back(wk);
            //terminal_vertices.push_back(wk);
            vertices.push_back(wk);
        } else if (type == "TerminalResistanceCoronary") {
            Edge * e;
            try {
                e = edge_map.at(mv.value()["edge"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its edge among Edges"));
            }
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            if (!e->test_connection(sv))
                throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));

            double Pout = mv.value()["P_out"].template get<double>();
            double R = terminalResistanceScaler * mv.value()["R"].template get<double>();
            auto * terminalCoronary = new TerminalResistanceCoronaryArteries(mv.key(), e, sv, Pout, R);
            terminal_vertices.push_back(terminalCoronary);
            vertices.push_back(terminalCoronary);
        } else if (type == "SimpleFlow") {
            Edge * e;
            try {
                e = edge_map.at(mv.value()["edge"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its edge among Edges"));
            }
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            if (!e->test_connection(sv))
                throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));

            double heart_period = 60.0 / mv.value()["Heart_rate"].template get<double>();
            heart =  new Simple_heart(mv.key(), e, sv,
                                      /*Period*/ heart_period,
                                      /*Stroke volume*/ mv.value()["Stroke_volume"].template get<double>());
            vertices.push_back(heart);
        } else if (type == "ConstantFlow") {
            Edge * e;
            try {
                e = edge_map.at(mv.value()["edge"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its edge among Edges"));
            }
            Simple_vertex * sv;
            try {
                sv = map_pop(s_vert_map, mv.value()["vertex"].template get<std::string>());
                //sv = s_vert_map.at(mv.value()["vertex"].template get<std::string>());
            } catch (const std::out_of_range &) {
                throw(GraphConfigError(mv.key() + ": cannot find its vertex among SimpleVertices"));
            }
            if (!e->test_connection(sv))
                throw(GraphConfigError(mv.key() + ": edge " + e->get_id() + " does not have a required Simple_vertex"));

            heart =  new Constant_flow_heart(mv.key(), e, sv,
                                             /*Flowrate*/ mv.value()["Flowrate"].template get<double>());
            vertices.push_back(heart);
        } else {
            throw(GraphConfigError(mv.key() + ": unsupported type of MetaVertex"));
        }
    }

    void find_max_abs_eigenvalue_div_dx();
    //void find_dt();
    /**
     * @brief run_step runs model for sim_timer seconds
     * @param sim_timer step time (in seconds)
     */
    void run_step(double sim_timer, bool detailed_dump);

    /**
     * @brief set_up
     * @param main_config_path
     * @param heart_additional_config_path
     */
    void set_up(const std::string & main_config_path, const std::string & heart_additional_config_path = std::string());
    void run_full(bool detailed_dump = false);
    bool is_valid();
    ~Task();

    int get_edge_data_size(const std::string & ename, double & len);
  //  template<typename C>
  //  void get_edge_data(const std::string & ename, value_type vt, C & array, int sz);//array of doubles

    template<typename C, typename Point>
    C get_edge_data(const std::string ename, quantity vt, C array, C end, double dx_add)
    {
        //double t1 = omp_get_wtime();
        Edge *edge = edge_map[ename];
        const double dx = edge->get_dx();
        const auto sz = edge->get_points_num();
        assert(sz > 0);
        switch (vt) {
            case area: {
                for (int j = 0; j < sz && array != end; j++, ++array) {
                    *array = Point(j*dx + dx_add, edge->get_S_unsafe(j));
                }
                break;
            }

            case speed: {
                for (int j = 0; j < sz && array != end; j++, ++array) {
                    *array = Point(j*dx + dx_add, edge->get_U_unsafe(j));
                }
                break;
            }

            case flow: {
                for (int j = 0; j < sz && array != end; j++, ++array) {
                    *array = Point(j*dx + dx_add, edge->get_Q_unsafe(j));
                }
                break;
            }

            case pressure: {// in mm Hg
                for (int j = 0; j < sz && array != end; j++, ++array) {
                    *array = Point(j*dx + dx_add, edge->get_P_unsafe(j) / 1333.22);
                }
                break;
            }

        }
        //eval_time += omp_get_wtime() - t1;
        return array;
    }
    double get_simTime() const;
    double get_computationTime() const;
    void get_LV_pv(double &p, double &v) const;
    void get_LA_pv(double &p, double &v) const;

    //void set_pump_LV_pressure_mmHg(const double &p);
    double get_strokeVolume();
    double get_aortic_reg();
    double get_mitral_reg();
    double get_aortic_valve();
    double get_mitral_valve();
    //double tprev;
    void generate_report_speed();
    void generate_report_test_LV_ESPVR();
    void generate_report_test_pulm_vein_pressure();
    void generate_report_test_afterload();
    void test_afterload(const std::string & output);
    void generate_report_basic();
    void generate_report(const std::string & outputPrefix, bool & is_pump_backflow,
                         bool & is_av_always_closed, bool & is_mv_always_opened,
                         double & sv_av, double & sv_pump, double & aortic_valve_opened_time,
                         double & mitral_valve_opened_time, double & pressure_pulsality,
                         double & max_pressure, double & min_pressure,
                         double & min_pump_flow, double & lv_work);

    void generate_report_edge_stats();
    void sensitivityAnalysis();
    void generate_report_wk_flows();
    void dump_data();

    std::string base_path;
    void set_base_path(std::string base_path) {
        this -> base_path.assign(base_path);
    }
    std::ofstream set_path_to_dump(const std::string &edge_id) {
       std::ofstream fout;
       fout.open(base_path + "/data/out/" + edge_id, std::ios::app);
       return fout;
    }
    void set_path_to_brachial_data(const std::string &base_path) {
        path_to_brachial_data = base_path + "/data/out/brachial.csv";
        fout_brachial.open(path_to_brachial_data);
    }

};
#endif
