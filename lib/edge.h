#ifndef EDGE_H
#define EDGE_H
#include <iostream>
#include <cmath>
#include <cassert>
#include <omp.h>
#include <algorithm>
#include "graph.h"
#include "calculator.h"
#include <fstream>

//important!
//Vertex doesn't know anything about its position on its edge
//We suppose its relative position is in the beginning
//Only its edge knows the absolute position of the both vertices
//So the edge manages everything
//Consider it as a change of variable


//right now we dont really need this class
class Meta_edge:
        public Simple_edge
{
    /*!
     * \brief Child class storing all common things for calculations at each time step.
     *        Right now not really used a lot.
     */
public:
    Meta_edge(const std::string & name, Simple_vertex *first, Simple_vertex *last);
    //virtual void solve() = 0;
    virtual ~Meta_edge() = 0;
};


template< template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Edge:
    public Meta_edge
{
    /*!
     * \brief This class implements explicit numerical scheme for
     *        1D hemodynamics. However it's an abstact class requiring
     *        methods with an equation of state (tube law). You need to specify them in child class.
     *
     */

protected:
    std::vector<std::reference_wrapper<BaseTimeCalculator>> timeCalculators;

    BinaryCalculator max_center_speed, min_center_speed;


    enum Vector_index {i_Sq, i_U, N, p_Sq = N, p_U = N + 1};//square of cross-section, velocity, and same from previous time step as well
    typedef Matrix<Scalar, N+N, Dynamic> GrandMatrix2ND;
    //typedef Matrix<Scalar, N, Dynamic> MatrixND;
    typedef Matrix<Scalar, N, N> MatrixNN;
    typedef Matrix<Scalar, N, 1> VectorN;
    typedef Matrix<Scalar, 1, Dynamic> RowVectorD;

    double len;
    double density;
    double viscosity;
    double zeta;
    double dx; //uniform step
    int points_num; //uniform partition, number of internal points on edge
    double gamma; // gamma = 0 for 2nd order scheme, gamma = 1 for 1st order monotonic scheme

    //matrix of size 2N x points_num
    GrandMatrix2ND GV; //matrix of column vectors, see Vector_index
    double mean_pressure_whole_vessel_at_t; //at one moment of time
    IntegralCalculator mean_pressure; //over cycle and the whole vessel
    IntegralCalculator mean_center_flow;
    IntegralCalculator mean_center_speed;

    /**
     * @brief Function to find pressure given area and position
     *
     * Usually it is specified by some tube law aka state equation
     *
     * @param S vascular area at a point
     * @param coord coordinate of a point
     * @return pressure at a point
     */
    virtual double get_pressure(const double & S, double coord) const = 0;

    /**
     * @brief get_pressure
     * @param S vascular area at a point
     * @param i node index of a point
     * @return pressure at a point
     */
    double get_pressure(const double & S, int i) const
    {
        // std::cout << "get_pressure(S, index) : " << S << std::endl;
        return get_pressure(S, i * dx);
    }

    /**
     * @brief get_pressure version with extra variables to find
     * @param S
     * @param speed_of_sound
     * @param P
     * @param dPdS
     * @param coord
     */
    virtual void get_pressure(const double & S, double & speed_of_sound, double & P, double & dPdS, double coord) const = 0;

    void get_pressure(const double & S, double & speed_of_sound, double & P, double & dPdS, int i) const
    {
        return get_pressure(S, speed_of_sound, P, dPdS, i * dx);
    }

    /**
     * @brief get_d_pressure_d_s
     * @param S
     * @param coord
     * @return
     */
    virtual double get_d_pressure_d_s(const double & S, double coord) const = 0;

    double get_d_pressure_d_s(const double & S, int i) const
    {
        return get_d_pressure_d_s(S, i * dx);
    }

    /**
     * @brief dPdx
     * @param S
     * @param coord
     * @return
     */
    virtual double dPdx(double S, double coord) const = 0;

    virtual void autoregulation()
    {
        return;
    }
public:
    int test_connection(Simple_vertex * v)
    {
        if (v == first_vertex || v == last_vertex)
            return 1;
        else
            return 0;
    }

    double get_mean_pressure() const
    {
        return mean_pressure.value();
    }
    double get_mean_pressure_whole_vessel_at_t() const
    {
        return mean_pressure_whole_vessel_at_t;
    }

    double get_max_center_speed()
    {
        return max_center_speed.value();
    }
    double get_min_center_speed()
    {
        return min_center_speed.value();
    }
    double get_center_speed() const
    {
        return get_U_unsafe(points_num / 2);
    }
    double get_center_pressure()
    {
        return get_P_unsafe(points_num / 2);
    }
    double get_center_area() const
    {
        return get_S_unsafe(points_num / 2);
    }
    double get_center_flow() const
    {
        return get_center_speed() * get_center_area();
    }
    double get_mean_center_flow() const
    {
        return mean_center_flow.value();
    }
    double get_mean_center_speed() const
    {
        return mean_center_speed.value();
    }

    virtual void print_info(std::ostream & os)
    {}

    virtual void print_info()
    {}

    const double * get_data(quantity vt, int & size, int & step) override
    {
        size = get_points_num();
        double * data = get_V_data();
        step = 2*N;
        switch (vt) {

        case area: {
            return data + i_Sq;
        }

        case speed: {
            return data + i_U;
        }

        default:
            break;
        }
        return nullptr;
    }

    double get_pressure(const double & S, const Simple_vertex * vertex, int i) const
    {
        if (vertex == first_vertex) {
            return get_pressure(S, i);
        } else {
            assert(vertex == last_vertex);
            return get_pressure(S, points_num - 1 - i);
        }
    }
    void get_pressure(const double & S, double & speed_of_sound, double & P, double & Pr, const Simple_vertex * vertex, int i) const
    {
        if (vertex == first_vertex) {
            return get_pressure(S, speed_of_sound, P, Pr, i);
        } else {
            assert(vertex == last_vertex);
            return get_pressure(S, speed_of_sound, P, Pr, points_num - 1 - i);
        }
    }

    double get_d_pressure_d_s(const double & S, const Simple_vertex * vertex, int i) const
    {
        if (vertex == first_vertex) {
            return get_d_pressure_d_s(S, i);
        } else {
            assert(vertex == last_vertex);
            return get_d_pressure_d_s(S, points_num - 1 - i);
        }
    }

    double get_P_unsafe(int i) override
    {
        return get_pressure(GV(i_Sq, i), i);
    }
    double get_Q_unsafe(int i) const override
    {
        return GV(i_Sq, i) * GV(i_U, i);
    }
    double get_S_unsafe(int i) const override
    {
        return GV(i_Sq, i);
    }
    double get_U_unsafe(int i) const override
    {
        return GV(i_U, i);
    }

    Edge(const std::string & name, Simple_vertex *first, Simple_vertex *last, double dx_step, double len, double density, double viscosity, double zeta, double gamma)
    : Meta_edge(name, first, last),
      max_center_speed(static_cast<double const&(*)(double const&, double const&)>(std::max), std::bind(&Edge::get_center_speed, this), -99999, 1),
      min_center_speed(static_cast<double const&(*)(double const&, double const&)>(std::min), std::bind(&Edge::get_center_speed, this), 99999, 1),
      len(len), density(density), viscosity(viscosity), zeta(zeta), gamma(gamma),
      mean_pressure(std::bind(&Edge::get_mean_pressure_whole_vessel_at_t, this), 1),
      mean_center_flow(std::bind(&Edge::get_center_flow, this), 1),
      mean_center_speed(std::bind(&Edge::get_center_speed, this), 1)

    {
        timeCalculators.push_back(max_center_speed);
        timeCalculators.push_back(min_center_speed);
        timeCalculators.push_back(mean_pressure);
        timeCalculators.push_back(mean_center_flow);
        timeCalculators.push_back(mean_center_speed);
        points_num = std::max(static_cast<int>(ceil(len / dx_step)), 3);
        dx = len / (points_num - 1);
        GV.resize(2*N, points_num);
    }

    double dt; // time step
    double T;
    double max_abs_eigenvalue_div_dx; // this is to find next dt according to Courant condition

    double get_max_abs_eigenvalue_div_dx() const
    {
        return max_abs_eigenvalue_div_dx;
    }
    double get_viscosity() const
    {
        return viscosity;
    }
    double get_density() const
    {
        return density;
    }
private:
    double *get_V_data()
    {
        return GV.data();
    }
public:
    int get_points_num() const override
    {
        return points_num;
    }

    // TODO : remove it, configure printing from the task
    double Period;
    void set_period(double new_period) {
        Period = new_period;
        for (BaseTimeCalculator & calc: timeCalculators)
            calc.set_period(new_period);
    }

private:
    /**
     * @brief right_hand_side with a force of friction
     * @param Vi (u,A) vector
     * @param i mesh node index
     * @return right hand side vector
     */
    VectorN right_hand_side(const VectorN & Vi, int i) const
    {
        VectorN res;
        res(0) = 0;
        res(1) = -2 * (zeta + 2) * viscosity * Vi(i_U) * M_PI / density / Vi(i_Sq);
        return res;
    }
    /**
     * @brief FV
     *
     * F term
     *
     * @param Vi (u, S) vector
     * @param P pressure
     * @return
     */
    VectorN FV(const VectorN & Vi, const double & P) const
    {
        VectorN F;
        F(0) = Vi(i_Sq) * Vi(i_U);
        F(1) = pow(Vi(i_U), 2) / 2 + P / density;
        return F;
    }
public:
    void set_dt(double value)
    {
        dt = value;
    }
    void set_T(double value)
    {
        T = value;
    }
    double get_len() const
    {
        return len;
    }


    VectorN get_right_hand_side(const Simple_vertex *vertex, const VectorN & D, int i) const
    {
        if (vertex == first_vertex) {
            return right_hand_side(D, i)                  - VectorN(0, dPdx(D(i_Sq),                   i*dx) / density);
        } else {
            assert(vertex == last_vertex);
            return right_hand_side(D, points_num - 1 - i) + VectorN(0, dPdx(D(i_Sq), (points_num - 1 - i)*dx) / density);
        }
    }
private:
    void Eigenvalues_Hybr(const VectorN & Vi, VectorN & eigenvalues,
                            const double & speed_of_sound, const double & P, MatrixNN & B_h,
                            MatrixNN & D_h) const
    {
        eigenvalues(0) = Vi(i_U) - speed_of_sound;
        eigenvalues(1) = Vi(i_U) + speed_of_sound;
	
        const VectorN sigma = eigenvalues * dt / dx;

        B_h(0,0) = fabs(sigma(0)) * (1.0 + 5.0 * (1.0 - gamma) * (1.0 - fabs(sigma(0))) / 19.0) / 2.0;
    	B_h(0,1) = 0;
        B_h(1,1) = fabs(sigma(1)) * (1.0 + 5.0 * (1.0 - gamma) * (1.0 - fabs(sigma(1))) / 19.0) / 2.0;
    	B_h(1,0) = 0;
	
        D_h(0,0) = 6.0 * (1.0 - gamma) * sigma(0) * (1.0 / fabs(sigma(0)) - 1.0) / 19.0;
    	D_h(0,1) = 0;
        D_h(1,1) = 6.0 * (1.0 - gamma) * sigma(1) * (1.0 / fabs(sigma(1)) - 1.0) / 19.0;
    	D_h(1,0) = 0;
    }
public:
    /**
     * @brief get_big_omega get big omega matrix, and a reverse to it
     * @param Vi
     * @param big_omega
     * @param inverse_big_omega
     * @param speed_of_sound
     * @param P
     */
    void get_big_omega(const VectorN & Vi, MatrixNN & big_omega, MatrixNN & inverse_big_omega,
                const double & speed_of_sound, const double & P) const
    {
        const double S = Vi(i_Sq);

        big_omega(0,0) = speed_of_sound/S;
        big_omega(0,1) = - 1;
        big_omega(1,0) = speed_of_sound/S;
        big_omega(1,1) = 1;

        inverse_big_omega(0,0) = S;
        inverse_big_omega(0,1) = S;
        inverse_big_omega(1,0) = - speed_of_sound;
        inverse_big_omega(1,1) = speed_of_sound;
        inverse_big_omega /= (2.0 * speed_of_sound);
    }

    void get_eigenvalues(const VectorN & Vi, VectorN & eigenvalues, MatrixNN & SL,
                     MatrixNN & SP, const double & speed_of_sound, const double & P) const
    {
        eigenvalues(0) = Vi(i_U) - speed_of_sound;
        eigenvalues(1) = Vi(i_U) + speed_of_sound;

        SL(0,0) = fabs(eigenvalues(0));
    	SL(0,1) = 0;
        SL(1,1) = fabs(eigenvalues(1));
    	SL(1,0) = 0;

        SP(0,0) = fabs(eigenvalues(0)) / eigenvalues(0);
    	SP(0,1) = 0;
        SP(1,1) = fabs(eigenvalues(1)) / eigenvalues(1);
    	SP(1,0) = 0;
    }
private:
    VectorN SH_Hybrid2nd(const VectorN & B_omV2, const VectorN & D_omV2, const VectorN & Vi, const double & coord)
    {
        double speed_of_sound, P, dPdS;
        VectorN eigenvalues;
        MatrixNN big_omega, inverse_big_omega, B_h, D_h;

        get_pressure(Vi(i_Sq), speed_of_sound, P, dPdS, coord);

        get_big_omega(Vi, big_omega, inverse_big_omega, speed_of_sound, P);

        Eigenvalues_Hybr(Vi, eigenvalues, speed_of_sound, P, B_h, D_h);
        const double m_ev1 = fabs(eigenvalues(0)),
                     m_ev2 = fabs(eigenvalues(1));
        max_abs_eigenvalue_div_dx = std::max(max_abs_eigenvalue_div_dx, std::max(m_ev1, m_ev2) / dx);
        return inverse_big_omega * (B_h * big_omega * B_omV2 + D_h * big_omega * D_omV2);
    }


private:
    VectorN get_Vint(const VectorN & Vm1, const VectorN & Vp1,
                     const double & Pm1, const double & Pp1)
    {
        return - (dt / (2*dx)) * (FV(Vp1, Pp1) - FV(Vm1, Pm1));
    }
    VectorN get_VintLeft(const VectorN & V, const VectorN & Vp1, const VectorN & Vp2,
                     const double & P, const double &Pp1, const double & Pp2)
    {
        return - (dt / (2*dx)) * (-3 * FV(V, P) + 4 * FV(Vp1, Pp1) - FV(Vp2, Pp2));
    }
    VectorN get_VintRight(const VectorN & Vm2, const VectorN & Vm1, const VectorN & V,
                     const double & Pm2, const double &Pm1, const double & P)
    {
        return - (dt / (2*dx)) * (3 * FV(V, P) - 4 * FV(Vm1, Pm1) + FV(Vm2, Pm2));
    }
public:
    /**
     * @brief test_80mmHg_L2_error
     *
     * get L2 error during test when input to the net is a constant pressure 80 mmHg
     */
    void test_80mmHg_L2_error()
    {
        const double expected_pressure = 80;
        double err = 0;
        for (int i = 0; i < points_num; i++) {
            err += pow(GV(i_U, i)* GV(i_Sq, i) - expected_pressure, 2) * dx;
        }
        std::cout << "L2-error: " << sqrt(err) << std::endl;
    }
    double solve()
    {
        autoregulation();
        mean_pressure_whole_vessel_at_t = 0;
        max_abs_eigenvalue_div_dx = 0;

        int i = 0;
        // save to previous time step
        // this scheme actually does not need previous time step values at all
        // but we still need them to find boundary conditions in nodes
        GV. template bottomRows<N>() = GV. template topRows<N>();

        double pressure_at_i        = get_pressure(GV(i_Sq, i    ), i    );
        double pressure_at_i_plus_1 = get_pressure(GV(i_Sq, i + 1), i + 1);
        double pressure_at_i_plus_2 = get_pressure(GV(i_Sq, i + 2), i + 2);

        mean_pressure_whole_vessel_at_t += pressure_at_i;

        VectorN VintI;
        VectorN VintIP1;

        //VintI.setZero();
        VintI = get_VintLeft(GV. template block<N,1>(0, i), GV. template block<N,1>(0, i + 1),
                             GV. template block<N,1>(0, i + 2), pressure_at_i, pressure_at_i_plus_1, pressure_at_i_plus_2);
        VintIP1 = get_Vint(GV. template block<N,1>(0, i), GV. template block<N,1>(0, i + 2), pressure_at_i, pressure_at_i_plus_2);

        VectorN B_omV2 = GV. template block<N,1>(0, i + 1) - GV. template block<N,1>(0, i);
        VectorN D_omV2 = VintI + VintIP1;

        VectorN V_at_i_plus_05 = (GV. template block<N,1>(0, i) + GV. template block<N,1>(0, i + 1)) / 2;

        VectorN BD_Sum1;
        VectorN BD_Sum2 = SH_Hybrid2nd(B_omV2, D_omV2, V_at_i_plus_05, dx * (i + 0.5));

        //V at i = 0 stays the same
        //GV.block<N,1>(0, i) = GV.block<N,1>(0, i)

        //now loop
        for (i = 1; i < points_num - 2; i++) {
            pressure_at_i        = pressure_at_i_plus_1;
            pressure_at_i_plus_1 = pressure_at_i_plus_2;
            pressure_at_i_plus_2 = get_pressure(GV(i_Sq, i + 2), i + 2);

            mean_pressure_whole_vessel_at_t += pressure_at_i;

            VintI = VintIP1;
            VintIP1 = get_Vint(GV.template block<N,1>(0, i), GV. template block<N,1>(0, i + 2), pressure_at_i, pressure_at_i_plus_2);

            //now a bit more complex procedure

            B_omV2 = GV. template block<N,1>(0, i + 1) - GV. template block<N,1>(0, i);
            D_omV2 = VintI + VintIP1;

            V_at_i_plus_05 = (GV. template block<N,1>(0, i) + GV. template block<N,1>(0, i + 1)) / 2;
            BD_Sum1 = BD_Sum2;
            BD_Sum2 = SH_Hybrid2nd(B_omV2, D_omV2, V_at_i_plus_05, dx * (i + 0.5));

            GV. template block<N,1>(0, i) += dt * right_hand_side(GV. template block<N,1>(0, i), i)
                    - BD_Sum1 + VintI + BD_Sum2;
        }
        //now i = points_num - 2

        double pressure_at_i_minus_1 = pressure_at_i;
        pressure_at_i                = pressure_at_i_plus_1;
        pressure_at_i_plus_1         = pressure_at_i_plus_2;

        mean_pressure_whole_vessel_at_t += pressure_at_i_plus_1;

        //VintIM1 = VintI;
        VintI = VintIP1;
        //VintIP1.setZero();
        VintIP1 = get_VintRight(GV. template block<N,1>(N, i - 1), GV. template block<N,1>(0, i),
                               GV. template block<N,1>(0, i + 1), pressure_at_i_minus_1, pressure_at_i, pressure_at_i_plus_1);

        //now a bit more complex procedure

        B_omV2 = GV. template block<N,1>(0, i + 1) - GV. template block<N,1>(0, i);
        D_omV2 = VintI + VintIP1;

        V_at_i_plus_05 = (GV. template block<N,1>(0, i) + GV. template block<N,1>(0, i + 1)) / 2;
        BD_Sum1 = BD_Sum2;
        BD_Sum2 = SH_Hybrid2nd(B_omV2, D_omV2, V_at_i_plus_05, dx * (i + 0.5));

        GV. template block<N,1>(0, i) += dt * right_hand_side(GV. template block<N,1>(0, i), i)
                - BD_Sum1 + VintI + BD_Sum2;

        //V at i = points_num - 1 stays the same

        mean_pressure_whole_vessel_at_t /= points_num;

        for (BaseTimeCalculator & calc: timeCalculators)
            calc.update(dt);

        return max_abs_eigenvalue_div_dx;
    }

    double get_dx() const
    {
       return dx;
    }

    /**
     * @brief get_V_s get current vascular area value at i-th node wrt to vertex position
     * @param v vertex
     * @param i node index
     * @return area
     */
    double get_V_s(const Simple_vertex * v, const int & i) const 
    {
        if (v == first_vertex) {
            return GV(i_Sq, i);
        } else {
            assert(v == last_vertex);
            return GV(i_Sq, points_num - i - 1);
        }
    }
    /**
     * @brief set_V_s set current vascular area value at i-th node wrt to vertex position
     * @param v vertex
     * @param i node index
     * @param value area value
     */
    void set_V_s(const Simple_vertex * v, const int & i, const double & value)
    {
        if (v == first_vertex) {
            GV(i_Sq, i) = value;
        } else {
            assert(v == last_vertex);
            GV(i_Sq, points_num - i - 1) = value;
        }
    }
    /**
     * @brief get_V_u get current speed value at i-th node wrt to vertex position
     * @param v vertex
     * @param i node index
     * @return speed
     */
    double get_V_u(const Simple_vertex * v, const int & i) const
    {
        if (v == first_vertex) {
            return GV(i_U, i);
        } else {
            assert(v == last_vertex);
            return -GV(i_U, points_num - i - 1);
        }
    }
    /**
     * @brief set_V_u set current speed value at i-th node wrt to vertex position
     * @param v vertex
     * @param i node index
     * @param value speed value
     */
    void set_V_u(const Simple_vertex * v, const int & i, const double & value)
    {
        if (v == first_vertex) {
            GV(i_U, i) = value;
        } else {
            assert(v == last_vertex);
            GV(i_U, points_num - i - 1) = -value;
        }
    }
    /**
     * @brief get_Vp_s get previous area value at i-th node wrt to vertex position
     * @param v vertex
     * @param i node index
     * @return area
     */
    double get_Vp_s(const Simple_vertex * v, const int & i) const 
    {
        if (v == first_vertex) {
            return GV(p_Sq, i);
        } else {
            assert(v == last_vertex);
            return GV(p_Sq, points_num - i - 1);
        }
    }
    /* you dont need it
    void set_Vp_s(const Simple_vertex * v, const int & i, const double & value)
    {
        if (v == first_vertex) {
            Vp(0, i) = value;
        } else {
            assert(v == last_vertex);
            Vp(0, points_num - i - 1) = value;
        }
    }
    */
    /**
     * @brief get_Vp_u get previous speed value at i-th node wrt to vertex position
     * @param v vertex
     * @param i node index
     * @return speed
     */
    double get_Vp_u(const Simple_vertex * v, const int & i) const
    {
        if (v == first_vertex) {
            return GV(p_U, i);
        } else {
            assert(v == last_vertex);
            return -GV(p_U, points_num - i - 1);
        }
    }
    /* you dont need it
    void set_Vp_u(const Simple_vertex * v, const int & i, const double & value)
    {
        if (v == first_vertex) {
            Vp(1, i) = value;
        } else {
            assert(v == last_vertex);
            Vp(1, points_num - i - 1) = -value;
        }
    }
    */
};


/**
 * @brief Edge class with a classic tube law
 */
template< template<typename, int, int...> typename Matrix, typename Scalar, auto Dynamic>
class EdgeClassicTubeLaw:
        public Edge<Matrix, Scalar, Dynamic>
{
    typedef Edge<Matrix, Scalar, Dynamic> Base;
    using Base::i_Sq;
    using Base::i_U;
    using Base::p_Sq;
    using Base::p_U;
    using Base::GV;
    using Base::len;
    using Base::density;
    using Base::points_num;
    using Base::T;
    using typename Base::RowVectorD;
    using typename Base::VectorN;

    const double bar_to_mmhg = 0.0007500637554192;
    int point_mean_ind;
    double width;
    double c; //velocity of small disturbances propagation in the vessel wall, characterizes vessel elasticity
    double S0; //cross-sectional area
    double c0; //for autoregulation
    double p0; //for autoregulation
    int first_time = 1;
    int autoregulation_enabled;
    double get_S0dx(double coord) const
    {
        double S1 = S0, S2 = S1 + 0;
        return (S2 - S1) / len;
    }
    double get_S0(double coord) const
    {
        double S1 = S0, S2 = S1 + 0;
        return S1 + (S2 - S1) * coord / len;
    }

    // std::ofstream fout;
    int k = -1; // for proper output in files

protected:
    virtual void autoregulation() final
    {
        /*
        if (T > 8) {
            //turn on autoregulation
            if (first_time) {
                //save p0
                first_time = 0;
                p0 = Base::get_mean_pressure();
            }
            c = c0 * sqrt(Base::get_mean_pressure() / p0);
        }
        */
        if (T > 3 && autoregulation_enabled)
            c = c0 * sqrt(Base::get_mean_pressure() / p0);
    }
public:
    /*
    void print_info(std::ostream & os) const
    {
        Base::print_info(os);
        DUMP_TO_OS(c);
        DUMP_TO_OS(S0);
        DUMP_TO_OS(c0);
    }
    */
    virtual void print_info(std::ostream &os) final {
        if (T >= 9.0 * this -> Period && T <= 10.0 * this -> Period) {
            os << T << "," << this -> get_center_pressure() * bar_to_mmhg <<
                "," << this -> get_center_speed() <<
                "," << this -> get_center_area() <<
                "," << this -> get_center_flow() <<
                "," << this -> get_mean_center_flow() << std::endl;
        }
    }

    EdgeClassicTubeLaw(const std::string & name, Simple_vertex *first, Simple_vertex *last,
                       double dx_step, double len, double density, double viscosity, double zeta, double gamma,
                   double width, double c, double p0_val)
    : Base(name, first, last, dx_step, len, density, viscosity, zeta, gamma) , width(width), c(c), c0(c)
    {
        S0 = M_PI * width * width / 4;
        GV.row(p_Sq) = GV.row(i_Sq) = RowVectorD::Constant(points_num, S0); //initial S
        GV.row(p_U) = GV.row(i_U) = RowVectorD::Constant(points_num, 0); //initial U
        if (p0_val == 0) {
            autoregulation_enabled = 0;
            // std::cout << "Disabled" << std::endl;
        } else {
            autoregulation_enabled = 1;
        }
        autoregulation_enabled = 0;
        p0 = p0_val;
        point_mean_ind = points_num / 2;
/*
        std::string path_to_file("/home/artem/Documents/work/sechenov/bloodflow-main/out/data/");
        path_to_file = path_to_file + name + ".csv";
        fout.open(path_to_file);
*/
    }

protected:
        double dPdx(double S, double coord) const
        {
            const double SS0 = get_S0(coord);
            const double eta = S / SS0;
            const double deta = - S / (SS0 * SS0) * get_S0dx(coord);
            if (eta > 1) {
                return density * c*c * std::exp(eta - 1) * deta;
            } else {
                return density * c*c / eta * deta;
            }
        }
        void get_pressure(const double & S, double & speed_of_sound, double & P, double & Pr, double coord) const
        {
            const double eta = S/get_S0(coord);
            if (eta > 1) {
                P = density * c*c * std::expm1(eta - 1);
                Pr = density * c*c * std::exp(eta - 1) / get_S0(coord);
            } else {
                P = density * c*c * std::log(eta);
                Pr = density * c*c / S;
            }
            speed_of_sound = std::sqrt(S * Pr / density);
    #ifndef NDEBUG
            if (speed_of_sound < 1e-8)
                std::cout << "Sound is too low!" << std::endl;
    #endif
        }
        double get_pressure(const double & S, double coord) const
        {
            const double eta = S/get_S0(coord);
            if (eta > 1) {
                return density * c*c * std::expm1(eta - 1);
            } else {
                return density * c*c * std::log(eta);
            }
        }
        double get_d_pressure_d_s(const double & S, double coord) const
        {
            const double eta = S/get_S0(coord);
            if (eta > 1) {
                return density * c*c * std::exp(eta - 1) / get_S0(coord);
            } else {
                return density * c*c / S;
            }
        }
};


/**
 * @brief Edge class with a tube law used in Adan56 model
 */
/*
template< template<typename, int, int...> typename Matrix, typename Scalar, auto Dynamic>
class EdgeAdan56TubeLaw:
        public Edge<Matrix, Scalar, Dynamic>
{
    typedef Edge<Matrix, Scalar, Dynamic> Base;
    using Base::i_Sq;
    using Base::i_U;
    using Base::p_Sq;
    using Base::p_U;
    using Base::GV;
    using Base::density;
    using Base::len;
    using Base::points_num;
    using Base::dx;
    using typename Base::Vector_index;
    using typename Base::VectorN;
    using typename Base::RowVectorD;

    double P_exterior, P_diastolic;
    double ca, cb, cc, cd;
    double E;//Young moduli

    double R1, R2; //cross-sectional radii (proximal and distal) to heart
    //R1 for first vertex, R2 for last

    double r0(double coord) const
    {
        return R1 + (R2 - R1) * coord / len;
    }
    double dr0dx(double coord) const
    {
        return (R2 - R1) / len;
    }
    double h(double r) const
    {//wall thickness
        double res = r * (ca * std::exp(cb * r) + cc * std::exp(cd * r));
        return res;
    }
    double dhdr(double r) const
    {
        double res = (r * cb  + 1) * ca * std::exp(cb * r) + (r * cd + 1) * cc * std::exp(cd * r);
        return res;
    }
    double betaC(double r) const
    {
        return 4.0/3 * std::sqrt(M_PI) * E * h(r);
    }
    double betaCdr(double r) const
    {
        return 4.0/3 * std::sqrt(M_PI) * E * dhdr(r);
    }
public:
    double mean_diameter() const
    {
        return R1 + R2;
    }
    void print_info(std::ostream & os) const
    {
        Base::print_info(os);
        DUMP_TO_OS(P_exterior);
        DUMP_TO_OS(P_diastolic);
        DUMP_TO_OS(ca);
        DUMP_TO_OS(cb);
        DUMP_TO_OS(cc);
        DUMP_TO_OS(cd);
        DUMP_TO_OS(E);
        DUMP_TO_OS(R1);
        DUMP_TO_OS(R2);
    }
    EdgeAdan56TubeLaw(const std::string & name, Simple_vertex *first,
                      Simple_vertex *last, double dx_step, double len,
                      double density, double viscosity, double zeta, double gamma,
                   double P_exterior, double P_diastolic, double ca,
                   double cb, double cc, double cd, double E, double R1, double R2)
        : Base(name, first, last, dx_step, len, density, viscosity, zeta, gamma),
          P_exterior(P_exterior), P_diastolic(P_diastolic),
          ca(ca), cb(cb), cc(cc), cd(cd), E(E), R1(R1), R2(R2)
    {
        for (int i = 0; i < points_num; i++) {
            const double rr = r0(i*dx);
            GV(p_Sq, i) = GV(i_Sq, i) = M_PI * rr*rr;
        }
        GV.row(p_U) = GV.row(i_U) = RowVectorD::Constant(points_num, 0); //initial U
    }

protected:
    double dPdx(double S, double coord) const
    {
        const double rr0 = r0(coord);
        const double r0dx = dr0dx(coord);
        const double A0 = M_PI*rr0*rr0;
        const double dA0dr = 2*M_PI*rr0;
        const double betaDivA0 = betaC(rr0) / A0;
        const double betaDivA0dr = (betaCdr(rr0) - betaDivA0 * dA0dr) / A0;

        return betaDivA0dr * r0dx * (std::sqrt(S) - std::sqrt(M_PI) * rr0) +
                betaDivA0 * ( - std::sqrt(M_PI) * r0dx);
    }
    void get_pressure(const double & S, double & speed_of_sound, double & P, double & Pr, double coord) const
    {
        const double rr0 = r0(coord);
        const double A0 = M_PI*rr0*rr0;

        const double betaDivA0 = betaC(rr0) / A0;
        P = P_exterior + P_diastolic + betaDivA0 * (std::sqrt(S) - std::sqrt(M_PI)*rr0);
        Pr = betaDivA0 / (2 * std::sqrt(S));
        speed_of_sound = std::sqrt(S * Pr / density);
        #ifndef NDEBUG
            if (speed_of_sound < 1e-8)
                std::cout << "Sound is too low!" << std::endl;
        #endif
    }
    double get_pressure(const double & S, double coord) const
    {
        const double rr0 = r0(coord);
        const double A0 = M_PI * rr0 * rr0;
        return P_exterior + P_diastolic + betaC(rr0) / A0 * (std::sqrt(S) - std::sqrt(M_PI)*rr0);
    }
    double get_d_pressure_d_s(const double & S, double coord) const
    {
        const double rr0 = r0(coord);
        const double A0 = M_PI * rr0 * rr0;
        return betaC(rr0) / (2 * A0 * std::sqrt(S));
    }
};
*/

#endif
