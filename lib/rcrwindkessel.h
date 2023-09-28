#ifndef RCRWINDKESSEL_H
#define RCRWINDKESSEL_H

#include <iostream>
#include <cmath>

#include "vertex.h"

template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class TerminalVertex:
    public Vertex<Edge, Matrix, Scalar, Dynamic>
{
    typedef Vertex<Edge, Matrix, Scalar, Dynamic> Base;
    IntegralCalculator strokeVOL;
    using Base::timeCalculators;
    virtual double get_flow() = 0;

public:
    TerminalVertex(const std::string & id)
        : Base(id),
          strokeVOL(std::bind(&TerminalVertex::get_flow, this), 1)
    {
        timeCalculators.push_back(strokeVOL);
    }
    double get_stroke_volume()
    {
        return strokeVOL.value();
    }
    virtual void print_info() {}
};

template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class RCR_Windkessel_vertex:
    public TerminalVertex<Edge, Matrix, Scalar, Dynamic>
{
    typedef TerminalVertex<Edge, Matrix, Scalar, Dynamic> Base;
    using Base::OutgoingCompatibilityCoeffs;
    using Base::dt;
    using Base::OCNFO2;
    using Base::id;
    // seems like it's edge before wk => get_Q_in might be correct
    Edge *e;
    Simple_vertex * sv;
    /*what is it for? i dont remember
    void strangeFunc(double &S, double &U)
    {
        const double Pout = 6666.12; //5 mmhg
        const double R = 800;

        const int Nmax = 30;
        const double eps = 1e-8;

        S = e->get_V_s(this, 1);
        double alpha, beta;
        OutgoingCompatibilityCoeffs(e, sv, alpha, beta);

        double S_norm = 1;

        for (int i = 0; i < Nmax && S_norm > eps; i++) {
            double F = - e->Equation_of_state_optimized(S) + Pout - R*S*(alpha*S + beta);
            S -= F / (- e->Equation_of_state_deriv(S) - R*(2*alpha*S + beta) );
            S_norm = fabs(F);
        }
        if (S_norm > eps) {
            std::cout << S_norm << std::endl;
        }
        U = alpha * S + beta;
    }
*/

    double Pout;
    double R1;
    double R2;
    double C;

    double PPG = 0.0;
    double previous_time;
    double Period;

    void WindKessel(double &S, double &U)
    {

        const int Nmax = 20;

        S = e->get_V_s(sv, 1);//initial guess
        const double Qp = /*!!*/ e->get_Vp_s(sv, 0) * e->get_Vp_u(sv, 0);

        const double Pp = e->get_pressure(e->get_Vp_s(sv, 0), sv, 0);
        double alpha, beta;
        OutgoingCompatibilityCoeffs(e, sv, alpha, beta);

        double S_norm = 1;

//        for (int i = 0; i < Nmax && S_norm > eps; i++) {
//            const double u = alpha * S + beta;
//           const double P = e->Equation_of_state_optimized(S);
//            const double G = dt * ( -/*!!*/ (S*u + Qp) / 2 * (R1 + R2) - ((P + Pp)/2 - Pout));
//           double F = -/*!!*/ R1*R2*C*(S*u - Qp) - (R2*C*(P - Pp) - G);
//            const double G_der = dt * ( -/*!!*/ 0.5*(R1 + R2)*(2*alpha*S + beta) - e->Equation_of_state_deriv(S)/2 );
//            const double dS = F / (-/*!!*/ R1*R2*C*(2*alpha*S + beta) - ( R2*C*e->Equation_of_state_deriv(S) - G_der) );
//            S -= dS;
//            S_norm = fabs(F);
//        }
        const double atol = 1e-7;
        for (int i = 0; i < Nmax && S_norm > atol; i++) {
            const double u = alpha * S + beta;
            const double P = e->get_pressure(S, sv, 0);
            const double G = dt * ( -/*!!*/ S*u * (R1 + R2) - (P - Pout));
            const double F = -/*!!*/ R1*R2*C*(S*u - Qp) - (R2*C*(P - Pp) - G);
            const double G_der = dt * ( -/*!!*/ (R1 + R2)*(2*alpha*S + beta) - e->get_d_pressure_d_s(S, sv, 0) );
            const double dS = F / (-/*!!*/ R1*R2*C*(2*alpha*S + beta) - ( R2*C*e->get_d_pressure_d_s(S, sv, 0) - G_der) );
            S -= dS;
            S_norm = fabs(F);
        }
        if (S_norm > atol) {
            std::cout << "windkessel failed" << std::endl;
            throw("Windkessel Newton failed");
        }
        U = alpha * S + beta;
    }
    double Q;
public:
    virtual double get_flow() final
    {
        return Q;
    }
    void set_R1(double new_R1)
    {
        R1 = new_R1;
    }
    void set_R2(double new_R2)
    {
        R2 = new_R2;
    }
    double get_R1() const
    {
        return R1;
    }
    double get_R2() const
    {
        return R2;
    }
    RCR_Windkessel_vertex(const std::string & id, Edge * e, Simple_vertex * sv, double Pout, double R1, double R2, double C)
    : Base(id), e(e), sv(sv), Pout(Pout), R1(R1), R2(R2), C(C)
    {
    }
    void update() final
    {
        double S, U;
        //strangeFunc(S, U); i dont remember what is it for
        WindKessel(S, U);
        Q = S * U;
        e->set_V_s(sv, 0, S);
        e->set_V_u(sv, 0, U);
#ifndef NDEBUG
        //std::cout << id << ": " << S << " " << U << " " << S*U << std::endl;
#endif
    }
    virtual void print_info(const double &virtual_time, std::ostream &os) final
    {
        if (virtual_time >= 9.0 * this -> Period && virtual_time <= 10.0 * this -> Period) {
            double Qout = get_Q_out_alt();
            double Qin = get_Q_in(); // ~ get_flow()
            PPG += (Qin - Qout) * (virtual_time - previous_time);
            previous_time = virtual_time;
            os << virtual_time << "," << Qin << "," << Qout << "," << PPG << std::endl;
        }
        // std::cout << "edgeName & R1 & R2 & C & Pout" << std::endl;
        // std::cout << e->get_id() << " & " << R1 << " & " << R2 << " & " << C << " & " << Pout << " \\" << std::endl;
    }
    double get_Q_in() {
        return e -> get_center_flow();
    }
    double get_Q_in_alt() {
        return e -> get_Vp_s(sv, 0) * e -> get_Vp_u(sv, 0);
    }
    double get_Q_out_alt() {
        return (e -> get_pressure(e->get_Vp_s(sv, 0), sv, 0) - Pout) / (R1 + R2);
    }
    void set_period(double new_period) {
        Period = new_period;
        previous_time = 9.0 * new_period;
    }
};


template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class TerminalResistanceCoronaryArteries:
    public TerminalVertex<Edge, Matrix, Scalar, Dynamic>
{
    typedef TerminalVertex<Edge, Matrix, Scalar, Dynamic> Base;
    using Base::OutgoingCompatibilityCoeffs;
    using Base::T;
    using Base::OCNFO2;
    using Base::id;
    using Base::get_id;
    Edge *e;
    Simple_vertex * sv;

    double Pout;
    double R, R_max;

    void newton(double &S, double &U)
    {
        const int Nmax = 20;

        S = e->get_Vp_s(sv, 0);//initial guess
        double alpha, beta;
        OutgoingCompatibilityCoeffs(e, sv, alpha, beta);

        double error = 1;

        const double atol = 1e-7;

        const double systole_time = 0.35;
        const double time = fmod(T, 1);
        const double resist = R + (R_max - R) * sin(M_PI * time / systole_time) * (time < systole_time);

        for (int i = 0; i < Nmax && error > atol; i++) {
            const double u = alpha * S + beta;
            const double P = e->get_pressure(S, sv, 0);
            const double F = u * S - (Pout - P) / resist;
            const double dFdS = 2 * alpha * S + beta + e->get_d_pressure_d_s(S, sv, 0) / resist;
            S -= F / dFdS;
            error = fabs(F);
            //std::cout << "iter: " << i << ", error: " << error << std::endl;
        }
        if (error > atol) {
            std::cout << "TerminalResistanceCoronaryArteries failed" << std::endl;
            std::cout << "Vertex: " << sv->get_id() << " " << get_id() << std::endl;
            std::cout << "Error: " << error << std::endl;
            throw("TerminalResistanceCoronaryArteries Newton failed");
        }
        U = alpha * S + beta;
    }
    double Q;
public:
    virtual double get_flow() final
    {
        return Q;
    }
    TerminalResistanceCoronaryArteries(const std::string & id, Edge * e, Simple_vertex * sv, double Pout, double Resistance)
    : Base(id), e(e), sv(sv), Pout(Pout)
    {
        R =  22950   *   Resistance / 22950;
        R_max = 3 * R;
    }
    void update() final
    {
        double S, U;
        newton(S, U);
        Q = S * U;
        e->set_V_s(sv, 0, S);
        e->set_V_u(sv, 0, U);
#ifndef NDEBUG
        //std::cout << id << ": " << S << " " << U << " " << S*U << std::endl;
#endif
    }
    virtual void print_info() final
    {
        std::cout << "edgeName & R & Rmax & Pout" << std::endl;
        std::cout << e->get_id() << " & " << R << " & " << R_max << " & " << Pout << " \\" << std::endl;
    }
};

#endif // RCRWINDKESSEL_H
