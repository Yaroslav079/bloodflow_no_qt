#ifndef HEART_ADVANCED_VALVES_H
#define HEART_ADVANCED_VALVES_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <Eigen/Dense>
#include <cmath>

#include "vertex.h"
#include "heart_part.h"


template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Heart_AdValves:
    public True_0d_heart<Edge, Matrix, Scalar, Dynamic>, public Heart_part
{
    typedef True_0d_heart<Edge, Matrix, Scalar, Dynamic> Base;
    using Base::OutgoingCompatibilityCoeffs;
    using Base::T;
    using Base::dt;
    using Base::id;

    typedef Eigen::Matrix<double, 15, 15> Matrix15d;
    typedef Eigen::Matrix<double, 11, 15> Matrix11x15;

    typedef Eigen::Matrix<double, 4, 1> Vector4d;
    typedef Eigen::Matrix<double, 15, 1> Vector15d;
    typedef Eigen::Matrix<double, 11, 1> Vector11d;
    typedef Eigen::Matrix<double, 4, 15> Matrix4x15;

    enum {v_v, v_v_d, v_a, v_a_d, ao_valve, ao_valve_d, mi_valve, mi_valve_d, q_av, q_mi, q_pu,
          a_aorta, p_vent, p_aort, p_atri};

    // y0 << initLVvolume, 0, 20, 0, 0.01, 0, 0.01, 0, 1, 1, 1, 7.9, 20, 80, P_pulmonary;

    //edge starting from this node
    Edge * e;

    //simple vertex of this node
    Simple_vertex * sv;

    //compatibility coeffs
    double alfa, beta;

    //right side of ODE system
    Vector11d F (const Vector15d & y, const double & t) {
        Vector11d d;

        d(v_v) = y(v_v_d);

        d(v_v_d) = (-R1 * y(p_vent) * y(v_v_d) - e1(t) * (y(v_v) - V1_0) + y(p_vent) )/ I1;

        d(v_a) = y(v_a_d);

        d(v_a_d) = (-R4 * y(p_atri) * y(v_a_d) - e4(t) * (y(v_a) - V4_0) + y(p_atri)) / I4;


        d(ao_valve) = y(ao_valve_d);

        d(ao_valve_d) = -Fr(y(ao_valve)) - Ff(y(ao_valve_d)) + Fp(y(p_vent), y(p_aort), y(ao_valve));

        d(mi_valve) = y(mi_valve_d);

        d(mi_valve_d) = -Fr(y(mi_valve)) - Ff(y(mi_valve_d)) + Fp(y(p_atri), y(p_vent), y(mi_valve));


        d(q_av) = (y(p_vent) - y(p_aort) - B_av(y(ao_valve), y(a_aorta)) * y(q_av) * fabs(y(q_av))) / L_av;


        d(q_mi) = (y(p_atri) - y(p_vent) - B_mi(y(mi_valve)) * y(q_mi) * fabs(y(q_mi))) / L_mi;

        d(q_pu) = (P_pulmonary - y(p_atri) - B_pu() * y(q_pu) * fabs(y(q_pu)) ) / L_pu;

        return d;
    }

    Matrix11x15 F_d(const Vector15d & y, const double & t) {
        Matrix11x15 d;
        d.setZero();

        d(v_v, v_v_d) = 1;

        d(v_v_d, p_vent) = (-R1 * y(v_v_d) + 1 )/ I1;
        d(v_v_d, v_v_d) = (-R1 * y(p_vent))/ I1;
        d(v_v_d, v_v) = (- e1(t) )/ I1;

        d(v_a, v_a_d) = 1;

        d(v_a_d, p_atri) = (-R4 * y(v_a_d) + 1) / I4;
        d(v_a_d, v_a_d) = (-R4 * y(p_atri) ) / I4;
        d(v_a_d, v_a) = ( - e4(t) ) / I4;

        d(ao_valve, ao_valve_d) = 1;

        d(ao_valve_d, ao_valve) = -Fr_d(y(ao_valve)) + Fp_dTet(y(p_vent), y(p_aort), y(ao_valve));
        d(ao_valve_d, ao_valve_d) = - Ff_tet_d(y(ao_valve_d));
        d(ao_valve_d, p_vent) = Fp_dPfrom(y(p_vent), y(p_aort), y(ao_valve));
        d(ao_valve_d, p_aort) = Fp_dPto(y(p_vent), y(p_aort), y(ao_valve));

        d(mi_valve, mi_valve_d) = 1;

        d(mi_valve_d, mi_valve) = -Fr_d(y(mi_valve)) + Fp_dTet(y(p_atri), y(p_vent), y(mi_valve));
        d(mi_valve_d, mi_valve_d) = - Ff_tet_d(y(mi_valve_d)) ;
        d(mi_valve_d, p_atri) = Fp_dPfrom(y(p_atri), y(p_vent), y(mi_valve));
        d(mi_valve_d, p_vent) = Fp_dPto(y(p_atri), y(p_vent), y(mi_valve));



        d(q_av, a_aorta) = L_av_inverse_d(y(a_aorta)) * (y(p_vent) - y(p_aort) - R_av(y(ao_valve)) * y(q_av)
                                         - B_av(y(ao_valve), y(a_aorta)) * y(q_av) * fabs(y(q_av)))

                           + L_av_inverse(y(a_aorta)) * (- B_av_dAorta_area(y(ao_valve), y(a_aorta)) * y(q_av) * fabs(y(q_av)));

        d(q_av, p_vent) = L_av_inverse(y(a_aorta));
        d(q_av, p_aort) = L_av_inverse(y(a_aorta)) * ( -1);
        d(q_av, ao_valve) = L_av_inverse(y(a_aorta)) * ( - R_av_d(y(ao_valve)) * y(q_av)
                                         - B_av_dAv_state(y(ao_valve), y(a_aorta)) * y(q_av) * fabs(y(q_av)));

        d(q_av, q_av) = L_av_inverse(y(a_aorta)) * (- R_av(y(ao_valve))
                                         - B_av(y(ao_valve), y(a_aorta)) * 2 * fabs(y(q_av)));



        d(q_mi, p_atri) = L_mi_inverse();
        d(q_mi, p_vent) = L_mi_inverse() * (-1);
        d(q_mi, mi_valve) = L_mi_inverse() * ( - R_av_d(y(mi_valve)) * y(q_mi)
                                              - B_mi_dMi_state(y(mi_valve)) * y(q_mi) * fabs(y(q_mi)));

        d(q_mi, q_mi) = L_mi_inverse() * ( - R_av(y(mi_valve))
                                          - B_mi(y(mi_valve)) * 2 * fabs(y(q_mi)));



        d(q_pu, p_atri) = L_pu_inverse() * (-1);
        d(q_pu, q_pu) = L_pu_inverse() * (- R_pu() - B_pu() * 2 * fabs(y(q_pu)) );

        return d;
    }


    //algebraic part
    Vector4d f(const Vector15d & y)
    {
        Vector4d ff;
        ff(0) = y(v_v_d) - y(q_mi) + y(q_av);
        ff(1) = y(v_a_d) - y(q_pu) + y(q_mi);
        ff(2) = y(p_aort) - e->get_pressure(y(a_aorta), sv, 0)/cc;
        ff(3) = y(q_av) - y(a_aorta) * (alfa * y(a_aorta) + beta);
        return ff;
    }
    Matrix4x15 f_d(const Vector15d & y)
    {
        Matrix4x15 ff;
        ff.setZero();

        ff(0, v_v_d) = 1;
        ff(0, q_mi) = -1;
        ff(0, q_av) = 1;

        ff(1, v_a_d) = 1;
        ff(1, q_pu) = -1;
        ff(1, q_mi) = 1;

        ff(2, p_aort) = 1;
        ff(2, a_aorta) = - e->get_d_pressure_d_s(y(a_aorta), sv, 0)/cc;

        ff(3, q_av) = 1;
        ff(3, a_aorta) = - (2 * alfa * y(a_aorta) + beta);

        return ff;
    }



    Vector15d y0, yn, R, y2prev;
    Matrix15d B;

    double LV_P, LV_V;
    double LA_P, LA_V;
    double aortic_valve, mitral_valve;
    double flow_val_previous = 0.0;
public:
    double get_ESPVR()
    {
        return Heart_part::get_ESPVR();
    }
    void set_ESPVR(double new_ESPVR)
    {
        Heart_part::set_ESPVR(new_ESPVR);
    }
    double get_PveinPressure()
    {
        return Heart_part::get_PveinPressure();
    }
    void set_PveinPressure(double new_PveinPressure)
    {
        Heart_part::set_PveinPressure(new_PveinPressure);
    }
    double get_LV_P()
    {
        return LV_P;
    }
    double get_LV_V()
    {
        return LV_V;
    }
    double get_LA_P()
    {
        return LA_P;
    }
    double get_LA_V()
    {
        return LA_V;
    }
    double get_aortic_valve()
    {
        return aortic_valve;
    }
    double get_mitral_valve()
    {
        return mitral_valve;
    }
    double get_flow_av()
    {
        return y0(q_av);
    }
    double get_aortic_root_pressure()
    {
        return y0(p_aort);
    }

    double get_HR() {
        return 60 / Period;
    }

    double get_Tsys() {
        return Ts1;
    }

    double get_flow_av_fwd() {
        if (y0(q_av) > 0.0)
            return y0(q_av);
        else
            return 0.0;
    }

    double get_flow_av_bwd() {
        if (y0(q_av) < 0.0)
            return y0(q_av);
        else
            return 0.0;
    }

    double get_flow_mv() {
        return y0(q_mi);
    }
    double get_flow_mv_fwd() {
        if (y0(q_mi) > 0.0)
            return y0(q_mi);
        else
            return 0.0;
    }
    double get_flow_mv_bwd() {
        if (y0(q_mi) < 0.0)
            return y0(q_mi);
        else
            return 0.0;
    }

    double get_flow_pv() {
        return y0(q_pu);
    }

    double get_lvet() {
        if (y0(q_av) > 0.0)
            return 1.0;
        else
            return 0.0;
    }

    double get_pft() {
        if (y0(q_av) > flow_val_previous) {
            flow_val_previous = y0(q_av);
            return 1.0;
        }
        else {
            flow_val_previous = y0(q_av);
            return 0.0;
        }
    }

    void set_params_av(const double &param_b_av, const double &param_l_av) {
        B_av_denomin = param_b_av;
        L_av = param_l_av;
    }

    void set_params_mi(const double &param_b_mi, const double &param_l_mi) {
        B_mi_denomin = param_b_mi;
        L_mi = param_l_mi;
    }

    void set_params_pu(const double &param_b_pu, const double &param_l_pu) {
        B_pu_const = param_b_pu;
        L_pu = param_l_pu;
    }

    void set_pulmVeinsPressure(const double &p) {
        P_pulmonary = p;
    }

    void set_I4(const double &param_I4) {
        I4 = param_I4;
    }

    void set_R4(const double &param_R4) {
        R4 = param_R4;
    }

    void set_I1(const double &param_I1) {
        I1 = param_I1;
    }

    void set_R1(const double &param_R1) {
        R1 = param_R1;
    }

    double get_param_b_av() {
        return B_av_denomin;
    }

    double get_param_l_av() {
        return L_av;
    }

    double get_param_b_mi() {
        return B_mi_denomin;
    }

    double get_param_l_mi() {
        return L_mi;
    }

    double get_param_b_pu() {
        return B_pu_const;
    }

    double get_param_l_pu() {
        return L_pu;
    }

    double get_pulmVeinsPressure() {
        return P_pulmonary;
    }

    double get_period() {
        return Period;
    }

    void print_info() {
        /*
        if ((T > 9.0 * Period) && (T < 10.0 * Period)) {
            os << T << "," << y0(q_av) << "," << y0(p_vent) << "," << y0(v_v) << std::endl;
        }
        */
    }

    void print_info(std::ostream &os) {
        if ((T > 9.0 * Period) && (T < 10.0 * Period)) {
            os << T << "," << y0(q_av) << "," << y0(p_vent) << "," << y0(v_v) << std::endl;
        }
    }


    Heart_AdValves(const std::string & id, Edge * e, Simple_vertex * sv,
              double density,
              double viscosity,

              double L_av,
              double L_pu,
              double L_mi,
              double B_pu_const,
              double B_av_denomin,
              double B_mi_denomin,

              double pulmVeinsPressure,
              double initLVvolume,
              double LV_V0, double LA_V0,
              double LV_ESPVR, double LV_EDPVR,
              double LA_ESPVR, double LA_EDPVR,
              double valvePressureForceCoeff,
              double valveFrictionalForceCoeff,
              double LV_inertiaCoeff,
              double LV_dynamicResistanceCoeff,
              double LA_inertiaCoeff,
              double LA_dynamicResistanceCoeff,
              double Period,
              double Ts1
              )
    : Base(id), Heart_part(
                    density,
                    viscosity,
                    L_av,
                    L_pu,
                    L_mi,
                    B_pu_const,
                    B_av_denomin,
                    B_mi_denomin,
                    pulmVeinsPressure,
                    LV_V0,  LA_V0,
                    LV_ESPVR, LV_EDPVR,
                    LA_ESPVR,  LA_EDPVR,
                    valvePressureForceCoeff,
                    valveFrictionalForceCoeff,
                    LV_inertiaCoeff,
                    LV_dynamicResistanceCoeff,
                    LA_inertiaCoeff,
                    LA_dynamicResistanceCoeff,
                    Period,
                    Ts1
                    ),
      e(e), sv(sv)
    {
        //see config try initLVvolume = 130 for healthy heart
        y0 << initLVvolume, 0, 20, 0, 0.01, 0, 0.01, 0, 1, 1, 1, 7.9, 20, 80, P_pulmonary;
        y2prev = y0;
        B.setZero();
        R.setZero();
        yn.setZero();
        Heart_part::set_heart_part_parameters(Base::parameterBlocks);

        for (BaseTimeCalculator & calc: this -> timeCalculators)
            calc.set_period(Period);
    }
protected:
    void update()
    {
        OutgoingCompatibilityCoeffs(e, sv, alfa, beta);

        double ts = T-dt;
        double tau = dt;
        while (tau > 1e-15) {

            double t = ts + tau;
            yn = y0 + (y0 - y2prev);

            int j = 0;
            double R0_norm = 0;

            int MaxIt = 20;
            while (1) {
                //Newton
                if (j > MaxIt) {
                    j = 0;
                    //yn = y0;
                    tau/=2;
                    t = ts + tau;

                    //std::cout << "!" << std::endl;
                    if (tau < 1e-15) {
                        std::cout << "LOL: Newton from Heart with Advanced valves failed" << std::endl;
                        std::cout << y0 << std::endl;
                        std::cout << yn << std::endl;
                        std::cout << "Res error: " << R.lpNorm<Eigen::Infinity>() << std::endl;

                        //break;
                        throw("Newton from Heart with Advanced valves failed");
                     }
                }

                R.head<11>() = yn.head<11>() - y0.head<11>() - tau * F(yn, t);
                R.tail<4>() = f(yn);
                if (R0_norm == 0) {
                    R0_norm = R.lpNorm<Eigen::Infinity>();
                }
                if (R.lpNorm<Eigen::Infinity>() < 1e-6 || R.lpNorm<Eigen::Infinity>() < 1e-14*R0_norm) {
                   // std::cout << R.lpNorm<Eigen::Infinity>() << std::endl;

                    break;
                }
                B.setZero();
                B.topLeftCorner<11, 11>().setIdentity();
                B.topRows<11>() -= tau * F_d(yn, t);
                B.bottomRows<4>() = f_d(yn);

                yn -= B.colPivHouseholderQr().solve(R);
                j++;
            }
            y2prev = y0;
            y0 = yn;
    /*
            const auto f_yynt = f_y(yn,t);
            const auto eigenvalues = f_yynt.eigenvalues();
            for (int i = 0; i < 8; i++) {
                if (eigenvalues(i).real() > 1)
                    cout << t << " " << eigenvalues(i) << endl;
            }
    */
            ts += tau;
            tau = T - ts;
        }
/*
        if (T >= 6 && T <= 7) {
            ff++;
            if (ff % 6 == 0)
                std::cout << y0(v1) << " " << P1(y0(s), y0(tet51)) << std::endl;
        }
        */
        if (y0(v_v) < 0 || y0(v_a) < 0) {
            std::cout << "time: " << T << std::endl;
            std::cout << "LV vol/pressure: " << y0(v_v) << " " << y0(p_vent) << std::endl;
            std::cout << "LA vol/pressure: " << y0(v_a) << " " << y0(p_atri) << std::endl;

            throw("OMG");
        }
        //save results
        LV_P = y0(p_vent);
        LV_V = y0(v_v);
        LA_P = y0(p_atri);
        LA_V = y0(v_a);
        aortic_valve = y0(ao_valve)/r2d;
        mitral_valve = y0(mi_valve)/r2d;
        const double S = y0(a_aorta);
        e->set_V_s(sv, 0, S);
        e->set_V_u(sv, 0, alfa*S + beta);
        if (y0(q_pu) < -1) {
            //std::cout << y0(q_pu) << std::endl;
        }
/*
        if (fmod(T, 1.0) > 0.999) {
            std::cout << T << " : " << this -> get_stroke_volume_av() << std::endl;
        }
*/
    }
};

#endif // HEART_ADVANCED_VALVES_H
