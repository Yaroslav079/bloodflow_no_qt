#include "heart_part.h"
#include <assert.h>

double Heart_part::e10 (const double & tcurr) {

    const double t = tnorm(tcurr);
    assert(t >= 0);

    if (t < Ts1) {

        return 1.0 - cos(t*M_PI/Ts1);

    } else if (t < Ts2) {

        return  1.0 - cos((t-Ts2)*M_PI/(Ts1-Ts2));

    } else {

        return 0;

    }
}


double Heart_part::e10_t (const double & tcurr) {

    const double t = tnorm(tcurr);
    assert(t >= 0);

    if (t < Ts1) {

        return sin(t*M_PI/Ts1)*M_PI/Ts1;

    } else if (t <= Ts2) {

        return sin((t-Ts2)*M_PI/(Ts1-Ts2))*M_PI/(Ts1-Ts2);

    } else {

        return 0;

    }
}

double Heart_part::e40 (const double & tcurr) {

    const double t = tnorm(tcurr);
    assert(t >= 0);

    if (t < Tpb || t > Tpb + Tpw) {

        return 0;

    } else {

        return 1.0 - cos((t-Tpb)*2*M_PI/Tpw);

    }
}

double Heart_part::e40_t (const double & tcurr) {

    const double t = tnorm(tcurr);
    assert(t >= 0);

    if (t < Tpb || t > Tpb + Tpw) {

        return 0;

    } else {

        return sin((t-Tpb)*2*M_PI/Tpw)*2*M_PI/Tpw;

    }
}

double Heart_part::e1 (const double & t) {

    return E1d + (E1s - E1d) * e10(t) / 2.0;

}

double Heart_part::e1_t (const double & t) {

    return  (E1s - E1d) * e10_t(t) / 2.0;

}
double Heart_part::e4 (const double & t) {

    return E4d + (E4s - E4d) * e40(t) / 2.0;

}


double Heart_part::e4_t (const double & t) {

    return (E4s - E4d) * e40_t(t) / 2.0;

}


double Heart_part::AR (const double & y) {

    if (y <= thetamin) {

        return veps;

    } else if (y >= thetamax) {

        return 1;

    } else {

        return (1-veps)*std::pow(1-cos(y), thetaPower)/AR0 + veps;

    }
}

double Heart_part::AR_aortic_reg (const double & y) {

    if (y <= thetamin_aortic_reg) {

        return veps;

    } else if (y >= thetamax) {

        return 1;

    } else {

        return (1-veps)*std::pow(1-cos(y), thetaPower)/AR0 + veps;

    }
}

double Heart_part::AR_mitral_reg (const double & y) {

    if (y <= thetamin_mitral_reg) {

        return veps;

    } else if (y >= thetamax) {

        return 1;

    } else {

        return (1-veps)*std::pow(1-cos(y), thetaPower)/AR0 + veps;

    }
}

double Heart_part::AR_y (const double & y) {

    if (y < thetamax && y > thetamin) {

        return (1-veps)*thetaPower*std::pow(1-cos(y), thetaPower - 1)*sin(y)/AR0;

    } else {

        return 0;

    }
}

double Heart_part::AR_aortic_reg_y (const double & y) {

    if (y < thetamax && y > thetamin_aortic_reg) {

        return (1-veps)*thetaPower*std::pow(1-cos(y), thetaPower - 1)*sin(y)/AR0;

    } else {

        return 0;

    }
}

double Heart_part::AR_mitral_reg_y (const double & y) {

    if (y < thetamax && y > thetamin_mitral_reg) {

        return (1-veps)*thetaPower*std::pow(1-cos(y), thetaPower - 1)*sin(y)/AR0;

    } else {

        return 0;

    }
}

double Heart_part::Kr (const double & y) {
    if (y >= thetamax) {
        return expm1(ec*(y-thetamax));
    } else if (y <= thetamin+thetaeps) {
        return -expm1(ec*(thetamin+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_y (const double & y) {
    if (y >= thetamax) {
        return ec*exp(ec*(y-thetamax));
    } else if (y <= thetamin+thetaeps) {
        return ec*exp(ec*(thetamin+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_yy (const double & y) {
    if (y >= thetamax) {
        return ec*ec*exp(ec*(y-thetamax));
    } else if (y <= thetamin+thetaeps) {
        return -ec*ec*exp(ec*(thetamin+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_aortic_reg (const double & y) {
    if (y >= thetamax) {
        return expm1(ec*(y-thetamax));
    } else if (y <= thetamin_aortic_reg+thetaeps) {
        return -expm1(ec*(thetamin_aortic_reg+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_aortic_reg_y (const double & y) {
    if (y >= thetamax) {
        return ec*exp(ec*(y-thetamax));
    } else if (y <= thetamin_aortic_reg+thetaeps) {
        return ec*exp(ec*(thetamin_aortic_reg+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_aortic_reg_yy (const double & y) {
    if (y >= thetamax) {
        return ec*ec*exp(ec*(y-thetamax));
    } else if (y <= thetamin_aortic_reg+thetaeps) {
        return -ec*ec*exp(ec*(thetamin_aortic_reg+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_mitral_reg (const double & y) {
    if (y >= thetamax) {
        return expm1(ec*(y-thetamax));
    } else if (y <= thetamin_mitral_reg+thetaeps) {
        return -expm1(ec*(thetamin_mitral_reg+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_mitral_reg_y (const double & y) {
    if (y >= thetamax) {
        return ec*exp(ec*(y-thetamax));
    } else if (y <= thetamin_mitral_reg+thetaeps) {
        return ec*exp(ec*(thetamin_mitral_reg+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Kr_mitral_reg_yy (const double & y) {
    if (y >= thetamax) {
        return ec*ec*exp(ec*(y-thetamax));
    } else if (y <= thetamin_mitral_reg+thetaeps) {
        return -ec*ec*exp(ec*(thetamin_mitral_reg+thetaeps-y));
    } else {
        return 0;
    }
}

double Heart_part::Ff (const double & y3) {
    return Kf*y3;
}

double Heart_part::Ff_tet_d (const double & y3) {
    return Kf;
}

double Heart_part::Fr (const double & y2) {
    return Kr(y2);
}

double Heart_part::Fr_d (const double & y2) {
    return Kr_y(y2);
}

double Heart_part::Fr_aortic_reg (const double & y2) {
    return Kr_aortic_reg(y2);
}

double Heart_part::Fr_aortic_reg_d (const double & y2) {
    return Kr_aortic_reg_y(y2);
}

double Heart_part::Fr_mitral_reg (const double & y2) {
    return Kr_mitral_reg(y2);
}

double Heart_part::Fr_mitral_reg_d (const double & y2) {
    return Kr_mitral_reg_y(y2);
}

double Heart_part::Fp (const double & Pfrom, const double & Pto, const double & tet) {
    return Kp*(Pfrom - Pto) * cos(tet);
}
double Heart_part::Fp_dPfrom(const double & Pfrom, const double & Pto, const double & tet) {
    return Kp*cos(tet);
}
double Heart_part::Fp_dPto(const double & Pfrom, const double & Pto, const double & tet) {
    return  - Kp * cos(tet);
}
double Heart_part::Fp_dTet (const double & Pfrom, const double & Pto, const double & tet) {
    return -Kp*(Pfrom - Pto) * sin(tet);
}
double Heart_part::L_av_inverse_d(const double & aorta_area)
{
    //const double Ku = 1, L_sten = 1;
//      return 1.0 / (Ku*density*L_sten) * cc;
    return 0;
}
double Heart_part::L_av_inverse(const double & aorta_area) { return 1.0 / L_av; }
double Heart_part::L_pu_inverse() { return 1.0 / L_pu; }
double Heart_part::L_mi_inverse() { return 1.0 / L_mi; }

double Heart_part::R_av(const double & av_state)
{
    return 0;
}

double Heart_part::R_av_d(const double & av_state)
{
    return 0;//TODO
}

double Heart_part::R_pu()
{
    return 0;
}

double Heart_part::valve_area(const double & valve_state)
{
    return 4 * AR(valve_state);
}
double Heart_part::valve_area_d(const double & valve_state)
{
    return 4 * AR_y(valve_state);
}

double Heart_part::valve_area_MI(const double & valve_state)
{
    return 4 * AR(valve_state);//5
}
double Heart_part::valve_area_MI_d(const double & valve_state)
{
    return 4 * AR_y(valve_state);//5
}

double Heart_part::valve_area_aortic_reg(const double & valve_state) {
    return 4 * AR_aortic_reg(valve_state);
}
double Heart_part::valve_area_aortic_reg_d(const double & valve_state) {
    return 4 * AR_aortic_reg_y(valve_state);
}

double Heart_part::valve_area_mitral_reg(const double & valve_state) {
    return 4 * AR_mitral_reg(valve_state);
}
double Heart_part::valve_area_mitral_reg_d(const double & valve_state) {
    return 4 * AR_mitral_reg_y(valve_state);
}

double Heart_part::B_pu()
{
   // const double Kt = 1;
   // const double pu_area = 3;
    //return density * Kt/2 * std::pow(1.0/pu_area - 1.0/5, 2)/4;
    return B_pu_const;
}

double Heart_part::B_av(const double & av_state, const double & aorta_area)
{
    const double Kt = 1;
    const double av_area = valve_area(av_state);
    return density * Kt/2 *
            std::pow(1.0/av_area - 1.0/aorta_area, 2) / B_av_denomin;
}
double Heart_part::B_av_dAv_state(const double & av_state, const double & aorta_area)
{
    const double Kt = 1;
    const double av_area = valve_area(av_state);
    return - density * Kt *
            (1.0/av_area - 1.0/aorta_area) / std::pow(av_area, 2)
            * valve_area_d(av_state) / B_av_denomin;
}
double Heart_part::B_av_dAorta_area(const double & av_state, const double & aorta_area)
{
    const double Kt = 1;
    const double av_area = valve_area(av_state);
    return density * Kt * (1.0/av_area - 1.0/aorta_area)
            / std::pow(aorta_area, 2) / B_av_denomin;
}

double Heart_part::B_av_aortic_reg(const double & av_state, const double & aorta_area) {
    const double Kt = 1;
    const double av_area = valve_area_aortic_reg(av_state);
    return density * Kt/2 *
            std::pow(1.0/av_area - 1.0/aorta_area, 2) / B_av_denomin;
}
double Heart_part::B_av_dAv_state_aortic_reg(const double & av_state, const double & aorta_area) {
    const double Kt = 1;
    const double av_area = valve_area_aortic_reg(av_state);
    return - density * Kt *
            (1.0/av_area - 1.0/aorta_area) / std::pow(av_area, 2)
            * valve_area_aortic_reg_d(av_state) / B_av_denomin;
}
double Heart_part::B_av_dAorta_area_aortic_reg(const double & av_state, const double & aorta_area) {
    const double Kt = 1;
    const double av_area = valve_area_aortic_reg(av_state);
    return density * Kt * (1.0/av_area - 1.0/aorta_area)
            / std::pow(aorta_area, 2) / B_av_denomin;
}

double Heart_part::B_mi(const double & mi_state)
{
    const double Kt = 1;
    const double mi_area = valve_area_MI(mi_state);
    return density * Kt/2 * std::pow(1.0/mi_area - 1.0/5, 2) / B_mi_denomin;
}
double Heart_part::B_mi_dMi_state(const double & mi_state)
{
    const double Kt = 1;
    const double mi_area = valve_area_MI(mi_state);
    return - density * Kt * (1.0/mi_area - 1.0/5) / std::pow(mi_area, 2)
            * valve_area_MI_d(mi_state) / B_mi_denomin;
}

double Heart_part::B_mi_mitral_reg(const double & mi_state) {
    const double Kt = 1;
    const double mi_area = valve_area_mitral_reg(mi_state);
    return density * Kt/2 * std::pow(1.0/mi_area - 1.0/5, 2) / B_mi_denomin;
}
double Heart_part::B_mi_dMi_state_mitral_reg(const double & mi_state) {
    const double Kt = 1;
    const double mi_area = valve_area_mitral_reg(mi_state);
    return - density * Kt * (1.0/mi_area - 1.0/5) / std::pow(mi_area, 2)
            * valve_area_mitral_reg_d(mi_state) / B_mi_denomin;
}



double Heart_part::get_ESPVR()
{
    return E1s;
}
void Heart_part::set_ESPVR(double new_ESPVR)
{
    E1s = new_ESPVR;
}
double Heart_part::get_PveinPressure()
{
    return P_pulmonary;
}
void Heart_part::set_PveinPressure(double new_PveinPressure)
{
    P_pulmonary = new_PveinPressure;
}
void Heart_part::set_heart_part_parameters(std::vector<ParameterBlock> & blocks)
{
    blocks.push_back({PAIR(E1d), PAIR(E1s), PAIR(V1_0), PAIR(I1), PAIR(R1)});
    blocks.push_back({PAIR(E4d), PAIR(E4s), PAIR(V4_0), PAIR(I4), PAIR(R4)});
    blocks.push_back({PAIR(Kf), PAIR(Kp)});
    blocks.push_back({PAIR(L_av), PAIR(B_av_denomin)});
    blocks.push_back({PAIR(L_pu), PAIR(B_pu_const), PAIR(P_pulmonary)});
    blocks.push_back({PAIR(L_mi), PAIR(B_mi_denomin)});
    blocks.push_back({PAIR(Period)});
}
