#ifndef HEART_PART_H
#define HEART_PART_H
#include <cmath>
#include "vertex.h"
class Heart_part
{
protected:
    const double cc = 1333.22;
    const double r2d = M_PI/180;

    double density;

    double viscosity;
    //coeffs for valve models
    double L_av, L_pu, L_mi;
    double B_pu_const;
    double B_av_denomin, B_mi_denomin;
protected:
    //pressure in pulmonary veins, mmHg
    double P_pulmonary;
private:
    //peak of systolic phase for ventricle
    // double Ts1;

    //end of systolic phase for ventricle
    double Ts2;

    double Tpb;

    double Tpw;
protected:
    //dead volume of lv
    double V1_0;

    //dead volume of la
    double V4_0;

    //elastance of ventricle in systolic phase, mmHg/ml
    //basically it is ESPVR
    double E1s;

    //elastance of ventricle in diastolic phase, mmHg/ml
    //this is EDPVR
    double E1d;

    //elastance of atrium in systolic phase, mmHg/ml
    double E4s;

    //elastance of atrium in diastolic phase, mmHg/ml
    double E4d;

    //coefficient of pressure force effect in valve
    double Kp;// = 10000;//was 5500

    //coefficient of frictional force effect in valve
    double Kf;// = 50;

    //max and min opening angles of both valves, degrees
    const double thetamax1_n = 75;
    const double thetamin1_n = 0;

    const double thetamax = thetamax1_n*r2d;
    const double thetamin = thetamin1_n*r2d;
    const double thetamin_mr = 60.0 * r2d;
    const int thetaPower = 4;
    const double AR0 = std::pow(1-cos(thetamax), thetaPower);

protected:
    //coefficient of inertia effect in ventricle
    double I1;

    //coefficient of resistance of ventricle
    double R1;

    //coefficient of inertia effect in atrium
    double I4;

    //cofficient of resistance of atrium
    double R4;

    //duration of heart period, s
    double Period;

    //peak of systole
    double Ts1;

protected:
    //time normalization
    inline double tnorm (const double & t) {
        return fmod(t, Period);
    }


    //ventricular activation function
    double e10 (const double & tcurr);
    //its time derivative
    double e10_t (const double & tcurr);

    //atrium activation function
    double e40 (const double & tcurr);
    double e40_t (const double & tcurr);

    //time-varying ventricular elastance
    double e1 (const double & t);

    double e1_t (const double & t);

    //time-varying atrium elastance
    double e4 (const double & t);

    double e4_t (const double & t);

    //valve opening function
    const double veps = 1e-6; //maybe doesnt work with 0? not sure

    double AR (const double & y);

    double AR_y (const double & y);

    double AR_mr(const double & y);

    double AR_mr_y (const double & y);


    //coefficient of resistance force in the valve
    const double ec = 1000;
    const double thetaeps = 0;

    double Kr (const double & y);

    double Kr_y (const double & y);

    double Kr_yy (const double & y);

    // for mitral valve regurgitation test
    double Kr_mr (const double & y);
    double Kr_mr_y (const double & y);
    double Kr_mr_yy (const double & y);



    //frictional force
    double Ff (const double & y3);
    double Ff_tet_d (const double & y3);

    //resistance force
    double Fr (const double & y2);
    double Fr_d (const double & y2);

    //resistance force for mitral regurgitation test
    double Fr_mr (const double & y2);
    double Fr_mr_d (const double & y2);


    //pressure force
    double Fp (const double & Pfrom, const double & Pto, const double & tet);
    double Fp_dPfrom(const double & Pfrom, const double & Pto, const double & tet);
    double Fp_dPto(const double & Pfrom, const double & Pto, const double & tet);
    double Fp_dTet (const double & Pfrom, const double & Pto, const double & tet);


    double L_av_inverse(const double & aorta_area);
    double L_av_inverse_d(const double & aorta_area);

    double L_pu_inverse();

    double L_mi_inverse();




    double valve_area(const double & valve_state);
    double valve_area_d(const double & valve_state);

    double valve_area_MI(const double & valve_state);
    double valve_area_MI_d(const double & valve_state);


    double R_av(const double & av_state);

    double R_av_d(const double & av_state);

    double R_pu();

    double B_pu();

    double B_av(const double & av_state, const double & aorta_area);
    double B_av_dAv_state(const double & av_state, const double & aorta_area);
    double B_av_dAorta_area(const double & av_state, const double & aorta_area);

    double B_mi(const double & mi_state);
    double B_mi_dMi_state(const double & mi_state);


    void set_heart_part_parameters(std::vector<ParameterBlock> & blocks);

public:
    Heart_part(
            double density,
            double viscosity,

            double L_av,
            double L_pu,
            double L_mi,
            double B_pu_const,
            double B_av_denomin,
            double B_mi_denomin,

            double pulmVeinsPressure,
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
        :  density(density), viscosity(viscosity),

          L_av(L_av),
          L_pu(L_pu),
          L_mi(L_mi),
          B_pu_const(B_pu_const),
          B_av_denomin(B_av_denomin),
          B_mi_denomin(B_mi_denomin),


          P_pulmonary(pulmVeinsPressure), V1_0(LV_V0), V4_0(LA_V0),
          E1s(LV_ESPVR),
          E1d(LV_EDPVR),
          E4s(LA_ESPVR),
          E4d(LA_EDPVR),
          Kp(valvePressureForceCoeff),
          Kf(valveFrictionalForceCoeff),
          I1(LV_inertiaCoeff),
          R1(LV_dynamicResistanceCoeff),
          I4(LA_inertiaCoeff),
          R4(LA_dynamicResistanceCoeff),
          Period(Period) //,
          // Ts1(Ts1)
    {
        this -> Ts1 = Ts1 * Period;
        Ts2 = this -> Ts1 + Period * 0.15;
        // Ts1 = 0.3 * Period;
        // Ts2 = 0.35 * Period;
        Tpb = 0.9 * Period;
        Tpw = 0.1 * Period;
    }
    double get_ESPVR();
    void set_ESPVR(double new_ESPVR);
    double get_PveinPressure();
    void set_PveinPressure(double new_PveinPressure);

    double get_Kp() {
        return Kp;
    }

    double get_Kf() {
        return Kf;
    }

    void set_Kp(const double &param_Kp) {
        Kp = param_Kp;
    }

    void set_Kf(const double &param_Kf) {
        Kf = param_Kf;
    }

    double get_I1() {
        return I1;
    }
    double get_R1() {
        return R1;
    }
    double get_I4() {
        return I4;
    }
    double get_R4() {
        return R4;
    }
    void set_chambers_coeffs(double& I1, double& R1, double& I4, double& R4) {
        this -> I1 = I1;
        this -> R1 = R1;
        this -> I4 = I4;
        this -> R4 = R4;
    }
};

#endif // HEART_PART_H
