#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>

#include "graph.h"
#include "edge.h"
#include "calculator.h"

typedef std::pair<std::string, double *> Parameter;
typedef std::vector<Parameter> ParameterBlock;
#define PAIR(x) Parameter((#x), &(x))

/*
 * ex vertex.cpp
Meta_vertex::Meta_vertex(const std::string & id)
    : id(id)
{}
Meta_vertex::~Meta_vertex()
{}

std::string Meta_vertex::get_id()
{
    return id;
}
*/

class Meta_vertex
{
public:
    std::string id;
    Meta_vertex(const std::string & id) : id(id) {}
    //virtual void update_boundary() = 0;
    virtual ~Meta_vertex() {}
    std::string get_id() {
        return id;
    }

    std::vector<ParameterBlock> parameterBlocks;
};


template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Vertex:
    public Meta_vertex
{
protected:
    typedef Meta_vertex Base;
    std::vector<std::reference_wrapper<BaseTimeCalculator>> timeCalculators;

    enum Vector_index {i_Sq, i_U};//square of cross-section, velocity
    const static int N = 2;

    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixDD;
    typedef Matrix<Scalar, N, Dynamic> MatrixND;
    typedef Matrix<Scalar, N, N> MatrixNN;
    typedef Matrix<Scalar, N, 1> VectorN;
    typedef Matrix<Scalar, Dynamic, 1> VectorD;
    typedef Matrix<Scalar, 1, Dynamic> RowVectorD;

    double dt;
    double T;//overall time
public:
    Vertex(const std::string & idd)
    : Base(idd)
    {
        dt = 0;
        T = 0;
    }

    void set_dt(double d)
    {
        dt = d;
    }
    void set_T(double d)
    {
        T = d;
    }
protected:
    void OCNFO2(const Edge * e, const Simple_vertex * sv, VectorN & OC, VectorN & OF, VectorN & OFF) const
    {
        OC(0) = e->get_Vp_s(sv, 0);
        OC(1) = e->get_Vp_u(sv, 0);
        OF(0) = e->get_Vp_s(sv, 1);
        OF(1) = e->get_Vp_u(sv, 1);
        OFF(0) = e->get_Vp_s(sv, 2);
        OFF(1) = e->get_Vp_u(sv, 2);
    }
    void OCNFO2_2nd(const Edge * e, const Simple_vertex * sv, VectorN & OC, VectorN & NC, VectorN & NF, VectorN & NFF) const
    {
        OC(0) = e->get_Vp_s(sv, 0);
        OC(1) = e->get_Vp_u(sv, 0);
        NC(0) = e->get_V_s(sv, 0);
        NC(1) = e->get_V_u(sv, 0);
        NF(0) = e->get_V_s(sv, 1);
        NF(1) = e->get_V_u(sv, 1);
        NFF(0) = e->get_V_s(sv, 2);
        NFF(1) = e->get_V_u(sv, 2);
    }
    /*
    void OutgoingCompatibilityCoeffs_old(const Edge * e, const Simple_vertex * sv, double & a, double & b) const
    {
        // Коэффициенты условия совместности (U = a * S + b) в начале i-ой ветви исходящей из узла
	    double sigm;
        double P, Ps, sound;
	    VectorN U, V, S, FPR;
	    VectorN OC, OF, OFF, NC, NF;
	    SqMatrixNN SL, SP, OM, OMO;

	    //Call OCNFO(br, Z, OC, OF, NC, NF)
        OCNFO2(e, OC, OF, OFF);
	
    	V = OC;
	
        e->Equation_of_state(V, sound, P, Ps);
	    e->Eigenvalues(V, S, SL, SP, sound, P);
        e->OMEGAB(V, OM, OMO, sound, P);
        FPR = e->get_FPRCH(sv, V, 0);
        const double dx = e->get_dx();

        sigm = S(0) * dt / dx;

    	a = - OM(0,0) / OM(0,1);
    
        //2nd order explicit without rhs
        //b = (1.0+3.0*sigm/2.0)*(-a*OC(0) + OC(1)) -
        //    2.0*sigm*(-a*OF(0) + OF(1)) + 0.5*sigm*(-a*OFF(0) + OFF(1));

        //with rhs
        b = (1+3*sigm/2)*(-a*OC(0) + OC(1)) -
            2*sigm*(-a*OF(0) + OF(1)) + 0.5*sigm*(-a*OFF(0) + OFF(1)) - (FPR(0)*a - FPR(1))*dt;
    }
*/
    void OutgoingCompatibilityCoeffs(const Edge * e, const Simple_vertex * sv, double & a, double & b) const
    {
        // Outgoing compatibility coefficients a, b for (U = a * S + b) at the beginning of Edge e
        double sigm;
        double P, Ps, sound;
        VectorN U, V, S, FPR;
        VectorN OC, NFF, NC, NF;
        MatrixNN SL, SP, OM, OMO;

        OCNFO2_2nd(e, sv, OC, NC, NF, NFF);
        //OC = NC;
        //V = NF;
        V = OC;

        e->get_pressure(V(0), sound, P, Ps, sv, 0);
        e->get_eigenvalues(V, S, SL, SP, sound, P);
        e->get_big_omega(V, OM, OMO, sound, P);
        FPR = e->get_right_hand_side(sv, V, 0); //not really implicit way
        const double dx = e->get_dx();
        sigm = S(0) * dt / dx;

        //implicit, actually not lol
        a = - OM(0,0) / OM(0,1);
        //b = (a * (sigm * NF(0) - OC(0)) + OC(1) - sigm * NF(1) + (FPR(1) - a * FPR(0)) * dt) / (1 - sigm);

        //b = (-a * NC(0)  + a * sigm * (NF(0) - NC(0)) + NC(1) - sigm * (NF(1) - NC(1)) + dt * (FPR(1) - a * FPR(0)));

         b = (a * (sigm * (2 * NF(0) - 0.5 * NFF(0)) - OC(0)) -
            (sigm * (2 * NF(1) - 0.5 * NFF(1)) - OC(1))
              + dt * (FPR(1) - a * FPR(0)))
                 / (1 - 3 * sigm / 2);

    }
    void setSU(Edge * e, const Simple_vertex * sv, const double & S, const double & U)
    {
        e->set_V_s(sv, 0, S);
        e->set_V_u(sv, 0, U);
    }
    double get_eos_deriv(const Edge * e, const double & S)
    {
        return e->Equation_of_state_deriv(S);
    }
    double get_pressure(const Edge * e, const double & S) const
    {
        return e->Equation_of_state_optimized(S);
    }
    virtual void update() = 0;
public:
    void update_boundary()
    {
        update();
        for (BaseTimeCalculator & calc: timeCalculators)
            calc.update(dt);
    }//need to think about the architecture

};


template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Heart:
    public Vertex<Edge, Matrix, Scalar, Dynamic>
{
    typedef Vertex<Edge, Matrix, Scalar, Dynamic> Base;

    IntegralCalculator strokeVOL;
    IntegralCalculator av_flow_fwd;
    IntegralCalculator av_flow_bwd;
    IntegralCalculator mv_flow_fwd;
    IntegralCalculator mv_flow_bwd;
    IntegralCalculator lvet; // Left Ventricular Ejection Time
    IntegralCalculator pft; // Peak Flow Time

protected:
    using Base::timeCalculators;
public:
    Heart(const std::string & id)
    : Base(id), strokeVOL(std::bind(&Heart::get_flow_av, this), 1),
      av_flow_fwd(std::bind(&Heart::get_flow_av_fwd, this), 1),
      av_flow_bwd(std::bind(&Heart::get_flow_av_bwd, this), 1),
      mv_flow_fwd(std::bind(&Heart::get_flow_mv_fwd, this), 1),
      mv_flow_bwd(std::bind(&Heart::get_flow_mv_bwd, this), 1),
      lvet(std::bind(&Heart::get_lvet, this), 1),
      pft(std::bind(&Heart::get_pft, this), 1)
    {
        timeCalculators.push_back(strokeVOL);
        timeCalculators.push_back(av_flow_fwd);
        timeCalculators.push_back(av_flow_bwd);
        timeCalculators.push_back(mv_flow_fwd);
        timeCalculators.push_back(mv_flow_bwd);
        timeCalculators.push_back(lvet);
        timeCalculators.push_back(pft);
    }
public:
    double get_stroke_volume_av()
    {
        return strokeVOL.value();
        //return stroke_volume_av;
    }
    virtual double get_stroke_volume_total()
    {
        return get_stroke_volume_av();
    }
    virtual double get_flow_av() = 0;

    double get_av_flow_bwd_val() {
        return av_flow_bwd.value();
    }
    virtual double get_flow_av_bwd() = 0;
    double get_av_flow_fwd_val() {
        return av_flow_fwd.value();
    }
    virtual double get_flow_av_fwd() = 0;

    double get_aortic_regurgitation_fraction() {
        return get_av_flow_bwd_val() / get_av_flow_fwd_val();
    }

    double get_mitral_regurgitation_fraction() {
        return get_mv_flow_bwd_val() / get_mv_flow_fwd_val();
    }

    double get_mv_flow_bwd_val() {
        return mv_flow_bwd.value();
    }
    virtual double get_flow_mv_bwd() = 0;
    double get_mv_flow_fwd_val() {
        return mv_flow_fwd.value();
    }
    virtual double get_flow_mv_fwd() = 0;

    virtual double get_lvet() = 0;
    double get_lvet_val() {
        return lvet.value();
    }

    virtual double get_pft() = 0;
    double get_pft_val() {
        return pft.value();
    }
};

template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class True_0d_heart:
    public Heart<Edge, Matrix, Scalar, Dynamic>
{
    typedef Heart<Edge, Matrix, Scalar, Dynamic> Base;

    IntegralCalculator aortic_valve_time_opened, mitral_valve_time_opened;
    BinaryCalculator max_aortic_pressure, min_aortic_pressure;
    AreaCalculator left_ventricular_work_pv;
    double integrate_aortic_valve_time_opened()
    {
        return (get_aortic_valve() > 65);
    }
    double integrate_mitral_valve_time_opened()
    {
        return (get_mitral_valve() > 65);
    }
protected:
    using Base::timeCalculators;
public:
    True_0d_heart(const std::string & id)
        : Base(id),
          aortic_valve_time_opened(std::bind(&True_0d_heart::integrate_aortic_valve_time_opened, this), 1),
          mitral_valve_time_opened(std::bind(&True_0d_heart::integrate_mitral_valve_time_opened, this), 1),
          max_aortic_pressure(static_cast<double const&(*)(double const&, double const&)>(std::max), std::bind(&True_0d_heart::get_aortic_root_pressure, this), 0, 1),
          min_aortic_pressure(static_cast<double const&(*)(double const&, double const&)>(std::min), std::bind(&True_0d_heart::get_aortic_root_pressure, this), 9999, 1),
          left_ventricular_work_pv(std::bind(&True_0d_heart::get_LV_V, this), std::bind(&True_0d_heart::get_LV_P, this), 1)
    {
        timeCalculators.push_back(aortic_valve_time_opened);
        timeCalculators.push_back(mitral_valve_time_opened);
        timeCalculators.push_back(max_aortic_pressure);
        timeCalculators.push_back(min_aortic_pressure);
        timeCalculators.push_back(left_ventricular_work_pv);
    }
    virtual double get_aortic_root_pressure() = 0;
    virtual double get_LV_P() = 0;
    virtual double get_LV_V() = 0;
    virtual double get_LA_P() = 0;
    virtual double get_LA_V() = 0;
    virtual double get_aortic_valve() = 0;
    double get_LV_work()
    {
        return left_ventricular_work_pv.value();
    }
    double get_time_aortic_valve_opened()
    {
        return aortic_valve_time_opened.value();
    }
    virtual double get_mitral_valve() = 0;
    double get_time_mitral_valve_opened()
    {
        return mitral_valve_time_opened.value();
    }
    double get_pulsePressure()
    {
        return max_aortic_pressure.value() - min_aortic_pressure.value();
    }
    double get_max_aortic_pressure()
    {
        return max_aortic_pressure.value();
    }
    double get_min_aortic_pressure()
    {
        return min_aortic_pressure.value();
    }
    virtual double get_ESPVR() = 0;
    virtual void set_ESPVR(double new_ESPVR) = 0;
    virtual double get_PveinPressure() = 0;
    virtual void set_PveinPressure(double new_PveinPressure) = 0;
    virtual void print_info() = 0;
    virtual void print_info(std::ostream &os) = 0;
    virtual double get_period() = 0;
};

template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Internal_vertex:
    public Vertex<Edge, Matrix, Scalar, Dynamic>
{
    typedef Vertex<Edge, Matrix, Scalar, Dynamic> Base;

    using Base::OutgoingCompatibilityCoeffs;
    using Base::OCNFO2;
    using typename Base::VectorN;
    using typename Base::VectorD;
    using typename Base::MatrixDD;

    int K;
    std::vector<Edge*> edges;
    Simple_vertex * sv;
    double density;
public:
    Internal_vertex(const std::string & id, const std::vector<Edge*> & edges, Simple_vertex * sv)
    : Base(id), edges(edges), sv(sv)
    {
        K = edges.size();
        density = edges[0]->get_density();
    }
private:
    /*
    void MakeRMatrix(MatrixDD & RMatrix, double & RDet) const
    {
        VectorD Resist(K);
        
        for (unsigned i = 0; i < K; i++) {
            Resist(i) = edges[i]->get_resistance();
        }
        // ! Подсчет определителя матрицы сопротивлений (было: "матрицы системы" ??)
        RDet = 0;

        for (unsigned i = 0; i < K; i++) {
            double W = 1;
            for (unsigned j = 0; j < K; j++) {
                if (j == i) continue;
                W *= Resist(j);
            }
            RDet += W;
        }

        // ! Подсчет диагональных элементов матрицы в системе для F
	    RMatrix.setZero();

        for (unsigned i = 0; i < K; i++) {
            double A = 0;
            for (unsigned j = 0; j < K; j++) {
                double W = 1;
                if (i == j) continue;
                for (unsigned m = 0; m < K; m++) {
                    if (m == i || m == j) continue;
                    W *= Resist(m);
                }
                A += W;
            }
            RMatrix(i, i) = -A;
        }

        // ! Подсчет недиагональных элементов матрицы в системе для F
        for (unsigned i = 0; i < K; i++) {
            for (unsigned j = 0; j < K; j++) {
                double W = 1;
                if (i == j) continue;
                for (unsigned m = 0; m < K; m++) {
                    if (m != i && m != j)
                        W *= Resist(m);
                }
                RMatrix(i, j) = W;
            }
        }
    }
*/
    void findSpeed(const VectorD & alf, const VectorD & bet, const VectorD & Xs, VectorD & Xu) const
    {
        for (int i = 0; i < K; i++) {
            Xu(i) = alf(i) * Xs(i) + bet(i);
        }
    }

    void find_F_Yac_Bern(VectorD & F, MatrixDD & YAC, const VectorD & Xs,
                         const VectorD & Xu, const VectorD & alf,
                         const VectorD & bet)
    {
        VectorD P(K), Ps(K);
        double sound;

        for (int i = 0; i < K; i++) {
          //  VectorN V;
           // V(0) = Xs(i);
            //V(1) = Xu(i);
            edges[i]->get_pressure(Xs(i), sound, P(i), Ps(i), sv, 0);
        }

        //find F
        F.setZero();
        F(0) =  Xu(0) * Xs(0);

        for (int i = 1; i < K; i++) {
            F(i) = P(i-1) - P(i) +
                    0.5 * density * (std::pow(Xu(i-1), 2) - std::pow(Xu(i), 2));
            F(0) += Xu(i) * Xs(i);
        }
        //F found


        //find YAC
        YAC.setZero();

        for (int i = 0; i < K; i++) {
            YAC(0,i) = 2 * alf(i) * Xs(i) + bet(i);
        }
        for (int i = 1; i < K; i++) {
            YAC(i,i-1) =  Ps(i-1) + alf(i-1) * density * Xu(i-1);
            YAC(i,i)   = -Ps(i) - alf(i) * density * Xu(i);
        }
        //done
    }

/*
    void FUZ(VectorD & F, VectorD & Xs, VectorD & Xu)
    {
        //	Подпрограмма вычисляет значение векторной функции определяющей 
        //	нелинейную систему в узле
	    VectorD P(K), Ps(K), alf(K), bet(K);
	    MatrixDD RM(K, K);
	    VectorN V;
	    double RDet, sound;

        for (unsigned i = 0; i < K; i++) {
            OutgoingCompatibilityCoeffs(edges[i], sv, alf(i), bet(i));
            V(0) = Xs(i);
		    V(1) = Xu(i);
            edges[i]->Equation_of_state(V, sound, P(i), Ps(i));
        }
        MakeRMatrix(RM, RDet);

        // ! Окончательное выражение для F

	    F.setZero();
	
        for (unsigned i = 0; i < K; i++) {
            for (unsigned j = 0; j < K; j++) {
			    if (i != j)
				    F(i) += RM(i,j) * P(j);
		    }
		    F(i) += RM(i,i) * P(i) - RDet * (alf(i) * Xs(i) + bet(i)) * Xs(i);
	    }
    }

    void TDYacobian(MatrixDD & YAC, const VectorD & Xs, const VectorD & Xu)
    {
        //	Подпрограмма вычисляет значение якобиана нелинейной системы в точке ветвления
	    VectorD P(K), Ps(K), alf(K), bet(K);
	    MatrixDD RM(K, K);
	    VectorN V;
	    double RDet, sound;
//this is the same code as in FUZ!!!!
        for (unsigned i = 0; i < K; i++) {
            OutgoingCompatibilityCoeffs(edges[i], sv, alf(i), bet(i));
            V(0) = Xs(i);
		    V(1) = Xu(i);
            edges[i]->Equation_of_state(V, sound, P(i), Ps(i));
        }
        MakeRMatrix(RM, RDet);
//end
        
        //! Подсчет элементов матрицы якобиана
	    YAC.setZero();
        for (unsigned i = 0; i < K; i++) {
	    	YAC(i,i) = -RDet * (2.0 * alf(i) * Xs(i) + bet(i)) + RM(i,i) * Ps(i);
            for (unsigned j = 0; j < K; j++) {
		    	if (i != j) YAC(i,j) = RM(i,j) * Ps(j);
            }
        }
    }
    */

public:
    void update() final
    {
        K = edges.size();
        MatrixDD YAC(K,K);
        VectorD F(K), Xs(K), Xu(K);
        //s - square
        //u - speed

        VectorD alf(K), bet(K); //compatibility coeffs

        for (int i = 0; i < K; i++) {
            OutgoingCompatibilityCoeffs(edges[i], sv, alf(i), bet(i));
            const Edge *e = edges[i];
            Xs(i) = e->get_Vp_s(sv, 0);
            Xu(i) = e->get_Vp_u(sv, 0);
        }

        //Solve system to find Xs with Netwons method
        const int Nmax = 20;//limit of iterations
        double F_norm = 1, F_normZero = 0;
        const double atol = 1e-6, rtol = 1e-13;
        for (int i = 0; F_norm > atol && F_norm > rtol*F_normZero; i++) {
	        if (i > Nmax) {
                for (int i = 0; i < K; i++) {
                    Edge *e = edges[i];
                    std::cout << "edge " << i << " " << e->get_id() << std::endl;
                }

                std::cout << "Xs: " << std::endl << Xs << std::endl;
                std::cout << "F: " << std::endl << F << std::endl;
                std::cout << sv->get_id() << std::endl;
                std::cout << "abs: " <<  F_norm << " rel: " << F_norm / F_normZero << std::endl;
                std::cout << "internalVertex Newton failed" << std::endl;
                throw("InternalVertex Newton failed");
            }

            find_F_Yac_Bern(F, YAC, Xs, Xu, alf, bet);

            //for large K dont use .inverse(), inverse implicitly!
            Xs -= YAC.inverse() * F;

            F_norm = F.template lpNorm<Eigen::Infinity>();
            if (i == 0)
                F_normZero = F_norm;

            findSpeed(alf, bet, Xs, Xu);
        }
        //done

        //send results to edges
        for (int i = 0; i < K; i++) {
#ifndef NDEBUG
         //   std::cout << Base::id << ": " << Xs(i) << " " << Xu(i) << " "<< Xs(i)*Xu(i) <<std::endl;
#endif
            Edge *e = edges[i];
            e->set_V_s(sv, 0, Xs(i));
            e->set_V_u(sv, 0, Xu(i));
        }
    }
};

template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Simple_heart:
    public Heart<Edge, Matrix, Scalar, Dynamic>
{
    typedef Heart<Edge, Matrix, Scalar, Dynamic> Base;
    using Base::OutgoingCompatibilityCoeffs;
    using Base::T;
    using Base::dt;
    using Base::setSU;

    Edge *e;
    Simple_vertex * sv;
    double Tc;
    double stroke_volume;
    double Q;
    void FlowToSU(double Qi, double & S, double & U)
    {
        double alfa, beta;
        OutgoingCompatibilityCoeffs(e, sv, alfa, beta);
        U = ( beta + std::sqrt(beta * beta + 4.0 * alfa * Qi) ) / 2;
        S = (-beta + std::sqrt(beta * beta + 4.0 * alfa * Qi) ) / (2 * alfa);
    }
public:
    Simple_heart(const std::string & id, Edge * e, Simple_vertex * sv, double period, double stroke_volume)
    : Base(id), e(e), sv(sv), Tc(period), stroke_volume(stroke_volume)
    {
        // std::cout <<  "Heart_period = " << Tc << std::endl;
    }
    double get_flow_av()
    {
        return Q;
    }

    // mock methods for correct reguratation fraction work with integrator
    double get_flow_av_fwd() {
        return 0.0;
    }
    double get_flow_av_bwd() {
         return 0.0;
    }
    double get_flow_mv_fwd() {
        return 0.0;
    }
    double get_flow_mv_bwd() {
         return 0.0;
    }

    // mock for LVET and pft
    double get_lvet() {
        return 0.0;
    }
    double get_pft() {
        return 0.0;
    }
private:
    void update()
    {
        double S, U;
        Q = stroke_volume / 58.75845 * (1.0/Tc) * 285 * (0.20617 + 0.37759*sin(2*M_PI*T/Tc + 0.59605) +
              0.2804*sin(4*M_PI*T/Tc - 0.35859) +  0.15337*sin(6*M_PI*T/Tc - 1.2509) -
              0.049889*sin(8*M_PI*T/Tc + 1.3921) +  0.038107*sin(10*M_PI*T/Tc - 1.1068) -
              0.041699*sin(12*M_PI*T/Tc + 1.3985) -  0.020754*sin(14*M_PI*T/Tc + 0.72921) +
              0.013367*sin(16*M_PI*T/Tc - 1.5394) -   0.021983*sin(18*M_PI*T/Tc + 0.95617) -
              0.013072*sin(20*M_PI*T/Tc - 0.022417) +   0.0037028*sin(22*M_PI*T/Tc - 1.4146) -
              0.013973*sin(24*M_PI*T/Tc + 0.77416) -   0.012423*sin(26*M_PI*T/Tc - 0.46511) +
              0.0040098*sin(28*M_PI*T/Tc + 0.95145) -   0.0059704*sin(30*M_PI*T/Tc + 0.86369) -
              0.0073439*sin(32*M_PI*T/Tc - 0.64769) +   0.0037006*sin(34*M_PI*T/Tc + 0.74663) -
              0.0032069*sin(36*M_PI*T/Tc + 0.85926) -   0.0048171*sin(38*M_PI*T/Tc - 1.0306) +
              0.0040403*sin(40*M_PI*T/Tc + 0.28009) -   0.0032409*sin(42*M_PI*T/Tc + 1.202) -
              0.0032517*sin(44*M_PI*T/Tc - 0.93316) +   0.0029112*sin(46*M_PI*T/Tc + 0.21405) -
              0.0022708*sin(48*M_PI*T/Tc + 1.1869) -  0.0021566*sin(50*M_PI*T/Tc - 1.1574) +
              0.0025511*sin(52*M_PI*T/Tc - 0.12915) -   0.0024448*sin(54*M_PI*T/Tc + 1.1185) -
              0.0019032*sin(56*M_PI*T/Tc - 0.99244) +   0.0019476*sin(58*M_PI*T/Tc - 0.059885) -
              0.0019477*sin(60*M_PI*T/Tc + 1.1655) -   0.0014545*sin(62*M_PI*T/Tc - 0.85829) +
              0.0013979*sin(64*M_PI*T/Tc + 0.042912) -   0.0014305*sin(66*M_PI*T/Tc + 1.2439) -
              0.0010775*sin(68*M_PI*T/Tc - 0.79464) +   0.0010368*sin(70*M_PI*T/Tc - 0.0043058) -
              0.0012162*sin(72*M_PI*T/Tc + 1.211) -   0.00095707*sin(74*M_PI*T/Tc - 0.66203) +
              0.00077733*sin(76*M_PI*T/Tc + 0.25642) -   0.00092407*sin(78*M_PI*T/Tc + 1.3954) -
              0.00079585*sin(80*M_PI*T/Tc - 0.49973));

        FlowToSU(Q, S, U);
        setSU(e, sv, S, U);
    }
};



template<typename Edge,
         template<typename, int, int ...> typename Matrix, typename Scalar, auto Dynamic>
class Constant_flow_heart:
    public Heart<Edge, Matrix, Scalar, Dynamic>
{
    typedef Heart<Edge, Matrix, Scalar, Dynamic> Base;
    using Base::OutgoingCompatibilityCoeffs;
    using Base::T;
    using Base::dt;
    using Base::setSU;

    Edge *e;
    Simple_vertex * sv;
    double flowrate;
    void FlowToSU(double Q, double & S, double & U)
    {
        double alfa, beta;
        OutgoingCompatibilityCoeffs(e, sv, alfa, beta);
        U = ( beta + std::sqrt(beta * beta + 4.0 * alfa * Q) ) / 2;
        S = (-beta + std::sqrt(beta * beta + 4.0 * alfa * Q) ) / (2 * alfa);
    }
public:
    Constant_flow_heart(const std::string & id, Edge * e, Simple_vertex * sv, double flowrate)
    : Base(id), e(e), sv(sv), flowrate(flowrate)
    {}

    double get_flow_av()
    {
        return flowrate;
    }
    // mock methods for correct reguratation fraction work with integrator
    double get_flow_av_fwd() {
        return 0.0;
    }
    double get_flow_av_bwd() {
         return 0.0;
    }
    double get_flow_mv_fwd() {
        return 0.0;
    }
    double get_flow_mv_bwd() {
         return 0.0;
    }

    // mock for LVET and pft
    double get_lvet() {
        return 0.0;
    }
    double get_pft() {
        return 0.0;
    }
protected:
    void update()
    {
        double S, U;
        FlowToSU(flowrate, S, U);
        setSU(e, sv, S, U);
    }
};

#endif
