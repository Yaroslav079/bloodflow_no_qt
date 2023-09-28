#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <functional>
#include <vector>

/**
 * @brief The BaseTimeCalculator class
 *
 * The idea of calculators is to evaluate some variables such as mean flow over period.
 * Calculator object should be bound to a specific edge object. Edge object provides calculator with some method
 * to obtain desired variable such as flow at time t in the edge center.
 * Calculator updates the value it tracks when update method is called.
 */
class BaseTimeCalculator
{
protected:
    double val;
    double init_value;
    double Period;
    double time;

    double new_val;
    double dt;
    virtual void action() = 0;
public:
    virtual void update(double dt);
    virtual ~BaseTimeCalculator();
    double value() const;
    BaseTimeCalculator(double initValue, double Period);
    BaseTimeCalculator();
    void set_period(double new_period);
};

class BinaryCalculator
    : public BaseTimeCalculator
{
protected:
    typedef BaseTimeCalculator Base;
    using Base::new_val;
    const std::function<double(double, double)> binary_op;
    const std::function<double()> next;
    virtual void action() final;
public:
    BinaryCalculator(std::function<double(double, double)> binary_op, std::function<double()> next, double init, double Period);
    virtual ~BinaryCalculator() override;
    BinaryCalculator();
};

class IntegralCalculator
    : public BaseTimeCalculator
{
protected:
    typedef BaseTimeCalculator Base;
    using Base::new_val;
    using Base::dt;
    const std::function<double()> integr_func;
    virtual void action() final;
public:
    IntegralCalculator(std::function<double()> interg_func, double Period);
    IntegralCalculator();
    virtual ~IntegralCalculator() override;
};

class AreaCalculator
    : public BaseTimeCalculator
{
protected:
    typedef BaseTimeCalculator Base;

    using Base::new_val;
    using Base::val;
    using Base::init_value;
    using Base::Period;
    using Base::time;
    using Base::dt;

    const std::function<double()> x, y;

    virtual void action() final;
    std::vector<std::pair<double, double>> points;
    void finalize();
public:
    virtual void update(double dt) override;
    AreaCalculator(std::function<double()> x, std::function<double()> y, double Period);
    AreaCalculator();
    virtual ~AreaCalculator() override;
};


#endif // CALCULATOR_H
