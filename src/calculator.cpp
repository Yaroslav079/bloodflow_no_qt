#include "calculator.h"

void BaseTimeCalculator::update(double new_dt)
{
    dt = new_dt;
    action();
    time += dt;
    if (time > Period) {
        time -= Period;
        val = new_val;
        new_val = init_value;
    }
}

BaseTimeCalculator::~BaseTimeCalculator()
{}

double BaseTimeCalculator::value() const
{
    return val;
}

BaseTimeCalculator::BaseTimeCalculator(double initValue, double Period)
    : val(initValue), init_value(initValue), Period(Period), time(0), new_val(initValue)
{}

BaseTimeCalculator::BaseTimeCalculator()
{}

void BaseTimeCalculator::set_period(double new_period) {
    this -> Period = new_period;
}


void BinaryCalculator::action()
{
    new_val = binary_op(new_val, next());
}

BinaryCalculator::BinaryCalculator(std::function<double(double, double)> binary_op,
                                   std::function<double ()> next, double init, double Period)
    : Base(init, Period), binary_op(binary_op), next(next)
{}

BinaryCalculator::~BinaryCalculator()
{}

BinaryCalculator::BinaryCalculator()
{}


void IntegralCalculator::action()
{
    new_val += integr_func() * dt;
}

IntegralCalculator::IntegralCalculator(std::function<double ()> interg_func, double Period)
    : Base(0, Period), integr_func(interg_func)
{}

IntegralCalculator::IntegralCalculator()
{}

IntegralCalculator::~IntegralCalculator()
{}




void AreaCalculator::action()
{
    points.push_back({x(), y()});
}

void AreaCalculator::finalize()
{
    //first, find the last point which is not ahead of very first point
    if (points.size() < 3)
        return;
    std::pair<double, double> vectorFirst = {points[1].first - points[0].first,
                                            points[1].second - points[0].second};
    while (points.size() > 3) {
        std::pair<double, double> vectorLast = {points.back().first - points[0].first,
                                                points.back().second - points[0].second};
        const double dot_prod = vectorFirst.first * vectorLast.first + vectorFirst.second * vectorLast.second;
        if (dot_prod > 0) {
            // last point is ahead of first point
            points.pop_back();
        } else {
            break;
        }
    }
    //now find area
    double res = 0;
    for (unsigned i = 0; i < points.size() - 1; i++) {
        res += (points[i].first + points[i+1].first) * (points[i+1].second - points[i].second);
    }
    //now close the curve
    res += (points.back().first + points[0].first) * (points[0].second - points.back().second);
    new_val = res / 2;
    //clear
    points.clear();
}

void AreaCalculator::update(double new_dt)
{
    dt = new_dt;
    action();
    time += dt;
    if (time > Period) {
        finalize();
        time -= Period;
        val = new_val;
        new_val = init_value;
    }
}

AreaCalculator::AreaCalculator(std::function<double ()> x, std::function<double ()> y, double Period)
    : Base(0, Period), x(x), y(y)
{}

AreaCalculator::AreaCalculator()
{}

AreaCalculator::~AreaCalculator()
{}


