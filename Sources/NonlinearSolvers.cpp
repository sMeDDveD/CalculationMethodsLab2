//
// Created by deem on 17.05.20.
//

#include <iostream>
#include "NonlinearSolvers.h"


double NewtonsMethod(double start, const Func1D &f, const Func1D &df, double eps, int &iterations)
{
    double curr = start;

    for (iterations = 0; iterations < defaultMaxIterations && std::abs(f(curr)) > eps; iterations++)
    {
        std::cout << std::abs(f(curr)) << " ";
        curr -= f(curr) / df(curr);
    }
    return curr;
}

double DiscreteNewtonsMethod(double start, const Func1D &f, double h, double eps, int &iterations)
{
    auto df = [h, f](double x)
    { return (f(x + h) - f(x)) / h; };

    return NewtonsMethod(start, f, df, eps, iterations);
}

std::pair<double, double> Bisection(std::pair<double, double> interval, const Func1D &f, double eps, int &iterations)
{
    auto[a, b] = interval;
    iterations = 0;
    while (std::abs(b - a) > eps)
    {
        iterations++;
        double mid = (b + a) / 2;
        if (f(mid) * f(a) < 0)
        {
            b = mid;
        }
        else
        {
            a = mid;
        }
    }

    return {a, b};
}


