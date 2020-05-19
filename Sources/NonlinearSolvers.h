#pragma once

#include "Functions.h"

constexpr double defaultEps = 1e-15;
constexpr int maxIterations = 1000;
static int it;

double
NewtonsMethod(double start, const Func1D &f, const Func1D &df, double eps = defaultEps, int &iteration = it);

double
DiscreteNewtonsMethod(double start, const Func1D &f, double h, double eps = defaultEps, int &iteration = it);

std::pair<double, double>
Bisection(std::pair<double, double> interval, const Func1D &f, double eps = defaultEps, int &iteration = it);
