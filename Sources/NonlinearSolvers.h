#pragma once

#include "Functions.h"
#include "Defs.h"

double
NewtonsMethod(double start, const Func1D &f, const Func1D &df, double eps = defaultEPS, int &iteration = it);

double
DiscreteNewtonsMethod(double start, const Func1D &f, double h, double eps = defaultEPS, int &iteration = it);

std::pair<double, double>
Bisection(std::pair<double, double> interval, const Func1D &f, double eps = defaultEPS, int &iteration = it);
