#pragma once

#include "Matrix.h"
#include "Utils.h"
#include "Defs.h"

struct EigenAnswer
{
    Vector vector;
    double value;
};

bool CheckEigen(const Matrix &m, std::initializer_list<EigenAnswer> c, double eps);

EigenAnswer RealEigen(double y, const Vector &u, const Matrix &m);

EigenAnswer PowerIterationMethods(const Matrix &m, int &iter = it);
