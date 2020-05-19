﻿#pragma once

#include <algorithm>
#include <complex>

#include "Matrix.h"

using Complex = std::complex<double>;

Vector operator*(const Vector &l, double lambda);

Vector operator/(const Vector &l, double lambda);

Vector operator*(double lambda, const Vector &l);

Vector operator+(const Vector &l, const Vector &r);

Vector operator-(const Vector &l, const Vector &r);

namespace Utils
{
    template<typename T>
    int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    Complex solveQuadratic(double a, double b, double c)
    {
        return (-b + std::sqrt(Complex(b * b - 4 * a * c))) / (2 * a);
    }

    double CubicNorm(const Matrix &m);

    double EuclideanNorm(const Vector &v);

    double CubicNorm(const Vector &v);

    double ScalarMultiply(const Vector &l, const Vector &r);

    Vector GenerateVector(int length, int n);

    Vector SubVectors(const Vector &l, const Vector &r);

    Vector SolveUpperTriangle(const Matrix &m, const Vector &b);

    Vector SolveLowerTriangle(const Matrix &m, const Vector &b);

    std::pair<int, int> FindMax(const Matrix &m, int start);
}