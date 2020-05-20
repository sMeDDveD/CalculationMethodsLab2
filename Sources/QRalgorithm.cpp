//
// Created by deem on 19.05.20.
//

#include <iostream>
#include "QRalgorithm.h"

struct Rotation
{
    double c;
    double s;
};

static Vector GetW(Vector a)
{
    const int n = a.size();
    Vector r(n, 0);
    r[0] = -Utils::sgn(a[0]) * Utils::EuclideanNorm(a);

    Vector w = Utils::SubVectors(a, r);
    const double norm = Utils::EuclideanNorm(w);

    for (auto &now: w)
    {
        now /= norm;
    }

    return w;
}

static Rotation GetRotation(Matrix &m, int i)
{
    double a = m(i, i);
    double b = m(i + 1, i);

    return {.c = a / hypot(a, b), .s = b / hypot(a, b)};
}

void ToHessenbergForm(Matrix &m)
{
    int n = m.GetCols();
    std::cout << m << std::endl;

    for (int j = 0; j < n - 2; ++j)
    {
        Vector w = GetW(m.GetColPart(j, j + 1, n));
        int offset = n - w.size();

        // Transforming columns
        for (int column = j; column < n; ++column)
        {
            Vector colPart = m.GetColPart(column, j + 1, n);
            double dot = Utils::ScalarMultiply(w, colPart);
            for (int i = j + 1; i < n; ++i)
            {
                m(i, column) -= 2 * dot * w[i - offset];
            }
        }

        // Transforming rows
        for (int row = 0; row < n; ++row)
        {
            Vector rowPart = m.GetRowPart(row, j + 1, n);
            double dot = Utils::ScalarMultiply(w, rowPart);
            for (int i = j + 1; i < n; ++i)
            {
                m(row, i) -= 2 * dot * w[i - offset];
            }
        }
    }
}

void UndoRotation(Matrix &m, int j, const Rotation &r)
{
    const int n = m.GetCols();

    Vector lCol = m.GetCol(j);
    Vector rCol = m.GetCol(j + 1);

    for (int k = 0; k < n; k++)
    {
        m(k, j) = r.c * lCol[k] + r.s * rCol[k];
        m(k, j + 1) = -r.s * lCol[k] + r.c * rCol[k];
    }
}

void Rotate(Matrix &m, int i, const Rotation &r)
{
    const int n = m.GetCols();

    Vector uRow = m.GetRow(i);
    Vector lRow = m.GetRow(i + 1);

    for (int k = 0; k < n; k++)
    {
        m(i, k) = r.c * uRow[k] + r.s * lRow[k];
        m(i + 1, k) = -r.s * uRow[k] + r.c * lRow[k];
    }
}

void IterationQR(Matrix &h)
{
    const int n = h.GetCols();
    std::vector<Rotation> rotations(n - 1);

    for (int i = 0; i < n - 1; ++i)
    {
        auto r = GetRotation(h, i);
        rotations[i] = r;
        Rotate(h, i, r);
    }

    for (int j = 0; j < n - 1; ++j)
    {
        UndoRotation(h, j, rotations[j]);
    }
}
