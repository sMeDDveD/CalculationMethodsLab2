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
                std::cout << m(i, column) << std::endl;
            }
            std::cout << m << std::endl;
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
