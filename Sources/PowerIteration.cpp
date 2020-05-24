#include "PowerIteration.h"

EigenAnswer RealEigen(const Vector &y, const Matrix &m)
{
    Vector u = Utils::Normalized(y);
    double lambda = Utils::ScalarMultiply(u, m * u);
    return {.vector = u, .value = lambda};
}

bool CheckEigen(const Matrix &m, std::initializer_list<EigenAnswer> c, double eps)
{
    double currentError = 0;

    for (const auto &now: c)
    {
        currentError += Utils::EuclideanNorm(m * now.vector - now.value * now.vector);
    }

    return eps > currentError;
}

EigenAnswer PowerIterationMethods(const Matrix &m, int &iter)
{
    const int n = m.GetRows();
    double currentError = std::numeric_limits<double>::max();
    Vector init(n);
    init[0] = 1;

    auto start = m * init;

    for (iter = 0; currentError >= defaultEPS; ++iter)
    {

    }
}
