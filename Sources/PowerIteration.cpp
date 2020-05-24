#include "PowerIteration.h"

constexpr double PIeps = 1e-13;

EigenAnswer RealValues(const Vector &x1, const Vector &x2)
{
    const Vector &ans = x2;
    double value = 0;

    // average?
    for (int i = 0; i < ans.size(); ++i)
    {
        if (x1[i] != 0)
        {
            value = std::max(value, x2[i] / x1[i]);
        }
    }

    return {.vector = ans, .value = value};
}

EigenAnswer OppositeRealValues(const Vector &x1, const Vector &x2, const Vector &x3)
{
    double value = 0;
    // average?
    for (int i = 0; i < x3.size(); ++i)
    {
        if (x1[i] != 0)
        {
            value = std::max(value, x3[i] / x1[i]);
        }
    }

    value = std::sqrt(value);

    // second?
    Vector ans = x3;

    for (int i = 0; i < x3.size(); ++i)
    {
        ans[i] += value * x2[i];
    }

    return {.vector = ans, .value = value};
}

bool TestEigen(const Matrix &m, std::initializer_list<EigenAnswer> c, double eps)
{
    double currentError = 0;

    for (const auto &now: c)
    {
        currentError += Utils::EigValueNorm(m, now.vector, now.value);
    }

    return eps > currentError;
}

EigenAnswer PowerIterationMethod(const Matrix &m, int &iter)
{
    const int n = m.GetRows();
    Vector initial(n);
    initial[0] = 1;

    EigenAnswer bad;
    for (iter = 0; iter < defaultMaxIterations; iter += 3)
    {
        auto first = initial;
        auto second = m * first;
        auto third = m * second;

        auto real = RealValues(second, third);
        if (TestEigen(m, {real}, PIeps))
        {
            return real;
        }

        auto opposite = OppositeRealValues(first, second, third);
        if (TestEigen(m, {opposite}, PIeps))
        {
            return opposite;
        }
        bad = opposite;
        initial = Utils::Normalized(m * third);
        iter++;
    }

    return bad;
}
