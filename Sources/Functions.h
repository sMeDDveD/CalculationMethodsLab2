#pragma once

#include <functional>
#include <cmath>

using Func1D = std::function<double(double)>;

namespace Functions
{
    const Func1D f = [](double x)
    { return std::exp(std::sin(x)) - (x * x - 8.0f * x + 4.0f) / (2.0f * x); };

    const Func1D df = [](double x)
    { return std::cos(x) * std::exp(std::sin(x)) + 2.0f / (x * x) - 1.0f / 2; };
}
