#include <iostream>
#include <iomanip>

#include "NonlinearSolvers.h"
#include "Functions.h"

template <typename U, typename V>
std::ostream& operator<<(std::ostream& out, const std::pair<U, V>& pair) {
    out << "[" << pair.first << ", " << pair.second << "]" << std::endl;
    return out;
}

int main()
{
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    auto[a, b] = Bisection({0, 1}, Functions::f, 1e-4);
    std::cout << a << ' ' << b << std::endl;
    double x = (a + b) / 2;
    std::cout << NewtonsMethod(x, Functions::f, Functions::df) << std::endl;
    std::cout << it << std::endl;
    std::cout << DiscreteNewtonsMethod(x, Functions::f, 1e-10) << std::endl;
    std::cout << it << std::endl;
}
