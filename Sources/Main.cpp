#include <iostream>
#include <iomanip>

#include "NonlinearSolvers.h"
#include "Functions.h"
#include "QRalgorithm.h"
#include "Matrix.h"

template<typename U, typename V>
std::ostream &operator<<(std::ostream &out, const std::pair<U, V> &pair)
{
    out << "[" << pair.first << ", " << pair.second << "]" << std::endl;
    return out;
}

int main()
{
    double array[] = {
            -309, 263, -421, 33,
            -83, 87, -131, 38,
            226, -187, 316, -33,
            178, -150, 246, -23};
    auto m = Matrix::FromArray(array, 4, 4);
    ToHessenbergForm(m);
    std::cout << m;
    for (int i = 0; i < 100; i++)
    {
        IterationQR(m);
    }
    std::cout << m;
}
