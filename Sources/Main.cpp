#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <fstream>

#include "NonlinearSolvers.h"
#include "Functions.h"
#include "QRalgorithm.h"
#include "Matrix.h"
#include "PowerIteration.h"

constexpr double mainEPS = 1e-15;
constexpr double bisectionEPS = 1e-4;

template<typename U, typename V>
std::ostream &operator<<(std::ostream &out, const std::pair<U, V> &pair)
{
    out << "[" << pair.first << ", " << pair.second << "]" << std::endl;
    return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v)
{
    out << "[";
    for (auto i = 0; i < v.size() - 1; i++)
    {
        out << v[i] << ", ";
    }
    out << v.back() << "]";

    return out;
}

std::pair<EigenAnswer, double> task2(const Matrix &A)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

    start = std::chrono::high_resolution_clock::now();
    auto ans = PowerIterationMethod(A);
    end = std::chrono::high_resolution_clock::now();

    return {ans, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()};
}

double task3(Matrix A)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

    start = std::chrono::high_resolution_clock::now();
    auto ans = EigQR(std::move(A));
    end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

std::pair<std::pair<double, double>, int>
task6(const std::pair<double, double> &interval, const Func1D &f, double eps)
{
    int iterations = 0;
    auto ans = Bisection(interval, f, eps, iterations);
    return {ans, iterations};
}

std::pair<double, int>
task7(const std::pair<double, double> &interval, const Func1D &f, const Func1D &df)
{
    int iterations = 0;
    auto ans = NewtonsMethod((interval.first + interval.second) / 2, f, df, mainEPS, iterations);
    return {ans, iterations};
}

std::pair<double, int>
task8(const std::pair<double, double> &interval, const Func1D &f, double h)
{
    int iterations = 0;
    auto ans = DiscreteNewtonsMethod((interval.first + interval.second) / 2, f, h, mainEPS, iterations);
    return {ans, iterations};
}

void loop(int count)
{
    std::ofstream fout("report.txt");
    std::ofstream vout("norms.txt");

    double minTimeQR, maxTimeQR;
    minTimeQR = std::numeric_limits<double>::max();
    maxTimeQR = std::numeric_limits<double>::min();
    double allTimeQR = 0;

    double minTimePI, maxTimePI;
    minTimePI = std::numeric_limits<double>::max();
    maxTimePI = std::numeric_limits<double>::min();
    double allTimePI = 0;
    double avgNorm = 0;

    auto[firstRootBisection, firstRootIt] = task6({0.1, 2}, Functions::f, bisectionEPS);
    auto[secondRootBisection, secondRootIt] =task6({8, 10}, Functions::f, bisectionEPS);

    auto[firstNewton, firstNewtonIt] = task7(firstRootBisection, Functions::f, Functions::df);
    auto[secondNewton, secondNewtonIt] = task7(secondRootBisection, Functions::f, Functions::df);

    auto[firstDNewton, firstDNewtonIt] = task8(firstRootBisection, Functions::f, 0.001);
    auto[secondDNewton, secondDNewtonIt] = task8(secondRootBisection, Functions::f, 0.001);

    for (int i = 0; i < count; ++i)
    {
        Matrix A = Matrix::GenerateMatrix(10, variant);
        auto timeQR = task3(A);

        allTimeQR += timeQR;
        minTimeQR = std::min(minTimeQR, timeQR);
        maxTimeQR = std::max(maxTimeQR, timeQR);

        auto[ansPI, timePI] = task2(A);
        allTimePI += timePI;
        minTimePI = std::min(minTimePI, timePI);
        maxTimePI = std::max(maxTimePI, timePI);
        avgNorm += Utils::EigValueNorm(A, ansPI.vector, ansPI.value);

        vout << "Value: " << ansPI.value << std::endl;
        vout << "Vector: " << ansPI.vector << std::endl;
        vout << "Norm: " << Utils::EigValueNorm(A, ansPI.vector, ansPI.value) << std::endl;
        vout << std::endl;
    }

    fout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    fout << "Power method: " << std::endl;
    fout << "Min time: " << minTimePI << std::endl;
    fout << "Max time: " << maxTimePI << std::endl;
    fout << "Avg time: " << allTimePI / count << std::endl;
    fout << "Avg norm: " << avgNorm / count << std::endl;
    fout << std::endl;

    fout << "QR-algorithm: " << std::endl;
    fout << "Min time: " << minTimeQR << std::endl;
    fout << "Max time: " << maxTimeQR << std::endl;
    fout << "Avg time: " << allTimeQR / count << std::endl;
    fout << std::endl;

    fout << std::endl;
    fout << "First root: " << std::endl;
    fout << "Bisection: " << firstRootBisection << "n = " << firstRootIt << std::endl;
    fout << "Discrete Newton's: " << firstDNewton << std::endl << "n = " << firstDNewtonIt << std::endl;
    fout << "Newton's: " << firstNewton << std::endl << "n = " << firstNewtonIt << std::endl;

    fout << std::endl;
    fout << "Second root: " << std::endl;
    fout << "Bisection: " << secondRootBisection << "n = " << secondRootIt << std::endl;
    fout << "Discrete Newton's: " << secondDNewton << std::endl << "n = " << secondDNewtonIt << std::endl;
    fout << "Newton's: " << secondNewton << std::endl << "n = " << secondNewtonIt << std::endl;
    fout << std::endl;
}

void stepping(auto interval)
{
    for (double h = 1e-13; h < 1; h *= 1.1)
    {
        auto[val, it] = task8(interval, Functions::f, h);
        std::cout << h << " " << it << std::endl;
    }
}

void comparing(double x) {
    NewtonsMethod(x, Functions::f, Functions::df, 1e-32);
    std::cout << std::endl;
    DiscreteNewtonsMethod(x, Functions::f, 1e-3, 1e-32);
    std::cout << std::endl;
    DiscreteNewtonsMethod(x, Functions::f, 1e-5, 1e-32);
    std::cout << std::endl;
    DiscreteNewtonsMethod(x, Functions::f, 1e-7, 1e-32);
    std::cout << std::endl;
}

int main()
{
    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    double n = variant + 3;
    double A[] = {
            (-7+15*n)/2, (1+n)/2, (7-n)/2, (1+n)/2, (11-13*n)/2, (-7+9*n)/2, 1+n, 2*(-1+n), (3-5*n)/2, (-5+11*n)/2, (1-3*n)/2, (1+n)/2, (-1+3*n)/2, (1-3*n)/2, (-1+3*n)/2, (1-3*n)/2, 0, 0, (-1+3*n)/2, (1-3*n)/2, -29.0/2+2*n, -5.0/2+2*n, 5.0/2+n, -5.0/2+2*n, 29.0/2-2*n, -21.0/2+2*n, 0, -4, 13.0/2-2*n, -21.0/2+2*n, (17+3*n)/2, (5+n)/2, (-5+n)/2, 5*(1+n)/2, -3*(7+n)/2, 3*(5+n)/2, -2, 5, -3*(5+n)/2, 3*(7+n)/2, (-55+9*n)/2, (-11+n)/2, (11-n)/2, (-11+n)/2, (63-11*n)/2, 9*(-5+n)/2,
            3-n, -11+2*n, (33-5*n)/2, 7*(-7+n)/2, -59.0/2+5*n, -11.0/2+n, 11.0/2-n, -11.0/2+n, 63.0/2-7*n, -43.0/2+6*n, 2-2*n, 2*(-6+n), 35.0/2-3*n, -55.0/2+3*n, -2+5*n/2, n/2, -n/2, n/2, 3-7*n/2, -2+5*n/2, -n, -1+n, 1-3*n/2, -1+3*n/2, -5.0/2-4*n, -5.0/2, -3.0/2, -5.0/2, 11.0/2+3*n, -11.0/2-2*n, 1-n, 0, 5.0/2+n, -1.0/2-2*n, 5*(5+2*n), 7+2*n, -3-2*n, 7+2*n, -30-11*n, 23+8*n, -3-n, 3*(3+n), -4*(4+n), 21+6*n, 8+3*n, 2+n, -2-n, 2+n, -11-4*n, 8+3*n, -1-n, 3+n, -5-2*n, 2*(3+n)};
    auto m = Matrix::FromArray(A, 10, 10);
    loop(300);
    std::cout << EigQR(m);
    std::cout << std::endl;
}
