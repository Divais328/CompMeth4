#include <iostream>
#include <cmath>

#include "CalcGrid.hpp"
#include "Poisson_Equation_Solvers.hpp"

long double pi = acos(-1.);

// Borders
//                   yM
//     |---------------------------|
//     |                           |
//     |                           |
//     |                           |
//     |                           |
//  x0 |                           |  xN
//     |                           |
//     |                           |
//     |                           |
//     |                           |
//     |---------------------------|
//                   y0

double x0 = 0., xN = pi, y0 = 0., yM = pi;

double Border_x0 (double y)
{
    return 0;
}
double Border_xN (double y)
{
    return 0;
}
double Border_y0 (double x)
{
    return 0;
}
double Border_yM (double x)
{
    return 0;
}

// 8th variant
double f (double x, double y)
{
    return (2 / (x + 1) - (x - pi) / std::pow(x + 1, 2)) * sin(2 * y) - 4 * (x - pi) * log(x + 1) * sin(2 * y);
}
double ExactSolution (double x, double y)
{
    return (x - pi) * log(x + 1) * sin(2 * y);
}


int main(int argc, char *argv[])
{
    Gauss_Seidel_Method GSM(x0, xN, y0, yM, Border_x0, Border_xN, Border_y0, Border_yM, f);
    int iter_num1 = GSM.Solve(10, 10, 1.e-2);

    std::cout << iter_num1 << std::endl;

    return 0;
}