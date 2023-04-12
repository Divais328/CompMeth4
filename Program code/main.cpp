#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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

double X0 = 0., XN = pi, Y0 = 0., YM = pi;

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
void SaveExactSolution()
{
    std::ofstream File;
    double dx = (XN - X0) / N, dy = (YM-Y0)/M;

    File.open("Exact_Solution_1024.txt", std::ofstream::out);
    if (!File.is_open())
    {
        std::cout << "Error at the file opening!" << std::endl;
    }
    else
    {
        for (int n = 0; n <= N;++n)
        {
            for (int m = 0; m <=M;++m)
            {
                File << X0 + n * dx << "\t" << Y0 + m * dy << "\t" << ExactSolution(X0 + n * dx, Y0 + m * dy) << std::endl;
            }
        }
        File.close();
    }
}
int N = 32, M = 32;
double e = 1.e-6, w = 1.7;

int main(int argc, char *argv[])
{
    int iter_num1 = 0;
    int iter_num2 = 0;
    int iter_num3 = 0;
    int iter_num4 = 0;
    std::string current_exec_name = argv[0]; // Name of the current exec program

    std::ofstream File;
    std::stringstream path;

    Jakobi_Method JM(X0, XN, Y0, YM, Border_x0, Border_xN, Border_y0, Border_yM, f);
    Gauss_Seidel_Method GSM(X0, XN, Y0, YM, Border_x0, Border_xN, Border_y0, Border_yM, f);
    SOR_Method SORM(X0, XN, Y0, YM, Border_x0, Border_xN, Border_y0, Border_yM, f);
    Gauss_Seidel_9Points_Method GS9PM(X0, XN, Y0, YM, Border_x0, Border_xN, Border_y0, Border_yM, f);

    if(argc == 1)
    {
        iter_num1 = JM.Solve(N, M, e);
        iter_num2 = GSM.Solve(N, M, e);
        iter_num3 = SORM.Solve(N, M, w, e);
        iter_num4 = GS9PM.Solve(N, M, e);

        std::cout << "Test mode: N = " << N << "M = " << M << "; e = " << e << "; w = " << w << std::endl
                  << "Exact value at pi/4 = " << ExactSolution(GSM.GetX(N / 4.), GSM.GetY(N / 4.)) << std::endl
                  << "Jacobi method: Number of iterrations  = " << iter_num1 << "; Value = " << JM.GetValue(N / 4., N / 4.) << std::endl
                  << "Gauss-Seidel method: Number of iterrations  = " << iter_num2 << "; Value = " << GSM.GetValue(N / 4., N / 4.) << std::endl
                  << "SOR method: Number of iterrations  = " << iter_num3 << "; Value = " << SORM.GetValue(N / 4., N / 4.) << std::endl
                  << "Gauss-Seidel 9 points method: Number of iterrations  = " << iter_num4 << "; Value = " << GS9PM.GetValue(N / 4., N / 4.) << std::endl;
    }
    else
    {
        if(atof(argv[1]) == 0)
        {
            if(argc > 2)
            {
                N = atof(argv[2]);
                M = N;
            }

            iter_num1 = JM.Solve(N, M, e);
            iter_num2 = GSM.Solve(N, M, e);
            iter_num3 = SORM.Solve(N, M, w, e);

            path << "JM_" << N << ".txt";
            JM.SaveSolution(path.str());
            path.clear();

            path << "GSM_" << N << ".txt";
            GSM.SaveSolution(path.str());
            path.clear();
            
            path << "SORM_" << N << "_" << w << ".txt";
            SORM.SaveSolution(path.str());
            path.clear();

            std::cout << "Save solution mode: N = " << N << "M = " << M << "; e = " << e << "; w = " << w << std::endl
                      << "Exact value at pi/4 = " << ExactSolution(GSM.GetX(N / 4.), GSM.GetY(N / 4.)) << std::endl
                      << "Jacobi method: Number of iterrations  = " << iter_num1 << "; Value = " << JM.GetValue(N / 4., N / 4.) << std::endl
                      << "Gauss-Seidel method: Number of iterrations  = " << iter_num2 << "; Value = " << GSM.GetValue(N / 4., N / 4.) << std::endl
                      << "SOR method: Number of iterrations  = " << iter_num3 << "; Value = " << SORM.GetValue(N / 4., N / 4.) << std::endl;
        }

        if(atof(argv[1]) == 1)
        {
            double w_opt = 0;
            double iter_num_min = 10000;
            path << "SORM_W_" << N << ".txt";
            File.open(path.str(), std::ofstream::out);
            if (!File.is_open())
            {
                std::cout << "Error at the file opening!" << std::endl;
            }
            else
            {
                w = 1.;
                while ((w<2&&iter_num3<10000)||iter_num_min>=10000)
                {
                    iter_num3 = SORM.Solve(N, M, w, e);
                    File << w << "\t" << iter_num3 << std::endl;
                    if(iter_num3<iter_num_min)
                    {
                        iter_num_min = iter_num3;
                        w_opt = w;
                    }
                    w = w + 0.01;
                }
                std::cout << "W optimization mode: N = " << N << "M = " << M << "; e = " << e << "; w = (1,2)" << std::endl
                          << "Optimal w = " << w_opt << "+-0.01" << std::endl;
                File.close();
            }
        }
    }
    return 0;
}