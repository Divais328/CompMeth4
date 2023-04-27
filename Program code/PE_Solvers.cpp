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
int N = 32, M = 32;
double e = 1.e-6, w = 1.69;

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
    N = 128;
    M = 128;
    double dx = (XN - X0) / N, dy = (YM - Y0) / M;

    File.open("Exact_Solution_128.txt", std::ofstream::out);
    if (!File.is_open())
    {
        std::cout << "Error at the file opening!" << std::endl;
    }
    else
    {
        for (int n = 0; n <= N; ++n)
        {
            for (int m = 0; m <=M; ++m)
            {
                File << X0 + n * dx << "\t" << Y0 + m * dy << "\t" << ExactSolution(X0 + n * dx, Y0 + m * dy) << std::endl;
            }
            File << std::endl;
        }
        File.close();
    }
}
double GetMaxError(Method &Meth)
{
    double error = 0, max_error = 0;
    for (int n = 0; n <= Meth.GetN(); ++n)
    {
        for (int m = 0; m <= Meth.GetM(); ++m)
        {
            error = abs(ExactSolution(Meth.GetX(n), Meth.GetY(m)) - Meth.GetValue(n, m));
            if (error > max_error)
            {
                max_error = error;
            }
        }
    }
    return max_error;
}
void SaveAllErrors(int N0, int k, Method &Meth, std::string path)
{
    std::ofstream File;
    
    File.open(path, std::ofstream::out);
    if (!File.is_open())
    {
        std::cout << "Error at the file opening!" << std::endl;
    }
    else
    {
        std::cout << typeid(Meth).name() << std::endl;
        for (int i = 0; i < k; ++i)
        {
            
        std::cout << "N  = " << N0 << "; nuber of iter = " << Meth.Solve(N0, N0, 1.e-9) << std::endl;

            File << N0 << "\t" << GetMaxError(Meth) << std::endl;
            N0 = N0 * 2;
        }
        File.close();
    }
}

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

        std::cout << std::endl
                  << "Test mode: N = " << N << ", M = " << M << "; e = " << e << "; w = " << w << std::endl
                  << "Exact value at pi/4 = " << ExactSolution(GSM.GetX(N / 4.), GSM.GetY(N / 4.)) << std::endl
                  << "Jacobi method: Number of iterrations  = " << iter_num1 << "; Value = " << JM.GetValue(N / 4., N / 4.) << std::endl
                  << "Gauss-Seidel method: Number of iterrations  = " << iter_num2 << "; Value = " << GSM.GetValue(N / 4., N / 4.) << std::endl
                  << "SOR method: Number of iterrations  = " << iter_num3 << "; Value = " << SORM.GetValue(N / 4., N / 4.) << std::endl
                  << "Gauss-Seidel 9 points method: Number of iterrations  = " << iter_num4 << "; Value = " << GS9PM.GetValue(N / 4., N / 4.) << std::endl
                  << std::endl;
    }
    else
    {
        //Save exact solution code
        if(atof(argv[1]) == 0)
        {
            N = 128;
            M = 128;
            std::cout << "Save exact solution mode: N = M = 128" << std::endl;
            SaveExactSolution();
        }
        //Save solution code
        if(atof(argv[1]) == 1)
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
            path.str("");

            path << "GSM_" << N << ".txt";
            GSM.SaveSolution(path.str());
            path.str("");
            
            path << "SORM_" << N << "_" << w << ".txt";
            SORM.SaveSolution(path.str());
            path.str("");

            std::cout << std::endl
                      << "Save solution mode: N = " << N << ", M = " << M << "; e = " << e << "; w = " << w << std::endl
                      << "Exact value at pi/4 = " << ExactSolution(GSM.GetX(N / 4.), GSM.GetY(N / 4.)) << std::endl
                      << "Jacobi method: Number of iterrations  = " << iter_num1 << "; Value = " << JM.GetValue(N / 4., N / 4.) << std::endl
                      << "Gauss-Seidel method: Number of iterrations  = " << iter_num2 << "; Value = " << GSM.GetValue(N / 4., N / 4.) << std::endl
                      << "SOR method: Number of iterrations  = " << iter_num3 << "; Value = " << SORM.GetValue(N / 4., N / 4.) << std::endl
                      << std::endl;
        }
        // W optimization code
        if(atof(argv[1]) == 2)
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
                std::cout << std::endl
                          << "W optimization mode: N = " << N << ", M = " << M << "; e = " << e << "; w = (1,2)" << std::endl
                          << "Optimal w = " << w_opt << "+-0.01" << std::endl
                          << std::endl;
                File.close();
            }
        }

        // Max error of N for Gauss-Seidel method (5 and 9 points)
        if(atof(argv[1]) == 3)
        {
            e = 1.e-12;

            std::cout << std::endl
                      << "Max error mode: N = M = (8,128); e = " << e << std::endl
                      << std::endl;

            SaveAllErrors(8, 5, GSM, "Errors_GSM.txt");
            SaveAllErrors(8, 5, GS9PM, "Errors_GS9PM.txt");
        }
    }
    return 0;
}