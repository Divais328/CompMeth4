#include "Poisson_Equation_Solvers.hpp"


// Method Realisation
Method::Method(double x0, double xN,
               double y0, double yM,
               double (&Border_x0)(double),
               double (&Border_xN)(double),
               double (&Border_y0)(double),
               double (&Border_yM)(double),
               double (&f)(double, double)) : Grid2D(1, x0, xN, 1, y0, yM),
                                              F(1, 1),
                                              Border_x0(Border_x0),
                                              Border_xN(Border_xN),
                                              Border_y0(Border_y0),
                                              Border_yM(Border_yM),
                                              f(f){}

Method::~Method(){}

void Method::SetBorders()
{
    for(int n  = 0; n < GetN(); ++n)
    {
        SetValue(n, 0, Border_y0(GetX(n)));
        SetValue(n, GetM(), Border_yM(GetX(n)));
    }
    for(int m  = 0; m < GetM(); ++m)
    {
        SetValue(0, m, Border_x0(GetY(m)));
        SetValue(GetN(), m, Border_xN(GetY(m)));
    }
}

void Method::CalcF()
{
    F.Resize(GetN(), GetM());
    for (int n = 0; n <= GetN(); ++n)
    {
        for(int m  = 0; m <= GetM(); ++m)
        {
            F.SetValue(n, m, f(GetX(n), GetY(m)));
        }
    }
}

double Method::GetF(int n, int m)
{
    return F.GetValue(n, m);
}

// Jakobi Method Realisation
Jakobi_Method::Jakobi_Method(double x0, double xN,
                             double y0, double yM,
                             double (&Border_x0)(double),
                             double (&Border_xN)(double),
                             double (&Border_y0)(double),
                             double (&Border_yM)(double),
                             double (&f)(double, double)) : Method(x0, xN, y0, yM,
                                                                   Border_x0,
                                                                   Border_xN,
                                                                   Border_y0,
                                                                   Border_yM,
                                                                   f){}

Jakobi_Method::~Jakobi_Method(){}

int Jakobi_Method::Solve(int N, int M, double e)
{
    double delta = 0;
    double delta_max = 1;
    double u0 = 0;
    iter_num = 0;

    Resize(N, M); // loose data
    SetBorders();

    do
    {
        for(int n  = 1; n < GetN(); ++n)
        {
            for(int m  = 1; m < GetM(); ++m)
            {
                u0 = GetValue(n, m);
                SetValue(n, m, 1. / 4. * (GetValue(n, m - 1) + GetValue(n - 1, m) + GetValue(n + 1, m) + GetValue(n, m + 1) - Getdx() * Getdy() * GetF(n, m)));
                delta = abs(u0 - GetValue(n, m));
                if (delta > delta_max)
                    delta_max = delta;
            }
        }
        ++iter_num;
    } while (delta_max > e);
    return iter_num;
}

// Gauss Seidel Method Realisatiom
Gauss_Seidel_Method::Gauss_Seidel_Method(double x0, double xN,
                                         double y0, double yM,
                                         double (&Border_x0)(double),
                                         double (&Border_xN)(double),
                                         double (&Border_y0)(double),
                                         double (&Border_yM)(double),
                                         double (&f)(double, double)) : Method(x0, xN, y0, yM,
                                                                               Border_x0,
                                                                               Border_xN,
                                                                               Border_y0,
                                                                               Border_yM,
                                                                               f){}

Gauss_Seidel_Method::~Gauss_Seidel_Method(){}

int Gauss_Seidel_Method::Solve(int N, int M, double e)
{
    double delta = 0;
    double delta_max = 1;
    double u0 = 0;
    iter_num = 0;

    Resize(N, M); // loose data
    SetBorders();

    do
    {
        for(int n  = 1; n < GetN(); ++n)
        {
            for(int m  = 1; m < GetM(); ++m)
            {
                u0 = GetValue(n, m);
                SetValue(n, m, 1. / 4. * (GetValue(n, m - 1) + GetValue(n - 1, m) + GetValue(n + 1, m) + GetValue(n, m + 1) - Getdx() * Getdy() * GetF(n, m)));
                delta = abs(u0 - GetValue(n, m));
                if (delta > delta_max)
                    delta_max = delta;
            }
        }
        ++iter_num;
    } while (delta_max > e);
    return iter_num;
}


// Successive OverRelaxation Method Realisatiom
SOR_Method::SOR_Method(double x0, double xN,
                       double y0, double yM,
                       double (&Border_x0)(double),
                       double (&Border_xN)(double),
                       double (&Border_y0)(double),
                       double (&Border_yM)(double),
                       double (&f)(double, double)) : Method(x0, xN, y0, yM,
                                                             Border_x0,
                                                             Border_xN,
                                                             Border_y0,
                                                             Border_yM,
                                                             f){};

SOR_Method::~SOR_Method(){}


int SOR_Method::Solve(int N, int M, double w, double e)
{
    
}