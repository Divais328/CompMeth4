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
                                              f(f),
                                              MatrixCopy(1, 1){};

Method::~Method(){}

void Method::InitialState()
{
    for (int n = 0; n <= GetN(); ++n)
    {
        for(int m  = 0; m <= GetM(); ++m)
        {
            SetValue(n, m, 0.);
        }
    }
}

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

void Method::GetMatrixCopy()
{
    MatrixCopy.Resize(GetN(), GetM());
    for (int i = 0; i <= GetN(); ++i)
    {
        for (int j = 0; j <= GetM(); ++j)
        {
            MatrixCopy.SetValue(i, j, GetValue(i, j));
        }
    }
}

double Method::GetMatrixCopyValue(int n, int m)
{
    return MatrixCopy.GetValue(n, m);
}

void Method::SaveSolution(std::string path)
{
    std::ofstream File;
    File.open(path, std::ofstream::out);
    if (!File.is_open())
    {
        std::cout << "Error at the file opening!" << std::endl;
    }
    else
    {
        for (int n = 0; n <= GetN();++n)
        {
            for (int m = 0; m <= GetM();++m)
            {
                File << GetX(n) << "\t" << GetY(m) << "\t" << GetValue(n, m)<<std::endl;
            }
            File << std::endl;
        }
        File.close();
    }
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
    double delta_max = 0;
    double u0 = 0;
    iter_num = 0;

    Reshape(N, M); // loose data
    InitialState();
    SetBorders();
    CalcF();

    do
    {
        delta = 0;
        delta_max = 0;
        GetMatrixCopy();
        for(int n  = 1; n < GetN(); ++n)
        {
            for(int m  = 1; m < GetM(); ++m)
            {
                u0 = GetValue(n, m);
                SetValue(n, m, 1. / 4. * (GetMatrixCopyValue(n, m - 1) + GetMatrixCopyValue(n - 1, m) + GetMatrixCopyValue(n + 1, m) + GetMatrixCopyValue(n, m + 1) - Getdx() * Getdy() * GetF(n, m)));
                delta = abs(u0 - GetValue(n, m));
                if (delta > delta_max)
                    delta_max = delta;
            }
        }
        ++iter_num;
    } while ((delta_max > e && iter_num < 100000) || iter_num < 5);
    return iter_num;
}


// Gauss Seidel Method Realisation
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
    double delta_max = 0;
    double u0 = 0;
    iter_num = 0;

    Reshape(N, M); // loose data
    InitialState();
    SetBorders();
    CalcF();

    do
    {
        delta = 0;
        delta_max = 0;
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
    } while ((delta_max > e && iter_num < 100000) || iter_num < 5);
    return iter_num;
}


// Successive OverRelaxation Method Realisation
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
    double delta = 0;
    double delta_max = 0;
    double u0 = 0;
    iter_num = 0;

    Reshape(N, M); // loose data
    InitialState();
    SetBorders();
    CalcF();

    do
    {
        GetMatrixCopy();
        delta = 0;
        delta_max = 0;
        for(int n  = 1; n < GetN(); ++n)
        {
            for(int m  = 1; m < GetM(); ++m)
            {
                SetValue(n, m, 1. / 4. * (GetValue(n, m - 1) + GetValue(n - 1, m) + GetValue(n + 1, m) + GetValue(n, m + 1) - Getdx() * Getdy() * GetF(n, m)));
            }
        }

        for(int n  = 1; n < GetN(); ++n)
        {
            for(int m  = 1; m < GetM(); ++m)
            {
                SetValue(n, m, (1. - w) * GetMatrixCopyValue(n, m) + w * GetValue(n,m));
                delta = abs(GetMatrixCopyValue(n, m) - GetValue(n, m));
                if (delta > delta_max)
                    delta_max = delta;
            }
        }
        
        ++iter_num;
    } while ((delta_max > e && iter_num < 100000) || iter_num < 5);
    return iter_num;
}


// Gauss Seidel 9 Points Method Realisation
Gauss_Seidel_9Points_Method::Gauss_Seidel_9Points_Method(double x0, double xN,
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

Gauss_Seidel_9Points_Method::~Gauss_Seidel_9Points_Method(){}

int Gauss_Seidel_9Points_Method::Solve(int N, int M, double e)
{
    double delta = 0;
    double delta_max = 0;
    double u0 = 0;
    iter_num = 0;

    Reshape(N, M); // loose data
    InitialState();
    SetBorders();
    CalcF();

    do
    {
        delta = 0;
        delta_max = 0;
        for(int n  = 1; n < GetN(); ++n)
        {
            for(int m  = 1; m < GetM(); ++m)
            {
                u0 = GetValue(n, m);
                SetValue(n, m, 1. / 20. * (GetValue(n - 1, m - 1) + GetValue(n + 1, m - 1) + GetValue(n - 1, m + 1) + GetValue(n + 1, m + 1) + 4. * (GetValue(n, m - 1) + GetValue(n - 1, m) + GetValue(n + 1, m) + GetValue(n, m + 1)) - 1. / 2. * Getdx() * Getdy() * (GetF(n - 1, m) + GetF(n + 1, m) + GetF(n, m - 1) + GetF(n, m + 1) + 8 * GetF(n, m))));
                delta = abs(u0 - GetValue(n, m));
                if (delta > delta_max)
                    delta_max = delta;
            }
        }
        ++iter_num;
    } while ((delta_max > e && iter_num < 100000) || iter_num < 5);
    return iter_num;
}