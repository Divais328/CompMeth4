#pragma once
#include <iostream>
#include "CalcGrid.hpp"

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

//Method difinition
class Method : public Grid2D
{
    Matrix2D F;
    double (&Border_x0)(double);
    double (&Border_xN)(double);
    double (&Border_y0)(double);
    double (&Border_yM)(double);
    double (&f)(double, double);

protected:
    void SetBorders();
    void CalcF();
    double GetF(int n, int m);

public:
    Method(double x0, double xN,
           double y0, double yM,
           double (&Border_x0)(double),
           double (&Border_xN)(double),
           double (&Border_y0)(double),
           double (&Border_yM)(double),
           double (&f)(double, double));
    ~Method();
};

// Jakobi Method Difinition
class Jakobi_Method : public Method //Grid2D
{
    // double delta = 0;
    int iter_num = 0;
    // double (&Border_x0)(double);
    // double (&Border_xN)(double);
    // double (&Border_y0)(double);
    // double (&Border_yM)(double);
    // double (&f)(double, double);

    // void SetBorders();

public:
    Jakobi_Method(double x0, double xN,
                  double y0, double yM,
                  double (&Border_x0)(double),
                  double (&Border_xN)(double),
                  double (&Border_y0)(double),
                  double (&Border_yM)(double),
                  double (&f)(double, double));
    ~Jakobi_Method();

    int Solve(int N, int M, double e);
};

// Gauss Seidel Method Difinition
class Gauss_Seidel_Method : public Method
{
    int iter_num = 0;

public:
    Gauss_Seidel_Method(double x0, double xN,
                        double y0, double yM,
                        double (&Border_x0)(double),
                        double (&Border_xN)(double),
                        double (&Border_y0)(double),
                        double (&Border_yM)(double),
                        double (&f)(double, double));
    ~Gauss_Seidel_Method();

    int Solve(int N, int M, double e);
};


// Successive OverRelaxation Method Difinition
class SOR_Method : public Method
{
    int iter_num = 0;
    double w = 1.5;
    
public:
    SOR_Method(double x0, double xN,
                                       double y0, double yM,
                                       double (&Border_x0)(double),
                                       double (&Border_xN)(double),
                                       double (&Border_y0)(double),
                                       double (&Border_yM)(double),
                                       double (&f)(double, double));
    ~SOR_Method();

    int Solve(int N, int M, double w, double e);
};