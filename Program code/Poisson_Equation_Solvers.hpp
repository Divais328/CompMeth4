#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
    Matrix2D MatrixCopy;
    Matrix2D F;
    double (&Border_x0)(double);
    double (&Border_xN)(double);
    double (&Border_y0)(double);
    double (&Border_yM)(double);
    double (&f)(double, double);

protected:
    void InitialState();
    void SetBorders();
    void CalcF();
    double GetF(int n, int m);
    void GetMatrixCopy();

public:
    Method(double x0, double xN,
           double y0, double yM,
           double (&Border_x0)(double),
           double (&Border_xN)(double),
           double (&Border_y0)(double),
           double (&Border_yM)(double),
           double (&f)(double, double));
    ~Method();
    double GetMatrixCopyValue(int n, int m);

    virtual int Solve(int N, int M, double e) = 0;
    virtual int Solve(int N, int M, double w, double e) = 0;

    void SaveSolution(std::string path);
};

// Jakobi Method Difinition
class Jakobi_Method : public Method //Grid2D
{
    int iter_num = 0;

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
    int Solve(int N, int M, double w, double e) { return 0; }
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
    int Solve(int N, int M, double w, double e) { return 0; }
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
    int Solve(int N, int M, double e) { return 0; }
};


// Gauss Seidel 9 Points Method Difinition
class Gauss_Seidel_9Points_Method : public Method
{
    int iter_num = 0;

public:
    Gauss_Seidel_9Points_Method(double x0, double xN,
                        double y0, double yM,
                        double (&Border_x0)(double),
                        double (&Border_xN)(double),
                        double (&Border_y0)(double),
                        double (&Border_yM)(double),
                        double (&f)(double, double));
    ~Gauss_Seidel_9Points_Method();

    int Solve(int N, int M, double e);
    int Solve(int N, int M, double w, double e) { return 0; }
};
