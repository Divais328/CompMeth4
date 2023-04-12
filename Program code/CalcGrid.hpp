#pragma once
#include <iostream>

//2D Matrix difinition
class Matrix2D 
{
    double **Matrix = NULL;
    int N, M;

public:
    Matrix2D(Matrix2D const &A);
    Matrix2D(int N, int M);
    ~Matrix2D();

    void Resize(int N, int M); //loose  data

    int GetN() const;
    int GetM() const;
    double GetValue(int n, int m) const;
    void SetValue(int n, int m, double const value);

    Matrix2D & operator=(Matrix2D const &B);
};

//3D Matrix difinition
class Matrix3D 
{
    double ***Matrix = NULL;
    const int P, N, M;
public:
    Matrix3D(Matrix3D const &A);
    Matrix3D(int P, int N, int M);

    ~Matrix3D();

    int GetP() const;
    int GetN() const;
    int GetM() const;

    double GetValue(int p, int n, int m) const;
    void SetValue(int p, int n, int m, double const value);
};

//2D Calculation Grid difinition
class CalcGrid2D
{
    Matrix3D Matrix;
    const double dx, dy, dt;

public:
    CalcGrid2D(int P, int N, int M);
    ~CalcGrid2D();

    int GetP() const;
    int GetN() const;
    int GetM() const;

    double Getdx() const;
    double Getdy() const;
    double Getdt() const;

    double GetX(int n) const;
    double GetY(int m) const;
    double GetT(int p) const;

    void SetValue(int p, int n, int m, double &value);
    double GetValue(int p, int n, int m) const;
};

//2D Grid difinition
class Grid2D : public Matrix2D
{
    // Matrix3D Matrix;
    const double x0, xN, y0, yM;
    double dx, dy;

public:
    Grid2D(int N, double x0, double xN, int M, double y0, double yM);
    ~Grid2D();

    double Getdx() const;
    double Getdy() const;

    double GetX(int n) const;
    double GetY(int m) const;

    void Reshape(int N, int M);
};
