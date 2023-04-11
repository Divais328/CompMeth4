#include "CalcGrid.hpp"

//2D Matrix realisation
Matrix2D::Matrix2D(Matrix2D const& A): N(A.N), M(A.M),Matrix(new double *[N+1])
{
    Matrix[0] = new double [(N+1) * (M+1)];
    for (int j = 0; j <= M; ++j)
        {
            Matrix[0][j] = A.Matrix[0][j];
        }
    for (int i = 1; i <= N; ++i)
    {
        Matrix[i] = Matrix[i-1] + M+1;
        for (int j = 0; j <= M; ++j)
        {
            Matrix[i][j] = A.Matrix[i][j];
        }
    }

}
Matrix2D::Matrix2D(int N, int M): N(N), M(M)
{
    Matrix = new double *[N+1];
    Matrix[0] = new double[(N+1) * (M+1)];
    for (int j = 0; j <= M; ++j)
    {
        Matrix[0][j] = 0;
        }
    for (int i = 1; i <= N; ++i)
    {
        Matrix[i] = Matrix[i - 1] + M+1;
        for (int j = 0; j <= M; ++j)
        {
            Matrix[i][j] = 0;
        }
    }
}
Matrix2D::~Matrix2D()
{
    delete[] Matrix[0];
    delete[] Matrix;
}

void Matrix2D::Resize(int N, int M)
{
    this->N = N;
    this->M = M;

    delete[] Matrix[0];
    delete[] Matrix;

    Matrix = new double *[N+1];
    Matrix[0] = new double[(N+1) * (M+1)];
    for (int j = 0; j <= M; ++j)
        {
            Matrix[0][j] = 0;
        }
    for (int i = 1; i <= N; ++i)
    {
        Matrix[i] = Matrix[i-1] + M+1;
        for (int j = 0; j <= M; ++j)
        {
            Matrix[i][j] = 0;
        }
    }

}

int Matrix2D::GetN() const
{
    return N;
}
int Matrix2D::GetM() const
{
    return M;
}

double Matrix2D::GetValue(int n, int m) const
{
    return Matrix[n][m];
}
void Matrix2D::SetValue(int n, int m, double const value)
{
    Matrix[n][m] = value;
}

Matrix2D & Matrix2D::operator=(Matrix2D const &B)
{
    if(this != &B)
    {
        Resize(B.GetN(), B.GetM());
        for (int i = 0; i <= N; ++i)
        {
            Matrix[i] = Matrix[i-1] + M;
            for (int j = 0; j <= M; ++j)
            {
                Matrix[i][j] = B.Matrix[i][j];
            }
        }
    }
    return *this;
}

//3D Matrix realisation
Matrix3D::Matrix3D(Matrix3D const &A): P(A.P), N(A.N), M(A.M)
{
    Matrix = new double **[P+1];
    Matrix[0] = new double *[(P+1) * (N+1)];
    Matrix[0][0] = new double[(P+1) * (N+1) * (M+1)];
    for (int j = 1; j <= N; ++j)
    {
        Matrix[0][j] = Matrix[0][j - 1] + M+1;
        for (int k = 0; k <= M; ++k)
        {
            Matrix[0][j][k] = A.Matrix[0][j][k];
        }
    }
    for (int i = 1; i <= P; ++i)
    {
        Matrix[i] = Matrix[i - 1] + N+1;
        for (int j = 1; j <= N; ++j)
        {
            Matrix[i][j] = Matrix[i][j - 1] + M+1;
            for (int k = 0; k <= M; ++k)
            {
                Matrix[i][j][k] = A.Matrix[i][j][k];
            }
        }
    }

}
Matrix3D::Matrix3D(int P, int N, int M): P(P), N(N), M(M)
{
    Matrix = new double **[P+1];
    Matrix[0] = new double *[(P+1) * (N+1)];
    Matrix[0][0] = new double[(P+1) * (N+1) * (M+1)];
    for (int j = 1; j <= N; ++j)
    {
        Matrix[0][j] = Matrix[0][j - 1] + M + 1;
    }
    for (int i = 1; i <= P; ++i)
    {
        Matrix[i] = Matrix[i - 1] + N + 1;
        Matrix[i][0] = Matrix[i - 1][N] + M + 1;
        for (int j = 1; j <= N; ++j)
        {
            Matrix[i][j] = Matrix[i][j - 1] + M + 1;
        }
    }
    for (int i = 0; i <= P; ++i)
    {
        for (int j = 0; j <= N; ++j)
        {
            for (int k = 0; k <= M; ++k)
            {
                Matrix[i][j][k] = 0;
            }
        }
    }
}
Matrix3D::~Matrix3D()
{
    delete[] Matrix[0][0];
    delete[] Matrix[0];
    delete[] Matrix;
}

int Matrix3D::GetP() const
{
    return P;
}
int Matrix3D::GetN() const
{
    return N;
}
int Matrix3D::GetM() const
{
    return M;
}

double Matrix3D::GetValue(int p, int n, int m) const
{
    return Matrix[p][n][m];
}
void Matrix3D::SetValue(int p, int n, int m, double const value)
{
    Matrix[p][n][m] = value;
}

void Matrix3D::SetLayer(int p, Matrix2D const &A)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            Matrix[p][i][j] = A.GetValue(i,j);
        }
    }
}
Matrix2D& Matrix3D::GetLayer(int p) const
{
    Matrix2D A(N, M);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
           A.SetValue(i,j, Matrix[p][i][j]);
        }
    }
    return A;
}

//2D Calculation Grid realisation
CalcGrid2D::CalcGrid2D(int P, int N, int M): Matrix(P,N,M), dt(1./P), dx(1./N), dy(1./M) {}
CalcGrid2D::~CalcGrid2D() {}

int CalcGrid2D::GetP() const
{
    return Matrix.GetP();
}
int CalcGrid2D::GetN() const
{
    return Matrix.GetN();
}
int CalcGrid2D::GetM() const
{
    return Matrix.GetM();
}

double CalcGrid2D::Getdx() const
{
    return dx;
}
double CalcGrid2D::Getdy() const
{
    return dy;
}
double CalcGrid2D::Getdt() const
{
    return dt;
}

double CalcGrid2D::GetX(int n) const
{
    if ((n <= GetN()) && (n >= 0))
        return n * dx;
    else
        return -1;
}
double CalcGrid2D::GetY(int m) const
{
    if ((m <= GetM()) && (m >= 0))
        return m * dy;
    else
        return -1;
}
double CalcGrid2D::GetT(int p) const
{
    if ((p <= GetP()) && (p >= 0))
        return p * dx;
    else
        return -1;
}

void CalcGrid2D::SetValue(int p, int n, int m, double &value)
{
    Matrix.SetValue(p, n, m, value);
}
double CalcGrid2D::GetValue(int p, int n, int m) const
{
    return Matrix.GetValue(p, n, m);
}




//2D Grid realisation
Grid2D::Grid2D(int N, double x0, double xN,
               int M, double y0, double yM) : Matrix2D(N, M),
                                              x0(x0), xN(xN), dx((xN - x0) / N),
                                              y0(y0), yM(yM), dy((yM - y0) / M){};
Grid2D::~Grid2D() {}

double Grid2D::Getdx() const
{
    return dx;
}
double Grid2D::Getdy() const
{
    return dy;
}


double Grid2D::GetX(int n) const



{
    if ((n <= GetN()) && (n >= 0))
        return n * dx;
    else
        return -1;
}
double Grid2D::GetY(int m) const
{
    if ((m <= GetM()) && (m >= 0))
        return m * dy;
    else
        return -1;
}

