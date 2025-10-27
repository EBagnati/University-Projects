#include "Matrix.hpp"

void Matrix::resize(int rows, int cols)
{
    n_rows = rows;
    n_cols = cols;
    v.resize(rows * cols, 0.0);
}

double Matrix::operator()(int i, int j) const
{
    return v[i*n_cols + j];
}

double& Matrix::operator() (int i, int j)
{
    return v[i*n_cols + j];
}

double* Matrix::data()
{
    return v.data();
}

const double* Matrix::data() const
{
    return v.data();
}

double* Matrix::getRow(int i)
{
    return &v[i*n_cols];
}

void Matrix::getRowSubMatrix(int start_row, int end_row, Matrix& submatrix) const
{
    submatrix.resize(end_row - start_row, n_cols);
    for(int i = start_row; i < end_row; ++i)
    {
        for(int j = 0; j < n_cols; ++j)
        {
            submatrix(i - start_row, j) = v[i*n_cols + j];
        }
    }
}