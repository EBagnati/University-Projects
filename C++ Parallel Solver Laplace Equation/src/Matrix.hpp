#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>


/*!
 * Class that implements a matrix with all elements stored contiguously in a single vector
 */
class Matrix
{
public:

    // Default constructor
    Matrix() = default;

    /*!
     * Constructor that takes number of rows and columns
     * @param rows Number of rows
     * @param cols Number of columns
     */
    Matrix(int rows, int cols) : n_rows(rows), n_cols(cols)
    {
        v.resize(rows * cols, 0.0);
    }

    /*! Const version of the call operator
     * @param i Row index
     * @param j Column index
     * @return The value
     */
    double operator()(int i, int j) const;

    /*! Non const version of the call operator
     * @param i Row index
     * @param j Column index
     * @return The value (which can be changed)
     */
    double& operator() (int i, int j);

    /*! Resizes the matrix
     * @param rows The new number of rows
     * @param cols The new number of columns
     */
    void resize(int rows, int cols);


    /*! Data method to access to the pointer to the internal vector
     * @return Pointer to the internal vector
     */
    double* data();

    /*! Const data method to access to the pointer to the internal vector
     * @return Pointer to the internal vector
     */
    const double* data() const;

    /*! Method to get a pointer to the i-th row of the matrix
     * @return The pointer to the 1st element of the i-th row
     */
    double* getRow(int i);

    /*! Method to get the number of rows of the matrix
     * @return The number of rows
     */
    inline int rows() const
    {
        return n_rows; 
    }

    /*! Method to get the number of columns of the matrix
     * @return The number of columns
     */
    inline int cols() const
    {
        return n_cols;
    }

    /*! Method to get a submatrix selecting a range of rows 
     * @param start_row The index of the first row
     * @param end_row The index of the last row (excluded)
     * @param submatrix The submatrix to be filled
     */
    void getRowSubMatrix(int start_row, int end_row, Matrix& submatrix) const;

    
private:
    std::vector<double> v;
    int n_rows = 0;
    int n_cols = 0;
};


#endif