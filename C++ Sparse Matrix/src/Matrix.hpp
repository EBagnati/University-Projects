#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>
#include <iostream>
#include <stdexcept>
#include "helper.hpp"
#include "CompareOperator.hpp"
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>

namespace algebra
{
// Forward declarations to avoid problems with template function operator*
template <class T, ORDERING StorageOrder> class Matrix;
template<class M_Type, class V_Type, ORDERING Order> std::vector<std::common_type_t<M_Type,V_Type>> operator* (Matrix<M_Type, Order> &A, const std::vector<V_Type> &v);

/*!
 * Class that implements a sparse matrix
 * @tparam T The type of the stored elements
 * @tparam StorageOrder The storage order of the elements: row-wise or column-wise
 */
template <class T, ORDERING StorageOrder> class Matrix
{
public:
    /*!
     * Constructor that takes number of rows and columns
     * @param rows Number of rows
     * @param cols Number of columns
     */
    Matrix(std::size_t rows, std::size_t cols) : nrows(rows), ncols(cols) {};

    // Generate also the default constructor (maybe useful to create an empty matrix)
    Matrix() = default;

    /*!
     * Returns the number of rows
     * @return nrows
     */
    inline std::size_t rows() const
    {
      return nrows;
    }

    /*!
     * Returns the number of columns
     * @return ncols
     */
    inline std::size_t cols() const
    {
      return ncols;
    }

    /*!
     * Converts the internal storage from the uncompressed format to CSR or CSC according to StorageOrder
    */
    void compress();

    /*!
     * Converts the internal storage from CSR or CSC to the uncompressed format
    */
    void uncompress();

    /*!
     * Checks if the matrix is compressed or not
     @return TRUE if compressed
    */
    inline bool isCompressed() const
    {
      return compressed;
    }

    /*! Const version of the call operator
     * @param i
     * @param j
     * @return The value
     */
    T operator() (std::size_t i, std::size_t j) const;

    /*! Non const version of the call operator
     * @param i
     * @param j
     * @return The value (which can be changed)
     */
    T& operator() (std::size_t i, std::size_t j);

    // Prints the matrix
    void print();

    /*! Resizes the matrix leaving it in uncompressed state
     * @param rows The new number of rows
     * @param cols The new number of columns
     */
    void resize(std::size_t rows, std::size_t cols);

    // Tell that operator* is a friend of my class Matrix
    template<class M_Type, class V_Type, ORDERING Order> friend std::vector<std::common_type_t<M_Type,V_Type>> operator* (Matrix<M_Type, Order> &A, const std::vector<V_Type> &v);
   
    /*! Reads the matrix written in Matrix Market Format from a file
     * @param fileName The file name where the matrix is stored 
     */
    void read(const std::string& fileName);

    /*! Computes the norm of the matrix
     * @tparam NormType The desired norm (One, Infinity or Frobenius)
     */
    template<NORM NormType> double norm() const;


private:
    // Internal storage for uncompressed version
    std::map<IndexType,T, CompareOperator<StorageOrder>> data;

    // Internal storage for compressed version
    std::vector<std::size_t> inner_idx;
    std::vector<std::size_t> outer_idx;
    std::vector<T> values;

    // Some auxiliary members 
    bool compressed = false;
    std::size_t nrows = 0;
    std::size_t ncols = 0;

    /*
     * Removes any zero element stored in the map 'data' (e.g. added by the user through the call operator).
     * This method works on an uncompressed matrix. If it is called on a compressed one, it does nothing
     */
    void clear_zeros();

};
}

#include "Matrix_imp.hpp"
#endif