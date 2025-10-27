#ifndef MATRIX_IMP_HPP
#define MATRIX_IMP_HPP

// ====================================================================================================================
// To improve readibility, I enclose in this file the definitions of the members of the template class Matrix
// ====================================================================================================================

#include "Matrix.hpp"
using algebra::Matrix;

// Compress method
template<class T, ORDERING StorageOrder> void Matrix<T, StorageOrder>::compress()
{
    if(isCompressed())
    {
        // The matrix is already compressed, I have nothing to do
        return; 
    }

    // Convert to CSR or CSC

    // First of all I clear any zero element accidentally inserted in the matrix
    clear_zeros();

    // Number of non zero stored elements
    std::size_t nnz = data.size();      

    if(nnz == 0)
    {
        // The matrix is empty, I have nothing to do
        compressed = true;
        return;
    }
    
    // Resize the vectors to the proper size: if I have Row Major ordering I consider nrows, otherwise ncols
    std::size_t nlines;
    if constexpr(StorageOrder == ORDERING::ROWMAJOR)
    {
        nlines = nrows;
    }

    else
    {
        nlines = ncols;
    }

    inner_idx.resize(nlines + 1);
    outer_idx.resize(nnz);
    values.resize(nnz);

    // Assign 1st value to inner_idx
    inner_idx[0] = 0;

    unsigned int count = 0;
    for(std::size_t i = 0; i < nlines; ++i)
    {
        // line_begin and line_end refer to the 1st and last entries of the current line 'i' (which can be a row or a column). 
        // These entries may be == 0 (so not stored in the map) or != 0 (so they are in the map)
        IndexType line_begin;
        IndexType line_end;

        if constexpr(StorageOrder == ORDERING::ROWMAJOR)
        {
            line_begin = {i, 0};
            line_end = {i, ncols - 1};
        }

        else
        {
            line_begin = {0, i};
            line_end = {nrows - 1, i};
        }

        // it1 = iterator to the 1st element != 0 of the current line (if it exists) or to the 1st element != 0 of the first non-empty next line (if it exists). 
        // If there isn't any element != 0 left, it will be == cend() 
        auto it1 = data.lower_bound(line_begin);
        // it2 = iterator to the 1st element != 0 of the first non-empty next line (if it exists). Otherwise it will be == cend()
        auto it2 = data.upper_bound(line_end);

        if(it1 == data.end())
        {
            // No more elements != 0
            break;
        }

        else if(it1 == it2)
        {
            // No elements != 0 in the current line
            inner_idx[i+1] = count;
        }

        else
        {
            // While it1 != it2 we have non zero elements in the current line
            while(it1 != it2)
            {
                values[count] = it1->second;        // Store the value
                
                if constexpr(StorageOrder == ORDERING::ROWMAJOR)
                {
                    outer_idx[count] = it1->first[1];   // Store column index
                }

                else
                {
                    outer_idx[count] = it1->first[0];   // Store row index
                }

                ++count;
                ++it1;
            }

            inner_idx[i+1] = count;
        }
    }

    // Once the compression is done, I clear the map to save memory
    data.clear();
    // And I set compressed to true
    compressed = true;
}

// Uncompress method
template<class T, ORDERING StorageOrder> void Matrix<T, StorageOrder>::uncompress()
{
    if(!isCompressed())
    {
        // The matrix is already uncompressed, I have nothing to do
        return;
    }

    // Convert to uncompressed format

    // Number of non zero stored elements
    std::size_t nnz = values.size();

    if(nnz == 0)
    {
        // The matrix is empty, I have nothing to do
        compressed = false;
        return;
    }

    // Number of rows (or of columns)
    std::size_t nlines = inner_idx.size() - 1;

    // Loop over all the rows (or all the cols)
    for(std::size_t i = 0; i < nlines; ++i)
    {
        // For each line, the elements are stored in the interval [ inner(i), inner(i+1) )
        for(std::size_t k = inner_idx[i]; k < inner_idx[i+1]; ++k)
        {
            IndexType idx;
            if constexpr (StorageOrder == ORDERING::ROWMAJOR)
            {
                idx = {i, outer_idx[k]};
            }

            else
            {
                idx = {outer_idx[k], i};
            }

            data.emplace(idx, values[k]);
        }
    }

    // Once the uncompression is done, I clear the vectors to save memory
    inner_idx.clear();
    outer_idx.clear();
    values.clear();
    // And I set compressed to false
    compressed = false;
}

// Const version of the call operator
template<class T, ORDERING StorageOrder> T Matrix<T,StorageOrder>::operator() (std::size_t i, std::size_t j) const
{    
    if(i >= nrows || j >= ncols)
    {
        // ERROR: at least one of the indeces is out of range -> I throw an exception        
        throw std::out_of_range("Indeces out of range");
    }

    if(!isCompressed())
    {
        // If not compressed, I can look for the element in the map
        auto it = data.find({i,j});
        if(it == data.end())
        {
            // Element not found
            return 0;
        }

        // Return the element
        return it->second;
    }

    else
    {
        // The matrix is compressed

        if constexpr(StorageOrder == ORDERING::COLUMNMAJOR)
        {
            // If the matrix is stored by columns, I just need to swap the indeces and then proceed as if it were stored by rows
            std::swap(i,j);
        }
            

        // Iterators to outer_idx[inner_idx[i]] and outer_idx[inner_idx[i+1]], i.e. the range in outer_idx where I could find 'j' if the element is stored
        const auto it1 = outer_idx.cbegin() + inner_idx[i];
        const auto it2 = outer_idx.cbegin() + inner_idx[i + 1];

        // I look for j in that range of outer_idx
        const auto res = std::find(it1, it2, j);

        if(res == it2)
        {
            // Element not found
            return 0;
        }

        // Otherwise that element exists and I can return its value. To get its position in 'values' I compute the distance between 'res' and the begin of 'outer_idx'
        return values[std::distance(outer_idx.cbegin(), res)];
        
    }
    
}

// Non const version of the call operator
template<class T, ORDERING StorageOrder> T& Matrix<T,StorageOrder>::operator() (std::size_t i, std::size_t j) 
{
    if(i >= nrows || j >= ncols)
    {
        // ERROR: at least one of the indeces is out of range -> I throw an exception        
        throw std::out_of_range("Indeces out of range");
    }

    if(!isCompressed())
    {
        // If not compressed, I can look for the element in the map
        auto it = data.find({i,j});
        if(it == data.end())
        {
            // Element not found -> I add it (with value 0)
            IndexType indeces = {i,j};
            auto ret = data.insert(std::make_pair(indeces, 0));
            // And I return the iterator to the element (so that it can be modified)
            return ret.first->second;
        }

        // Otherwise the element is already present
        return it->second;
    }

    else
    {
        // The matrix is compressed 

        if constexpr(StorageOrder == ORDERING::COLUMNMAJOR)
        {
            // If the matrix is stored by columns, I just need to swap the indeces and then proceed as if it were stored by rows
            std::swap(i,j);
        }

        // Iterators to outer_idx[inner_idx[i]] and outer_idx[inner_idx[i+1]], i.e. the range in outer_idx where I could find 'j' if the element is stored
        const auto it1 = outer_idx.cbegin() + inner_idx[i];
        const auto it2 = outer_idx.cbegin() + inner_idx[i + 1];

        // I look for j in that range of outer_idx
        const auto res = std::find(it1, it2, j);

        if(res == it2)
        {
            // Element not found -> ERROR: the user is trying to add a new element while the matrix is compressed
            throw std::invalid_argument("You cannot add new elements while the matrix is compressed");
        }

        // Otherwise that element exists and I can return its value. To get its position in 'values' I compute the distance between 'res' and the begin of 'outer_idx'
        return values[std::distance(outer_idx.cbegin(), res)];  
    }
}

// Clear zeros method
template<class T, ORDERING StorageOrder> void Matrix<T,StorageOrder>::clear_zeros() 
{
    if(isCompressed())
    {
        // I don't do anything because the matrix is compressed 
        return;
    }

    // Looks for the first zero element in the map (if any)
    auto it = std::find_if(data.begin(), data.end(), [](auto const& i) {return i.second == static_cast<T>(0);});
    while(it != data.end())
    {
        // Erase that element
        data.erase(it);
        // Looks for the next zero element (if any)
        it = std::find_if(data.begin(), data.end(), [](auto const& i) {return i.second == static_cast<T>(0);});
    }
}

// Print method
template<class T, ORDERING StorageOrder> void Matrix<T,StorageOrder>::print() 
{
    if(!isCompressed())
    {

        // First of all I clear any zero element accidentally inserted in the matrix
        clear_zeros();

        for(auto it = data.cbegin(); it != data.cend(); ++it)
        {
            std::cout << "[" << it->first[0] << ", " << it->first[1] << "] " << it->second << "\n";
        }
    }

    else
    {
        print_vector(values, "Values");
        print_vector(inner_idx, "Inner indeces");
        print_vector(outer_idx, "Outer indeces");
    }
}

// Resize method
template<class T, ORDERING StorageOrder> void Matrix<T,StorageOrder>::resize(std::size_t rows, std::size_t cols)
{
    // If the matrix is compressed, I uncompress it
    if(isCompressed())
    {
        uncompress();
    }

    if(StorageOrder == ORDERING::ROWMAJOR)
    {
        // If I am reducing the number of rows, I need to delete some rows 
        if(nrows > rows)
        {   
            // With row major ordering I can access a whole row easily
            for(std::size_t r = rows; r < nrows; ++r)
            {
                IndexType row_begin = {r, 0};
                IndexType row_end = {r, ncols - 1};

                auto it1 = data.lower_bound(row_begin);
                auto it2 = data.upper_bound(row_end);

                data.erase(it1, it2);
            }
        }

        // Then if I am reducing the number of cols, I need to delete some cols
        if(ncols > cols)
        {
            // With Row Major ordering I have no special access to the columns so I proceed in this way
            for(std::size_t r = 0; r < nrows; ++r)
            {
                for(std::size_t c = cols; c < ncols; ++c)
                {
                    data.erase({r,c});
                }
            }
        }

    }

    else
    {
        if(ncols > cols)
        {
            // With column major ordering I can access a whole column easily
            for(std::size_t c = cols; c < ncols; ++c)
            {
                IndexType col_begin = {0, c};
                IndexType col_end = {nrows - 1, c};

                auto it1 = data.lower_bound(col_begin);
                auto it2 = data.upper_bound(col_end);

                data.erase(it1, it2);
            }
        }

        if(nrows > rows)
        {
            // With column Major ordering I have no special access to the rows so I proceed in this way
            for(std::size_t c = 0; c < ncols; ++c)
            {
                for(std::size_t r = rows; r < nrows; ++r)    
                {
                    data.erase({r,c});
                }
            }
        }
    }

    nrows = rows;
    ncols = cols;
} 

/*! Multiplication between a matrix (A) and a vector (v)
 * @param v The vector
 * @return The result of A*v
 */
template<class M_Type, class V_Type, ORDERING Order> std::vector<std::common_type_t<M_Type,V_Type>> algebra::operator* (Matrix<M_Type, Order> &A, const std::vector<V_Type> &v)
{
    // Check if the dimensions are correct
    if(A.ncols != v.size())
    {
        // ERROR: at least one of the indeces is out of range -> I throw an exception
        throw std::invalid_argument("Incompatible dimensions");
    }    
    
    // I clear any zero element accidentally inserted in the matrix to avoid useless computations (if the matrix is compressed this won't do anything)
    A.clear_zeros();

    // The vector where I will store the result. It must be of the common type between the type of elements stored in the vector (V_Type) and of those stored in the matrix (M_Type)
    std::vector<std::common_type_t<M_Type,V_Type>> res(A.nrows);

    if(A.compressed)
    {
        if constexpr(Order == ORDERING::ROWMAJOR)
        {
            // The matrix is in CSR storage -> the most efficient way is the classical algorithm 'row-times-vector' with a double loop
            for(std::size_t i = 0; i < A.nrows; ++i)
            {
                // In CSR format, the elements of the i-th row correspond to the positions [inner_idx[i], inner_idx[i+1]) of values and outer_idx 
                std::size_t row_begin = A.inner_idx[i];
                std::size_t row_end = A.inner_idx[i+1];

                for(std::size_t j = row_begin; j < row_end; ++j)
                {
                    // I multiply only the non-zero elements of A with the corresponding elements of v
                    res[i] += static_cast<std::common_type_t<M_Type, V_Type>>(A.values[j]) * static_cast<std::common_type_t<M_Type, V_Type>>(v[A.outer_idx[j]]);
                }
            }

        }

        else
        {
            // The matrix is in CSC storage -> the most efficient way is to do a linear combination of the columns
            for(std::size_t j = 0; j < A.ncols; ++j)
            {
                // In CSC format, the elements of the i-th col correspond to the positions [inner_idx[i], inner_idx[i+1]) of values and outer_idx 
                std::size_t col_begin = A.inner_idx[j];
                std::size_t col_end = A.inner_idx[j+1];

                for(std::size_t i = col_begin; i < col_end; ++i)
                {
                    // I update the proper component of the result 
                    res[A.outer_idx[i]] += static_cast<std::common_type_t<M_Type, V_Type>>(v[j]) * static_cast<std::common_type_t<M_Type, V_Type>>(A.values[i]);
                }
            }

        }

    }

    else
    {
        // The matrix is not compressed -> traverse the map (both when it is stored row wise or column wise)
        for(auto it = A.data.cbegin(); it != A.data.cend(); ++it)
        {
            res[it->first[0]] += static_cast<std::common_type_t<M_Type, V_Type>>(v[it->first[1]]) * static_cast<std::common_type_t<M_Type, V_Type>>(it->second);
        }
    }

    return res;
}

template<class T, ORDERING StorageOrder> void Matrix<T,StorageOrder>::read(const std::string& fileName)
{
    std::ifstream fileStream(fileName);

    if(fileStream.fail())
    {
        // The file has not been opened correclty
        std::cerr << "ERROR: Cannot open the file " << fileName << std::endl;
        return;
    }
    
    // Reads the 1st line of the file
    std::string line;
    getline(fileStream, line);

    if(line != "%%MatrixMarket matrix coordinate real general")
    {
        // The matrix is not in the correct format (general real matrices in matrix market format)
        std::cerr << "ERROR: Matrix not in the correct format" << std::endl;
        std::cout << line << std::endl;
        return;
    }

    if(isCompressed())
    {
        // If the matrix is compressed, I uncompress it to be able to add easily the new elements
        uncompress();
    }


    // I discard the initial lines which contain comments (I recognize them because they start with '%')
    while(line.starts_with('%'))
    {
        getline(fileStream, line);
    }

    // According to Matrix Market Format, the 1st line without '%' is like "M N K" where M = nrows, N = ncols, K = number of non zero elements
    std::istringstream iss{line};
    std::size_t M, N, K;
    iss >> M >> N >> K;

    // I clear my matrix (in case it was already filled with some values)
    data.clear();

    // I resize the matrix to the proper size
    resize(M,N);

    for(std::size_t i = 0; i < K; ++i)
    {
        // I read all the K following lines with the non zero elements: they are like "r c val" where 'r' is the row, 'c' the column and 'val' the value of the current element 
        getline(fileStream, line);
        std::istringstream iss{line};
        std::size_t r,c;
        T val;
        iss >> r >> c >> val;
        // I add the element in the matrix: since in Matrix Market format indeces starts from 1, I need to subtract 1 both to 'r' and 'c'
        operator()(r-1,c-1) = val;
    }

    return;
} 

template<class T, ORDERING StorageOrder>
template<NORM NormType> double Matrix<T,StorageOrder>::norm() const
{
    // The result to be returned
    double res = 0.;

    if constexpr(NormType == NORM::One)
    {
        // For the ONE norm, we must compute the sum of all the columns and take the maximum

        // Create a vector of ncols component to store the sum of the abs values of each column
        std::vector<double> cols_sum(ncols);

        if constexpr(StorageOrder == ORDERING::ROWMAJOR)
        {
            if(isCompressed())
            {
                // CSR format
                // =================================================================
                for(std::size_t i = 0; i < values.size(); ++i)
                {
                    // I exploit the fact that the element 'values[i]' is stored in 'outer_idx[i]' column
                    cols_sum[outer_idx[i]] += std::abs(values[i]);
                }
            }

            else
            {
                // Uncompressed - Row Major matrix
                // =================================================================
                for(auto it = data.cbegin(); it != data.cend(); ++it)
                {
                    cols_sum[it->first[1]] += std::abs(it->second);
                }
            }
        }

        else if constexpr(StorageOrder == ORDERING::COLUMNMAJOR)
        {
            if(isCompressed())
            {
                // CSC format
                // =================================================================

                // Loop over all the columns
                for(std::size_t j = 0; j < inner_idx.size() - 1; ++j)
                {
                    // Loop over all the elements of that column
                    for(std::size_t k = inner_idx[j]; k < inner_idx[j+1]; ++k)
                    {
                        cols_sum[j] += std::abs(values[k]);
                    }
                }
            }

            else
            {
                // Uncompressed - Column Major matrix
                // =================================================================
                for(auto it = data.cbegin(); it != data.cend(); ++it)
                {
                    cols_sum[it->first[1]] += std::abs(it->second);
                }
            }
        }
   
        auto max_it = std::max_element(cols_sum.begin(), cols_sum.end());
        // Checks if the vector is non empty (otherwise the norm is 0)
        if(max_it != cols_sum.end())
        {
            res = *max_it;
        }
    
        return res;
    }

    else if constexpr(NormType == NORM::Infinity)
    {
        // For the INFINITY norm, we must compute the sum of all the rows and take the maximum

        // Create a vector of nrows component to store the sum of the abs values of each row
        std::vector<double> rows_sum(nrows);

        if constexpr(StorageOrder == ORDERING::ROWMAJOR)
        {
            if(isCompressed())
            {
                // CSR format
                // =================================================================

                // Loop over all the rows
                for(std::size_t i = 0; i < inner_idx.size(); ++i)
                {
                    // Loop over all the elements of that row
                    for(std::size_t j = inner_idx[i]; j < inner_idx[j+1]; ++j)
                    {
                        rows_sum[i] += std::abs(values[j]);
                    }
                }
            }

            else
            {
                // Uncompressed - Row Major matrix
                // =================================================================
                for(auto it = data.cbegin(); it != data.cend(); ++it)
                {
                    rows_sum[it->first[0]] += std::abs(it->second);
                }
            }
        }

        else if constexpr(StorageOrder == ORDERING::COLUMNMAJOR)
        {
            if(isCompressed())
            {
                // CSC format
                // =================================================================

                for(std::size_t i = 0; i < values.size(); ++i)
                {
                    // I exploit the fact that the element 'values[i]' is stored in 'outer_idx[i]' row
                    rows_sum[outer_idx[i]] += std::abs(values[i]);
                }
            }

            else
            {
                // Uncompressed - Column Major matrix
                // =================================================================
                for(auto it = data.cbegin(); it != data.cend(); ++it)
                {
                    rows_sum[it->first[0]] += std::abs(it->second);
                }
            }
        }

        auto max_it = std::max_element(rows_sum.begin(), rows_sum.end());
        // Checks if the vector is non empty (otherwise the norm is 0)
        if(max_it != rows_sum.end())
        {
            res = *max_it;
        }
        return res;
    }

    else
    {
        // For the FROBENIUS norm, we must compute the sum of all the entries squared and take the square root

        // Variable to store the sum of all elements squared
        double sum = 0.;

        if(isCompressed())
        {
            for(auto val : values)
            {
                sum += pow(std::abs(val), 2);
            }
        }

        else
        {
            for(auto it = data.cbegin(); it != data.cend(); ++it)
            {
                sum +=  pow(std::abs(it->second), 2);
            }
        }
        
        return sqrt(sum);
    }
}

#endif
