# A sparse matrix

In this project we create a dynamic matrix class template `Matrix<T, StorageOrder>` to handle sparse matrices, and we test the time required to perform a matrix * vector product using different storage orders and compressed or uncompressed formats.

## The class Matrix

The class, since it is a template, is declared in the file `Matrix.hpp` and its methods are implemented in `Matrix_imp.hpp` to improve readibility. The file `helper.hpp` contains some utilities used in the code. More in detail the class has:

- Two template parameters:
    - `T`: the type of values stored in the matrix (e.g. double, float, int, ...) 
    - `StorageOrder`: an enumerator which indicates how elements are stored in memory (row wise or column wise)

- Two constructors: the default one and a constructor which allows to specify the number of rows and the number of columns 

- Two const methods `rows()` and `cols()` which return the number of rows and of columns

- Two methods `compress()` and `uncompress()` which allow to pass from the uncompressed internal storage (COOmap) to the compressed one (CSR or CSC according to `StorageOrder`) and viceversa

- A method `isCompressed()` to test if the matrix is in the compressed format

- Two **const** and **non const** versions of the call operator (`operator()`) coded as requested. If the user tries to do a non admissible operation (e.g. access to an element out of range or add a new element in the compressed format) an exception is thrown. 
    
- A method `print()` to display the matrix

- A method `resize()` to change the matrix dimensions in the uncompressed format. It takes care to delete all the rows/columns that must be discarded if we reduce the matrix to a smaller size.

- A method `read()` to read general real matrices written in *Matrix Market Format* 

- A template method `norm()` to compute the norm of the matrix. The template parameter `NormType` is an enumerator which defines which is the norm to compute (One, Infinity or Frobenius)

In the private part of the class the most relevant members are:

- A map `data` to store matrix in uncompressed format. It relies on the compare operator defined in `CompareOperator.hpp`

- Three vectors `inner_idx` `outer_idx` and `values`to store matrix in compressed format

- A method `clear_zeros()` that, if the matrix is uncompressed, removes any zero element stored (e.g. added by mistake by the user with the call operator). This is employed at the beginning of some other methods, for example to avoid printing or adding to the compressed format some useless zero entries. 
     

Finally we have a friend function `operator*` to perform the product matrix * vector. The matrix and the vector may also store two different types (e.g. a Matrix of double and a vector of std::complex) and the method returns a vector of *std::common_type_t* between these two types.

## Timing
In `main.cpp` we create two matrices (one with row major ordering and the other with column major ordering) reading them from `lnsp_131.mtx` and we multiply them (both in compressed and uncompressed format) by a random vector of proper size. We measure the elapsed time thanks to the utility `chrono.hpp` provided.
As expected, we observe that using the compressed formats is more efficient.

## Compilation and Documentation
When we are inside the `src` folder:
 - Typing `make` we can compile the code and get the executable (then to run the program `./main`)
 - Typing `make doc` we can generate the Doxygen documentation (it will be put in the `doc` folder)



