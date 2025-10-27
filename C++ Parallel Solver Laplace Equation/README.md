# A matrix–free parallel solver for the Laplace equation

In this project we implement the Jacobi iteration method to solve the Laplace equation on a uniform Cartesian decomposition of the domain (0,1)x(0,1)

The solution is represented as a dense matrix U where each entry U(i,j) is the value of the approximate solution in the node (x_i, y_j)

## Implementation
The code (contained in the `src` folder) is structured as follows:

 - in the file `parameters.txt` the user can set the parameters of the problem, i.e. the tolerance for the convergence criterion, the maximum number of iterations, the forcing term, the exact solution and the Dirichlet boundary condition. Although it would have been more intuitive to assign here also the number of points of the Cartesian grid `n`, I have decided to leave this as a parameter of the `solveJacobi()` function (described later): in this way is easier to change `n` automatically in the main function while performing the test.

 - in `parameters.hpp` and `parameters.cpp` I have defined a struct called `Parameters` to read parameters from file. In particular to read the functions *f*, *g* and *u_ex* I rely on  `MuParserInterface` taken from the Examples of the PACS repository. These 3 functions are read as strings, processed to replace the string *pi* with its value 3.1415... and finally fed to the parser

 - in the files `solveJacobi.hpp` and `solveJacobi.cpp` the Jacobi method is implemented. The function `solveJacobi()` takes as input the parameters of the problem and the number of points `n`, and returns the matrix `U` with the approximate solution.
 Essentially it works as follows: 
    - Rank 0 reads the parameters from file, fills the matrix F with the evaluations of the forcing term in all the points of the cartesian grid and fills the boundary cells of matrix U with the Dirichlet B.C.
    
    -  Then some parameters are broadcasted to the other ranks, whereas matrices F and U are scattered so that each rank gets its `local_F` and `local_U` to process.

    - At this point each rank enters in a loop, that will stop when all ranks have converged or when the maximum number of iterations is reached. First of all each rank exchanges data with its adjacent processors by using `MPI_Sendrecv()` to avoid deadlocks. Then each rank updates its local solution and computes the increment norm w.r.t the previous iteration to check its own local convergence criterion. Finally the informations about local convergence of all the ranks are put together with `MPI_Allreduce()` to verify if we have reached a global convergence. (Unfortunately for this I had to use `char` variables because I didn't manage to use bools with MPI)

    - Once we have reached convergence (or the max number of iterations), rank 0 prints some informations (such as the elapsed time and the L2 error w.r.t. the exact solution) and assembles the global matrix U with the entire solution. Then the solution is exported on a VTK file called `solution.vtk` and located in the `VTK` folder, which can be opened with Paraview

- in the files `Matrix.hpp` and `Matrix.cpp` I have implemented a class to realize a matrix with all elements stored contiguously in a single vector. I have found this useful because using vectors I can resize my matrices at run time, and I have decided to store elements contiguously to make it easier to call MPI functions like `MPI_Scatterv()` and  `MPI_Sendrecv()` (I tried before to use a `std::vector<std::vector<double>>` but I found it much more complicated).

- in the files `helper.hpp` and `helper.cpp` I have enclosed some helper functions, for example those for the computation of norms and the function to export matrices in VTK format

- finally the `main()` test the performance of the code as the grid size increases (n = 2^k, k = 4,...,8)

## Note on OpenMP

Once I have implemented my program using MPI, I have tried to further parallelize the local computations adding an OpenMP directive.
As it can be seen in the code, where I left this directive commented, I have tried to add `#pragma omp parallel for` at line 181, before the loop where each rank updates its local solution. Anyway, performing several tests with 2, 4, or 8 threads and with a different number of MPI ranks, I have noticed that on my machine the elapsed time increases (I have printed some results in `test/data/result_OpenMP.txt` for the case with 4 MPI ranks and 4 OpenMP threads)

## Compilation and Documentation
When we are inside the `src` folder:
 - Typing `make` we can compile the code and get the executable (the code with best performance, i.e. only MPI, will be compiled. To compile with the further OpenMP parallelization, we just need to remove comments from line 181).

 - To run the code with just MPI: `mpiexec -np [MPI_RANKS] ./main` (the user can choose the number of desired ranks replacing [MPI_RANKS])

 - To run the code with also OpenMP (after compiling without comments in line 181): `mpiexec -np [MPI_RANKS] -x OMP_NUM_THREADS=[OPENMP_THREADS] ./main` (the user can write the number of desired ranks and of OpenMP threads in place of [MPI_RANKS] and [OPENMP_THREADS] )

 - Typing `make doc` we can generate the Doxygen documentation (it will be put in the `doc` folder)

 - Since I am using `MuParserInterface` from the PACS repository, in the Makefile I have set `PACS_ROOT` with my path to the PACS Examples

## Testing
Once the program (with MPI or with also OpenMP) is compiled, in the `test` folder we find 2 tests:
 - `scalability_test.sh` executes the `main()` with 1, 2 and 4 MPI ranks (the code should be compiled only with MPI)
 - `scalability_test_OpenMP.sh` executes the `main()` with 1, 2 and 4 MPI ranks and 4 OpenMP threads (the code should be compiled with also OpenMP)

 Before running the scripts, on my machine I had to authorize them with `chmod u+r+x [SCRIPT NAME]`.
 Then we can run the script (for example the 'only MPI' one) with `./scalability_test.sh` (in this case all the outputs will be printed on the terminal) or we can write the outputs in a text file stored in the `data` folder, so that it will be easier to analyze them later, with `./scalability_test.sh > data/result_onlyMPI.txt`

## Discussion of the results

From the file `result_only_MPI.txt` located in `test/data`, we can see that:
 - if we fix the number of ranks and we increase the grid size, as expected the elapsed time and the number of iterations will increase (obviously)

 - if we fix the grid size (for example we consider the grid 128x128) we can see that in most cases the elapsed time decreases as the number of ranks increases (e.g. it goes from 9.08e+05 microsec to 4.36e+05 microsec when doubling the number of ranks from 1 to 2). Anyway we can see that this is not true expecially for small grids (like 16x16) where using 4 processes can slow down the computations (probably because the communication overheads are more significant than the speedup we get parallelizing the computations)

 If we look to the file `result_OpenMP.txt` in the same folder we see that in general the elapsed time for a given number of processes and a given grid size is bigger than the time needed with just MPI. Unfortunately I didn't manage to find an efficient OpenMP implementation 

## Extra
Dirichlet boundary conditions of non homogeneous type are implemented. The user can assign them through the function *g* in the file `parameters.txt` 

