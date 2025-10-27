#include "solveJacobi.hpp"

Matrix solveJacobi(const std::string& paramFileName, int n)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
   
    // Rank 0 reads the parameters from file and broadcasts them to all the other ranks
    Matrix F;                           // Matrix F to store the forcing term evaluated at each point
    Matrix U;                           // Matrix U to store the solution
    double h;                           // Step size
    std::vector<int> counts_send;       // Number of elements to send to each rank
    std::vector<int> displacements;     // Displacements to apply to the message of each rank
    Parameters p;                       // Struct to store the parameters
    
    if (rank == 0)
    {
        p.readFile("parameters.txt");
        F.resize(n,n);
        U.resize(n,n);

        // Step size
        h = 1.0 / (n - 1); 

        // Fill the matrix F with the forcing term evaluated at each point
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                F(i,j) = p.f(0, i*h, j*h);
            }
        }

        // Fill the boundary cells of U with Dirichlet boundary condition
        for (int i = 0; i < n; ++i)
        {
            U(i,0) = p.g(0, i*h, 0*h);
            U(i,n-1) = p.g(0, i*h, (n-1)*h);
        }

        for (int j = 1; j < n-1; ++j)
        {
            U(0,j) = p.g(0, 0, j*h);
            U(n-1,j) = p.g(0, (n-1)*h, j*h);
        }

        // Prepare the counts and displacements for the scatterv
        counts_send.resize(size);
        displacements.resize(size,0);

        for(int i = 0; i < size; ++i)
        {
            counts_send[i] = (n % size) > i ? (n / size + 1)*n : (n / size)*n;
            if(i > 0)
            {
                displacements[i] = displacements[i-1] + counts_send[i-1];
            }
        }            
    }

    // Broadcast the number of points to each rank
    MPI_Bcast(&(n), 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast the tolerance to each rank
    MPI_Bcast(&(p.tol), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Broadcast the maximum number of iterations to each rank
    MPI_Bcast(&(p.max_iter), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // Each rank gets its local number of rows to process
    int n_local;
    MPI_Scatter(counts_send.data(), 1, MPI_INT, &n_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
    n_local /= n;

    // Each rank gets its local matrix F
    Matrix local_F(n_local, n);
    MPI_Scatterv(F.data(), counts_send.data(), displacements.data(), MPI_DOUBLE, local_F.data(), n_local*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each rank gets its local matrix U (with the boundary cells already filled with the Dirichlet boundary condition)
    Matrix local_U(n_local, n);
    MPI_Scatterv(U.data(), counts_send.data(), displacements.data(), MPI_DOUBLE, local_U.data(), n_local*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each rank defines:
    // - a matrix local_U to store the local solution at the previous iteration
    // - a matrix comp_U to compute the local solution at the current iteration, taking into account also the rows of the adjacent ranks
    // - a matrix comp_U_old to store the local solution at the previous iteration, taking into account also the rows of the adjacent ranks
    Matrix comp_U;
    Matrix comp_U_old;
    Matrix local_U_old = local_U;

    // Rank 0 initializes comp_U and comp_U_old with the value of local_U
    if(rank == 0)
    {
        // Rank 0 leaves the last row of comp_U free to store the row coming from rank 1
        comp_U.resize(n_local+1, n);
        comp_U_old.resize(n_local+1, n);

        for(int i = 0; i < n_local; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                comp_U(i,j) = local_U(i,j);
                comp_U_old(i,j) = local_U(i,j);
            }
        }
    }

    else if (rank == size-1)
    {
        // Rank 'size-1' leaves the first row of comp_U free to store the row coming from rank 'size-2'
        comp_U.resize(n_local+1, n);
        comp_U_old.resize(n_local+1, n);

        for(int i = 1; i < n_local+1; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                comp_U(i,j) = local_U(i-1,j);
                comp_U_old(i,j) = local_U(i-1,j);
            }
        }
    }
    
    else
    {
        // All the other ranks leave the first and last row of comp_U free to store the rows coming from the adjacent ranks
        comp_U.resize(n_local+2, n);
        comp_U_old.resize(n_local+2, n);

        for(int i = 1; i < n_local+1; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                comp_U(i,j) = local_U(i-1,j);
                comp_U_old(i,j) = local_U(i-1,j);
            }
        }
    }

    char local_convergence = 0;     // Flag which becomes != 0 when the current rank has converged (I use char to be able to use MPI datatypes)
    char global_convergence = 0;    // Flag which becomes != 0 when ALL the ranks have converged (I use char to be able to use MPI datatypes)    
    h = 1.0 / (n - 1);

    // Measure the time elapsed to solve the problem
    Timings::Chrono chrono;
    chrono.start();
    unsigned k;

    for(k = 1; k < p.max_iter && global_convergence < size; ++k)
    {
        if(size > 1)
        {
            // Each rank exchanges the rows it needs with adjacent ranks
            if(rank == 0)
            {
                // Rank 0 needs to exchange its last row with the first row of rank 1                
                MPI_Sendrecv(comp_U.getRow(n_local-1), n, MPI_DOUBLE, 1, 0, comp_U_old.getRow(n_local), n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            else if(rank == size-1)
            {
                // The last rank needs to exchange its first row with the last row of rank size-2
                MPI_Sendrecv(comp_U.getRow(1), n, MPI_DOUBLE, size-2, 0, comp_U_old.getRow(0), n, MPI_DOUBLE, size-2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            else
            {
                // All other ranks need to exchange their first row with the last row of rank 'rank-1'
                // and their last row with the first row of rank 'rank+1'
                MPI_Sendrecv(comp_U.getRow(1), n, MPI_DOUBLE, rank-1, 0, comp_U_old.getRow(0), n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Sendrecv(comp_U.getRow(n_local), n, MPI_DOUBLE, rank+1, 0, comp_U_old.getRow(n_local+1), n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // The last row to process according to the rank
        int last_row = (rank > 0 && rank < size-1) ? n_local+1 : n_local;
        // An offset to consider the correct row in the local_F matrix
        int offset = (rank == 0) ? 0 : 1;
        
        // Each rank updates the local solution
        // UNCOMMENT THE FOLLOWING LINE TO PARALLELIZE THE LOOP WITH OPENMP
        //#pragma omp parallel for
        for(int i = 1; i < last_row; ++i)
        {
            for(int j = 1; j < n-1; ++j)
            {
                comp_U(i,j) = 0.25*(comp_U_old(i-1,j) + comp_U_old(i+1,j) + comp_U_old(i,j-1) + comp_U_old(i,j+1) + h*h*local_F(i-offset,j));

            }
        }

        // Each rank extracts the local solution from comp_U (ignoring the rows used to exchange information with adjacent ranks)
        if(rank == 0)
        {
            comp_U.getRowSubMatrix(0, n_local, local_U);
        }
            
        else
        {
            comp_U.getRowSubMatrix(1, n_local+1, local_U);
        }

        // Each rank checks its own local convergence criterion
        increment_norm(local_U, local_U_old) < p.tol ? local_convergence = 1 : local_convergence = 0;       

        // Each rank checks if we have convergence in all the ranks
        MPI_Allreduce(&local_convergence, &global_convergence, 1, MPI_CHAR, MPI_SUM, MPI_COMM_WORLD);

        // Update local_U_old
        local_U_old = local_U;   

        // $$ NEW $$
        comp_U_old = comp_U;
    }

    chrono.stop();

    // Rank 0 assemble the final matrix U
    MPI_Gatherv(local_U.data(), n_local*n, MPI_DOUBLE, U.data(), counts_send.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        // Rank 0 prints the results 
        printResults(p, k, U, h, chrono);
        // And writes the solution in VTK format on a file
        exportToVTK("../VTK/solution.vtk", U, n, h);
    }

    return U;
}


void printResults(const Parameters& p, unsigned n_iter, const Matrix& U, double h, const Timings::Chrono& chrono)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Solution computed with " << size << " rank(s) in " << n_iter << " iterations" << std::endl;

    if(n_iter >= p.max_iter)
    {
        std::cout << "WARNING: The loop ended after reaching maximum number of iterations. The solution may not be correct" << std::endl;
    }

    std::cout << "Error between approximate and exact solution (L2 norm): " << L2_norm(U, p.u_ex, h) << std::endl;

    std::cout << "Time elapsed (on rank 0). "<< std::scientific << chrono.wallTime() << " microsec" << std::endl;
  
}