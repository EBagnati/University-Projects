#include <mpi.h>
#include "solveJacobi.hpp"

int main(int argc, char **argv)
{   
    
    // Initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Test the performance of the code as the grid size increases (n = 2^k, k = 4,...,8)
    for(int k = 4; k <= 8; ++k)
    {
        int n = std::pow(2,k);
        if(rank == 0)
        {
            std::cout << "===========================================================" << std::endl;
            std::cout << "Grid size: " << n << "x" << n << std::endl;
            std::cout << "===========================================================" << std::endl;
        }
        
        solveJacobi("parameters.txt", n);
    }


    MPI_Finalize();

    
    return 0;
}