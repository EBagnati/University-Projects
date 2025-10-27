#ifndef SOLVEJACOBI_HPP
#define SOLVEJACOBI_HPP

#include "Matrix.hpp"
#include <string>
#include <mpi.h>
#include "Parameters.hpp"
#include "chrono.hpp"
#include <omp.h>

/*! Method to solve Laplace Equation by Jacobi method
 * @param paramFileName The name of the file with the parameters (that will be read only by rank 0)
 * @param n The number of points in the Cartesian grid
 * @return The matrix with the solution evaluated at each point of the Cartesian grid
 */
Matrix solveJacobi(const std::string& paramFileName, int n);

/*! Method to print results of Jacobi method
 * @param p The struct with the parameters of the problem
 * @param n_iter The number of iterations performed
 * @param U The matrix with the solution evaluated at each point of the Cartesian grid
 * @param h The spacing between points in the Cartesian grid
 * @param chrono The object used to evaluate the elapsed time
 */
void printResults(const Parameters& p, unsigned n_iter, const Matrix& U, double h, const Timings::Chrono& chrono);

#endif
