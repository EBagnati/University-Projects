#ifndef HELPER_HPP
#define HELPER_HPP

#include <string>
#include <fstream>
#include <iostream>
#include "Matrix.hpp"
#include "muParserInterface.hpp"

// ================================================================
// Here I enclose some helper functions that I will use in my code
// ================================================================

/*!
 * Helper function to replace all occurrences of a substring in a string (it returns a new string without modifying the input one)
 * @param input_string The string to be processed
 * @param replace_word The substring to be replaced
 * @param replace_by The substring to be put in place of 'replace_word'
 * @return The new string with all occurrences of 'replace_word' replaced by 'replace_by'
 */
std::string replace_all(const std::string& input_string, const std::string& replace_word, const std::string& replace_by);

/*!
 * Helper function to compute the norm of the increment for the convergence criterion
 * @param U The matrix U at the new iteration
 * @param U_old The matrix U at the previous iteration
 * @return The norm of the increment
 */
double increment_norm(const Matrix& U, const Matrix& U_old);

/*!
 * Helper function to export a matrix in VTK format on a file
 * @param fileName The name of the file where to store the matrix
 * @param U The matrix to be stored
 * @param n Number of points along each direction for Cartesian decomposition of the domain
 * @param h Step size
 */
void exportToVTK(const std::string& fileName, const Matrix& U, int n, double h);

/*!
 * Helper function to compute the L2 norm between the exact and the approximate solution
 * @param U The matrix U storing the approximate solution
 * @param u_ex The exact solution
 * @param h Step size
 * @return The L2 norm between the exact and the approximate solution
 */
double L2_norm(const Matrix& U, const MuParserInterface::muParserInterface& u_ex, double h);


#endif