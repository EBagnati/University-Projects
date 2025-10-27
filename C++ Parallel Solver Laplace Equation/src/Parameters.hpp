#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include "muParserInterface.hpp"
#include "helper.hpp"
#include "GetPot"


/*!
 * Struct to hold all the parameters of my problem
 */
struct Parameters
{        
    /*!
     * Method to read parameters from a file
     * @param filename The name of the file containing the parameters
     */
    void readFile(const std::string &filename);
    
    double tol;                                 // Tolerance for the convergence criterion
    unsigned max_iter;                          // Maximum number of iterations
    std::string f_string;                       // Forcing term stored as string (to be printed)
    std::string g_string;                       // Dirichlet boundary condition stored as string (to be printed)
    std::string u_ex_string;                    // Exact solution stored as string (to be printed)
    MuParserInterface::muParserInterface f;     // Forcing term
    MuParserInterface::muParserInterface g;     // Dirichlet boundary condition
    MuParserInterface::muParserInterface u_ex;  // Exact solution
    
    /*!
     * Prints all the parameters
     */
    void print() const;
};

#endif