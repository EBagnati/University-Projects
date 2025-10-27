#include "Parameters.hpp"

void Parameters::readFile(const std::string& filename)
{
    GetPot file(filename.c_str());

    tol = file("tol", 1.e-6);
    max_iter = file("max_iter", 100000);
    f_string = file("f", "0");    
    g_string = file("g", "0");
    u_ex_string = file("u_ex", "0");

    // Before feeding 'f', 'g' and 'u_ex' to set_expression, replaces all occurences of 'pi' with its value
    f.set_expression(replace_all(f_string, "pi", std::to_string(M_PI)));
    g.set_expression(replace_all(g_string, "pi", std::to_string(M_PI)));
    u_ex.set_expression(replace_all(u_ex_string, "pi", std::to_string(M_PI)));
}

void Parameters::print() const
{
    std::cout << "==============================" << std::endl;
    std::cout << "PARAMETERS" << std::endl;
    std::cout << "==============================" << std::endl;

    std::cout << "tol = " << tol << std::endl;
    std::cout << "max_iter = " << max_iter << std::endl;
    std::cout << "f = " << f_string << std::endl;
    std::cout << "g = " << g_string << std::endl;
    std::cout << "u_ex = " << u_ex_string << std::endl;
}