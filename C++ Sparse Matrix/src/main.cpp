#include "Matrix.hpp"
#include <iostream>
#include "helper.hpp"
#include "chrono.hpp"
#include <random>   

int main()
{
    // Create 2 matrices, one with row major ordering and the other with column major ordering
    algebra::Matrix<double, ORDERING::ROWMAJOR> m_row;
    algebra::Matrix<double, ORDERING::COLUMNMAJOR> m_col;
    
    // Fill them with the values read from file
    m_row.read("lnsp_131.mtx");
    m_col.read("lnsp_131.mtx");

    // Generate a random vector of the proper size with doubles between 0 and 100
    std::default_random_engine gen;
    std::uniform_real_distribution<>dist{0, 100};
    std::vector<double> v(m_row.cols());
    for(std::size_t i = 0; i < v.size(); ++i)
    {
        v[i] = dist(gen);
    }

    Timings::Chrono myChrono;

    // Matrix * vector product: ROW MAJOR ordering - UNCOMPRESSED matrix
    myChrono.start();
    std::vector<double> res = m_row * v;
    myChrono.stop();
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "Matrix * vector product: ROW MAJOR ordering - UNCOMPRESSED matrix" << std::endl;
    std::cout << myChrono << std::endl;

    // Matrix * vector product: COLUMN MAJOR ordering - UNCOMPRESSED matrix
    myChrono.start();
    res = m_col * v;
    myChrono.stop();
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "Matrix * vector product: COLUMN MAJOR ordering - UNCOMPRESSED matrix" << std::endl;
    std::cout << myChrono << std::endl;

    // Matrix * vector product: ROW MAJOR ordering - COMPRESSED matrix
    m_row.compress();
    myChrono.start();
    res = m_row * v;
    myChrono.stop();
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "Matrix * vector product: ROW MAJOR ordering - COMPRESSED matrix" << std::endl;
    std::cout << myChrono << std::endl;

    // Matrix * vector product: COLUMN MAJOR ordering - COMPRESSED matrix
    m_col.compress();
    myChrono.start();
    res = m_col * v;
    myChrono.stop();
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "Matrix * vector product: COLUMN MAJOR ordering - COMPRESSED matrix" << std::endl;
    std::cout << myChrono << std::endl;
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "Matrix norms" << std::endl;
    std::cout << "Norm ONE = " << m_row.norm<NORM::One>() << std::endl;
    std::cout << "Norm INF = " << m_row.norm<NORM::Infinity>() << std::endl;
    std::cout << "Norm FROBENIUS = " << m_row.norm<NORM::Frobenius>() << std::endl;

    return 0;
}