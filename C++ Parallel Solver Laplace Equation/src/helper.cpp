#include "helper.hpp"
#include <cmath>

std::string replace_all(const std::string& input_string, const std::string& replace_word, const std::string& replace_by)
{   
    std::string output_string = input_string;
    // Find the 1st occurrence
    size_t pos = output_string.find(replace_word); 

    // Repeat as long as we have occurrences
    while (pos != output_string.npos) 
    { 
        output_string.replace(pos, replace_word.size(), replace_by); 

        // Find the next occurrence of the substring 
        pos = output_string.find(replace_word, pos + replace_by.size()); 
    } 

    return output_string;
}

double increment_norm(const Matrix& U, const Matrix& U_old)
{
    double sum = 0.0;
    for (int i = 0; i < U.rows(); ++i)
    {
        for (int j = 0; j < U.cols(); ++j)
        {
            sum += (U(i,j) - U_old(i,j))*(U(i,j) - U_old(i,j));
        }
    }

    double h = 1.0 / (U.cols() - 1);
    return std::sqrt(sum*h);
}

void exportToVTK(const std::string& fileName, const Matrix& U, int n, double h)
{
    // open the file
    std::ofstream file(fileName);

    // check if the file was opened
    if (!file.is_open())
    {
        std::cerr << "ERROR: could not open file " << fileName << std::endl;
        return;
    }

    // Write VTK header
    file <<  "# vtk DataFile Version 3.0\n";
    file << "Scalar Field Data\n";
    file << "ASCII\n";                               

    // Write grid data
    file << "DATASET STRUCTURED_POINTS\n";                                // format of the dataset
    file << "DIMENSIONS " << n << " " << n << " " << 1 << "\n";           // number of points in each direction
    file << "ORIGIN 0 0 0\n";                                             // lower-left corner of the structured grid
    file << "SPACING" << " " << h << " " << h << " " << 1 << "\n";        // spacing between points in each direction
    file << "POINT_DATA " << n * n << "\n";                               // number of points

    // Write scalar field data
    file << "SCALARS scalars double\n";               // description of the scalar field
    file << "LOOKUP_TABLE default\n";                 // color table

    // Write vector field data
    for(int j = 0; j < n; ++j)
    {
        for(int i = 0; i < n; ++i)
        {
            file << U(i,j) << "\n";
        }
    }
}

double L2_norm(const Matrix& U, const MuParserInterface::muParserInterface& u_ex, double h)
{
    double sum = 0.0;
    for (int i = 0; i < U.rows(); ++i)
    {
        for (int j = 0; j < U.cols(); ++j)
        {
            sum += (U(i,j) - u_ex(0, i*h, j*h))*(U(i,j) - u_ex(0, i*h, j*h));
        }
    }

    return std::sqrt(sum*h);
}
