#ifndef HELPER_HPP
#define HELPER_HPP

// ==================================================================================
// Here I enclose some useful things that I will use in my code
// ==================================================================================

/*!
 * Enumerator to select storage order of the matrix
 */
enum class ORDERING
{
  ROWMAJOR = 0,
  COLUMNMAJOR = 1
};

/*!
 * Enumerator to select the norm
 */
enum class NORM
{
  One = 0,
  Infinity = 1,
  Frobenius = 2 
};

using IndexType = std::array<std::size_t,2>;

/*!
 * Helper function to print a vector
 * @param v The vector to be printed
 * @param name The name of the vector to be displayed
 */
template<class T> void print_vector(const std::vector<T>& v, const std::string& name)
{
   std::cout << name << " = [";
    for(size_t i = 0; i < v.size(); ++i)
    {
        std::cout << v[i];
        if(i < v.size() - 1)
        {
            std::cout << ", ";
        }
        else
        {
            std::cout << "]" << std::endl;
        }
    }
}

#endif