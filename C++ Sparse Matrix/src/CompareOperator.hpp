#ifndef COMPARE_OPERATOR_HPP
#define COMPARE_OPERATOR_HPP

// ====================================================================================================================
// Here I define the compare operator to be used in the map depending on the requested ordering
// ====================================================================================================================

#include "helper.hpp"

/*!
 * Template functor that implements the comparison operator according to the selected storage order
 * @tparam StorageOrder The storage order of the elements: row-wise or column-wise
 */
template<ORDERING StorageOrder> struct CompareOperator
{
    bool operator() (const IndexType& lhs, const IndexType& rhs) const
    {
        if constexpr(StorageOrder == ORDERING::ROWMAJOR)
        {
            // The defaulted comparison operator is fine
            return lhs < rhs; 
        }

        else
        {
            return (lhs[1] < rhs [1]) || (lhs[1] == rhs[1] && lhs[0] < rhs[0]);
        }
    }
};

#endif