#ifndef UTILS_COMPOSITE_OPERATOR_HPP
#define UTILS_COMPOSITE_OPERATOR_HPP

#include "sSOP.hpp"

namespace utils
{

template <typename T>
class composite_operator
{
public:
    composite_operator(){}

protected:
    //the different coefficients in the sum of these terms
    std::vector<T> m_coeff;

    //the physical modes that this is acting on
    std::vector<size_t> m_modes;

    //any child composite operators that this mode is constructed from
    std::vector<composite_operator<T> > m_lcoeff;


    //the operators that this object actually contains.  This will only ever be referenced if m_modes.size() == 1, in which case the composite operator is a primitive operator
    std::vector<sOP> m_ops;
};


template <typename T>
class compressedSOP
{

};

}//namespace utils

#endif

