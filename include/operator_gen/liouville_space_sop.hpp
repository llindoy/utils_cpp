#ifndef UTILS_OPERATOR_GEN_LIOUVILLE_SPACE_SOP_HPP
#define UTILS_OPERATOR_GEN_LIOUVILLE_SPACE_SOP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>

#include "../timing_macro.hpp"
#include "sSOP.hpp"

namespace utils
{

//TODO: Come up with more efficient structure for handling JW strings.  
class liouvilleSpace
{
public:
    template <typename T>
    static inline sSOP<T> commutator(const sSOP<T>& sop)
    {
        //generate a commutator sop given the standard sop
        sSOP<T> csop;
        for(const auto& pop : sop)
        {
            csop += pop;
            csop -= tildefy(pop);
        }
        return csop;
    }

    template <typename T>
    static inline sSOP<T> anticommutator(const sSOP<T>& sop)
    {
        //generate a commutator sop given the standard sop
        sSOP<T> csop;
        for(const auto& pop : sop)
        {
            csop += pop;
            csop += tildefy(pop);
        }
        return csop;
    }

    template <typename T>
    static inline sSOP<T> left(const sSOP<T>& sop)
    {
        return sop;
    }

    template <typename T>
    static inline sSOP<T> right(const sSOP<T>& sop)
    {
        //generate a commutator sop given the standard sop
        sSOP<T> csop;
        for(const auto& pop : sop)
        {
            csop += tildefy(pop);
        }
        return csop;
    }

    template <typename T>
    static inline sSOP<T> outer(const sSOP<T>& lsop, const sSOP<T>& rsop)
    {
        //generate a commutator sop given the standard sop
        sSOP<T> csop;
        for(const auto& rpop : rsop)
        {
            auto rtpop = tildefy(rpop);
            for(const auto& lpop : lsop)
            {
                csop += lpop*rtpop;
            }
        }
        return csop;
    }

protected:
    sNBO<T> tildefy(const sNBO<T>& pop)
    {
        sNBO<T> term;    term.coeff() = pop.coeff();

        for(const auto& op : pop)
        {
            term *= sOP({op.op()+std::string("~"), op.mode()});
        }
        return term;
    }

};



}   //namespace utils

#endif

