#ifndef UTILS_OPERATOR_GEN_COMPRESSED_SOP_HPP
#define UTILS_OPERATOR_GEN_COMPRESSED_SOP_HPP

#include <linalg/linalg.hpp>

#include <memory>
#include <list>
#include <vector>
#include <algorithm>

#include <linalg/linalg.hpp>

#include "sSOP.hpp"

#ifdef CEREAL_LIBRARY_FOUND
#include <cereal/types/vector.hpp>
#endif

namespace utils
{

class opTerms
{
public:
    using op_type = sOP;
    using pointer_type = std::shared_ptr<op_type>;

    using container_type = std::vector<size_t>;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using reverse_iterator = typename container_type::reverse_iterator;
    using const_reverse_iterator = typename container_type::const_reverse_iterator;

public:
    opTerms(){}
    opTerms(size_t nterms)
    {
        try {m_r.reserve(nterms);}  
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct opTerms object.");
        }
    }
    opTerms(const container_type& r) 
    try : m_r(r) {CALL_AND_HANDLE(std::sort(m_r.begin(), m_r.end()), "Failed to sort r array.");}
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to construct opTerms object.");
    }
        
    opTerms(container_type&& r) 
    try : m_r(std::move(r)) {CALL_AND_HANDLE(std::sort(m_r.begin(), m_r.end()), "Failed to sort r array.");}
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to construct opTerms object.");
    }

    opTerms(const opTerms& o) = default;
    opTerms(opTerms&& o) = default;

    opTerms& operator=(const opTerms& o) = default;
    opTerms& operator=(opTerms&& o) = default;

    const container_type& r() const{return m_r;}
    size_t nterms() const{return m_r.size();}

    void append_indices(const container_type& _r)
    {
        try
        {
            for(auto& r : _r){m_r.push_back(r);}
            std::sort(m_r.begin(), m_r.end());
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to append indices.");
        }
    }

    iterator begin() {  return iterator(m_r.begin());  }
    iterator end() {  return iterator(m_r.end());  }
    const_iterator begin() const {  return const_iterator(m_r.begin());  }
    const_iterator end() const {  return const_iterator(m_r.end());  }

    reverse_iterator rbegin() {  return reverse_iterator(m_r.rbegin());  }
    reverse_iterator rend() {  return reverse_iterator(m_r.rend());  }
    const_reverse_iterator rbegin() const {  return const_reverse_iterator(m_r.rbegin());  }
    const_reverse_iterator rend() const {  return const_reverse_iterator(m_r.rend());  }

    sOP& op(){return m_op;}
    const sOP& op() const{return m_op;}

    void clear()
    {
        m_r.clear();
        m_op.clear();
    }

    bool contains_index(size_t r) const{return (std::find(m_r.begin(), m_r.end(), r) != m_r.end());}
protected:
    container_type m_r;
    sOP m_op;

#ifdef CEREAL_LIBRARY_FOUND
public:
    template <typename archive>
    void serialize(archive& ar) 
    {
    }
#endif
};

//a generic sum of product operator object.
template <typename T>
class compressedSOP
{
public:
    using element_type = opTerms;
    using op_type = opTerms::op_type;

    using element_container_type = typename element_type::container_type;

    using mode_terms_type = std::vector<element_type>;
    using container_type = std::vector<mode_terms_type>;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using reverse_iterator = typename container_type::reverse_iterator;
    using const_reverse_iterator = typename container_type::const_reverse_iterator;

    compressedSOP() : m_nterms(0){}
    compressedSOP(size_t nterms, size_t nmodes) 
    try : m_mode_operators(nmodes), m_nterms(nterms), m_coeff(nterms)
    {
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to construct sum of product operator object.");
    }    


    compressedSOP(const compressedSOP& o) = default;
    compressedSOP(compressedSOP&& o) = default;

    compressedSOP& operator=(const compressedSOP& o) = default;
    compressedSOP& operator=(compressedSOP&& o) = default;

    void resize(size_t nterms)
    {
        try
        {
            clear();
            m_coeff.resize(nterms); std::fill(m_coeff.begin(), m_coeff.end(), T(1));
            m_nterms = nterms;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to resize sp hamiltonian object.");
        }
    }
   
    void bind(const opTerms& op, size_t nu) 
    {
        //first we check that non of the r-indices in this object have already been bound
        ASSERT(valid_rvals(op.r(), nu), "Unable to bind operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");

        //we add in the operator if it isn't 
        CALL_AND_HANDLE(m_mode_operators[nu].push_back(op), "Failed to push mode operator term to list.");
    }

    void bind(opTerms&& op, size_t nu) 
    {
        ASSERT(valid_rvals(op.r(), nu), "Unable to bind operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");
        CALL_AND_HANDLE(m_mode_operators[nu].push_back(std::move(op)), "Failed to push mode operator term to list.");
    }

    void bind(const sOP& op, const element_container_type& r, size_t nu)
    {
        //first we check that non of the r-indices in this object have already been bound
        ASSERT(valid_rvals(r, nu), "Unable to bind operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");

        CALL_AND_HANDLE(m_mode_operators[nu].push_back(opTerms(r)), "Failed to push mode operator term to list.");
        m_mode_operators[nu].back().op() = op;
    } 
    
    void bind(sOP&& op, const element_container_type& r, size_t nu)
    {
        //first we check that non of the r-indices in this object have already been bound
        ASSERT(valid_rvals(r, nu), "Unable to bind operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");
        CALL_AND_HANDLE(m_mode_operators[nu].push_back(opTerms(r)), "Failed to push mode operator term to list.");
        m_mode_operators[nu].back().op() = std::move(op);
    } 

    size_t index(size_t r, size_t nu) const
    {
        size_t count = 0;
        for(const auto& li : m_mode_operators[nu])
        {
            if(li.contains_index(r)){return count;}
            ++count;
        }
        return count;
    }

    void clear()
    {
        m_mode_operators.clear();
        m_coeff.clear();
        m_nterms = 0;
    }

    const mode_terms_type& operators(size_t nu) const
    {
        return m_mode_operators[nu];
    }

    const mode_terms_type& operator[](size_t nu) const
    {
        return m_mode_operators[nu];
    }

    const mode_terms_type& operator()(size_t nu) const
    {
        return m_mode_operators[nu];
    }

    const element_type& operator()(size_t nu, size_t k) const
    {
        return m_mode_operators[nu][k];
    }

    const std::vector<T>& coeff() const{return m_coeff;}
    const T& coeff(size_t r)const{return m_coeff[r];}
    T& coeff(size_t r){return m_coeff[r];}

    size_t nterms() const{return m_nterms;}
    size_t nterms(size_t nu) const{return m_mode_operators[nu].size();}
    size_t nmodes() const{return m_mode_operators.size();}

    iterator begin() {  return iterator(m_mode_operators.begin());  }
    iterator end() {  return iterator(m_mode_operators.end());  }
    const_iterator begin() const {  return const_iterator(m_mode_operators.begin());  }
    const_iterator end() const {  return const_iterator(m_mode_operators.end());  }

    reverse_iterator rbegin() {  return reverse_iterator(m_mode_operators.rbegin());  }
    reverse_iterator rend() {  return reverse_iterator(m_mode_operators.rend());  }
    const_reverse_iterator rbegin() const {  return const_reverse_iterator(m_mode_operators.rbegin());  }
    const_reverse_iterator rend() const {  return const_reverse_iterator(m_mode_operators.rend());  }

protected:
    //check whether any of the new r values are already bound for this mode.
    bool valid_rvals(const element_container_type& r, size_t nu)
    {
        for(auto& combop : m_mode_operators[nu])
        {
            for(const auto& ri : combop)
            {
                if(std::find(r.begin(), r.end(), ri) != r.end()){return false;}
            }
        }
        return true;
    }

protected:
    container_type m_mode_operators;
    size_t m_nterms;
    std::vector<T> m_coeff;

#ifdef CEREAL_LIBRARY_FOUND
public:
    template <typename archive>
    void serialize(archive& ar)
    {
    }
#endif
};  //class compressedSOP


}   //namespace utils

#endif  //UTILS_OPERATOR_GEN_COMPRESSED_SOP_HPP

