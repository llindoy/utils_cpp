#ifndef UTILS_OPERATOR_GEN_COMPRESSED_SOP_HPP
#define UTILS_OPERATOR_GEN_COMPRESSED_SOP_HPP

#include <linalg/linalg.hpp>

#include <memory>
#include <list>
#include <vector>
#include <algorithm>
#include <map>

#include <linalg/linalg.hpp>

#include "sSOP.hpp"

#ifdef CEREAL_LIBRARY_FOUND
#include <cereal/types/vector.hpp>
#endif

namespace utils
{
template <typename T> class cSSOP;
template <typename T> std::ostream& operator<<(std::ostream& os, const cSSOP<T>& op);

//a generic sum of product operator object.
template <typename T>
class cSSOP
{
public:
    using element_type = std::map<std::string, std::vector<size_t>>;
    using element_container_type = std::vector<size_t>;
    using container_type = std::vector<element_type>;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using reverse_iterator = typename container_type::reverse_iterator;
    using const_reverse_iterator = typename container_type::const_reverse_iterator;

public:

    cSSOP() : m_nterms(0){}
    cSSOP(size_t nterms, size_t nmodes) 
    try : m_mode_operators(nmodes), m_nterms(nterms), m_coeff(nterms)
    {
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to construct sum of product operator object.");
    }    

    cSSOP(const cSSOP& o) = default;
    cSSOP(cSSOP&& o) = default;

    cSSOP(const sSOP<T>& op) 
    {
        try
        {
            initialise(op);
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sum of product operator object.");
        }    
    }

    cSSOP& operator=(const cSSOP& o) = default;
    cSSOP& operator=(cSSOP&& o) = default;

    cSSOP& operator=(const sSOP<T>& o)
    {
        try
        {
            initiailise(o);
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sum of product operator object.");
        }    
    }


    void initialise(const sSOP<T>& o, size_t nmodes = 0)
    {
        try
        {
            m_label = o.label();
            resize(o.nterms(), nmodes > o.nmodes() ? nmodes : o.nmodes());
            size_t r = 0;
            for(const auto& nbo : o)
            {
                std::map<size_t, std::list<std::string>> terms;
                m_coeff[r] = nbo.coeff();
                for(const auto& op : nbo)
                {
                    terms[op.mode()].push_back(op.op());
                }

                for(const auto& kv : terms)
                {
                    std::string val;

                    for(const std::string& sv : kv.second)
                    {
                        val += sv;
                    }
                    this->insert(val, r, kv.first);
                }
                ++r;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to initialise compressed sum of product operator object.");
        }    
    }

    void resize(size_t nterms, size_t nmodes)
    {
        try
        {
            clear();
            m_coeff.resize(nterms); std::fill(m_coeff.begin(), m_coeff.end(), T(1));
            m_mode_operators.resize(nmodes);
            m_nterms = nterms;
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to resize sp hamiltonian object.");
        }
    }
   
    void insert(const sOP& op, const element_container_type& r)
    {
        size_t nu = op.mode();
        try
        {
            //first we check that non of the r-indices in this object have already been bound
            ASSERT(valid_rvals(r, nu), "Unable to insert operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");
            auto& v = m_mode_operators[nu][op.op()];
            v.insert(std::end(v), std::begin(r), std::end(r));
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to insert operator into cSSOP");
        }
    } 
    
    void insert(const std::string& op, const element_container_type& r, size_t nu)
    {
        try
        {
            //first we check that non of the r-indices in this object have already been bound
            ASSERT(valid_rvals(r, nu), "Unable to insert operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");
            auto& v = m_mode_operators[nu][op];
            v.insert(std::end(v), std::begin(r), std::end(r));
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to insert operator into cSSOP");
        }
    } 
    

    void insert(const sOP& op, size_t r)
    {
        try
        {
            size_t nu = op.mode();
            //first we check that non of the r-indices in this object have already been bound
            ASSERT(valid_rvals(r, nu), "Unable to insert operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");
            CALL_AND_HANDLE(m_mode_operators[nu][op.op()].push_back(r), "Failed to push mode operator term to list.");
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to insert operator into cSSOP");
        }
    } 
    
    void insert(const std::string& op, size_t r, size_t nu)
    {
        try
        {
            //first we check that non of the r-indices in this object have already been bound
            ASSERT(valid_rvals(r, nu), "Unable to insert operator to sum of product operator.  At least one of the r indices specified has previously been bound for this mode.");

            CALL_AND_HANDLE(m_mode_operators[nu][op].push_back(r), "Failed to push mode operator term to list.");
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to insert operator into cSSOP");
        }
    } 



    size_t index(size_t r, size_t nu) const
    {
        size_t count = 0;
        for(const auto& li : m_mode_operators[nu])
        {
            if(std::find(li.second.begin(), li.second().end(), r) != li.second.end()){return count;}
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

    const element_type& operators(size_t nu) const
    {
        return m_mode_operators[nu];
    }

    const element_type& operator[](size_t nu) const
    {
        return m_mode_operators[nu];
    }

    const element_type& operator()(size_t nu) const
    {
        return m_mode_operators[nu];
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


    friend std::ostream& operator<< <T>(std::ostream& os, const cSSOP<T>& op);

    const std::string& label() const{return m_label;}
    std::string& label(){return m_label;}
protected:
    //check whether any of the new r values are already bound for this mode.
    bool valid_rvals(const element_container_type& r, size_t nu)
    {
        //iterate over all elements in the map associated with mode nu
        for(const auto& combop : m_mode_operators[nu])
        {
            for(const auto& ri : combop.second)
            {
                if(std::find(r.begin(), r.end(), ri) != r.end()){return false;}
            }
        }
        return true;
    }

    bool valid_rvals(size_t r, size_t nu)
    {
        //iterate over all elements in the map associated with mode nu
        for(const auto& combop : m_mode_operators[nu])
        {
            const auto& vec = combop.second;
            if(std::find(vec.begin(), vec.end(), r) != vec.end()){return false;}
        }
        return true;
    }

protected:
    container_type m_mode_operators;
    size_t m_nterms;
    std::vector<T> m_coeff;
    std::string m_label;

#ifdef CEREAL_LIBRARY_FOUND
public:
    template <typename archive>
    void serialize(archive& ar)
    {
    }
#endif
};  //class cSSOP


template <typename T>
std::ostream& operator<<(std::ostream& os, const utils::cSSOP<T>& op)
{
    if(!op.label().empty()){os << op.label() << ": " << std::endl;}
    for(size_t nu = 0; nu < op.m_mode_operators.size(); ++nu)
    {
        os << "mode: " << nu << std::endl;
        for(const auto& combop : op.m_mode_operators[nu])
        {
            os << combop.first << ": ";
            for(const auto& ri : combop.second)
            {
                os << ri << " ";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    return os;
}

}   //namespace utils

#endif  //UTILS_OPERATOR_GEN_COMPRESSED_SOP_HPP

