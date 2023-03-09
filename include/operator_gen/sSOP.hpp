#ifndef UTILS_OPERATOR_GEN_LIB_SSOP_HPP
#define UTILS_OPERATOR_GEN_LIB_SSOP_HPP

#include <iostream>
#include <list>
#include <regex>
#include <io/input_wrapper.hpp>
#include <tmp_funcs.hpp>

namespace utils
{


class sOP
{
public:
    sOP() : m_fermionic(false){}
    sOP(const std::string& op, size_t mode) : m_op_data(op), m_mode(mode), m_fermionic(false){} 
    sOP(const std::string& op, size_t mode, bool fermionic) : m_op_data(op), m_mode(mode), m_fermionic(fermionic){} 

    sOP(const sOP& o) = default;
    sOP(sOP&& o) = default;

    sOP& operator=(const sOP& o) = default;
    sOP& operator=(sOP&& o) = default;

    
    const std::string& op() const{return m_op_data;}
    std::string& op(){return m_op_data;}

    const size_t& mode() const{return m_mode;}
    size_t& mode(){return m_mode;}

    const bool& fermionic() const{return m_fermionic;}
    bool& fermionic(){return m_fermionic;}

    friend std::ostream& operator<<(std::ostream& os, const sOP& op);
    friend bool operator==(const sOP& A, const sOP& B);
    friend bool operator!=(const sOP& A, const sOP& B);
protected:
    std::string m_op_data;
    size_t m_mode;
    bool m_fermionic;
};

std::ostream& operator<<(std::ostream& os, const utils::sOP& op)
{
    if(op.m_fermionic)
    {
        os << "fermi_" << op.m_op_data << "_" << op.m_mode;
    }
    else
    {
        os << op.m_op_data << "_" << op.m_mode; 
    }
    return os;
}

bool operator==(const sOP& A, const sOP& B)
{
    return (A.m_op_data == B.m_op_data && A.m_mode == B.m_mode && A.m_fermionic == B.m_fermionic);
}

bool operator!=(const sOP& A, const sOP& B){return !(A == B);}

static inline sOP fermion_operator(const std::string& op, size_t mode)
{
    std::string label(op);
    io::remove_whitespace_and_to_lower(label);

    return sOP(op, mode, true);
};


}


namespace utils
{

//product of string operators
class sPOP
{
public:
    using iterator = typename std::list<sOP>::iterator;
    using const_iterator = typename std::list<sOP>::const_iterator;
    using reverse_iterator = typename std::list<sOP>::reverse_iterator;
    using const_reverse_iterator = typename std::list<sOP>::const_reverse_iterator;

public:
    sPOP(){}


    sPOP(const sOP& ops) {m_ops.push_back(ops);}
    sPOP(sOP&& ops)  {m_ops.push_back(std::forward<sOP>(ops));}

    template <typename ... Args>
    sPOP(const sOP& ops, Args&& ... args) 
    {
        m_ops.push_back(ops);
        unpack_args(std::forward<Args>(args)...);
    }

    template <typename ... Args>
    sPOP(sOP&& ops, Args&& ... args) 
    {
        m_ops.push_back(std::forward<sOP>(ops));
        unpack_args(std::forward<Args>(args)...);
    }

    sPOP( const std::list<sOP>& ops) : m_ops(ops){}
    sPOP( std::list<sOP>&& ops) : m_ops(std::forward<std::list<sOP>>(ops)){}

    sPOP(const sPOP& o) = default;
    sPOP(sPOP&& o) = default;

    sPOP& operator=(const sPOP& o) = default;
    sPOP& operator=(sPOP&& o) = default;

    const std::list<sOP>& ops() const{return m_ops;}
    std::list<sOP>& ops(){return m_ops;}
    void clear(){m_ops.clear();}
    void append(const sOP& o){m_ops.push_back(o);}
    void prepend(const sOP& o){m_ops.push_front(o);}

    sPOP& operator*=(const sOP& b)
    {
        m_ops.push_back(b);
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& os, const sPOP& op);

    friend bool operator==(const sPOP& A, const sPOP& B);
    friend bool operator!=(const sPOP& A, const sPOP& B);
protected: 
    template <typename ... Args>
    void unpack_args(const sOP& ops, Args&& ... args)
    {
        m_ops.push_back(ops);
        unpack_args(std::forward<Args>(args)...);
    }

    template <typename ... Args>
    void unpack_args(sOP&& ops, Args&& ... args)
    {
        m_ops.push_back(std::forward<sOP>(ops));
        unpack_args(std::forward<Args>(args)...);
    }

    void unpack_args(){}


public:
    iterator begin() {  return iterator(m_ops.begin());  }
    iterator end() {  return iterator(m_ops.end());  }
    const_iterator begin() const {  return const_iterator(m_ops.begin());  }
    const_iterator end() const {  return const_iterator(m_ops.end());  }

    reverse_iterator rbegin() {  return reverse_iterator(m_ops.rbegin());  }
    reverse_iterator rend() {  return reverse_iterator(m_ops.rend());  }
    const_reverse_iterator rbegin() const {  return const_reverse_iterator(m_ops.rbegin());  }
    const_reverse_iterator rend() const {  return const_reverse_iterator(m_ops.rend());  }
protected:
    std::list<sOP> m_ops;
};




std::ostream& operator<<(std::ostream& os, const utils::sPOP& op)
{
    const auto separator = " ";    const auto* sep = "";
    for(const auto& t : op.m_ops)
    {
        os << sep << t;
        sep = separator;
    }
    return os;
}

bool operator==(const sPOP& A, const sPOP& B)
{
    if(A.m_ops.size() == B.m_ops.size())
    {
        for(auto z : zip(A.m_ops, B.m_ops))
        {
            if(std::get<0>(z) != std::get<1>(z)){return false;}
        }
        return true;
    }
    return false;
}

bool operator!=(const sPOP& A, const sPOP& B){return !(A == B);}

}
  
utils::sPOP operator*(const utils::sOP& a, const utils::sOP& b)
{
    utils::sPOP ret{{a, b}};
    return ret;
}

utils::sPOP operator*(const utils::sPOP& a, const utils::sOP& b)
{
    utils::sPOP ret(a);
    ret.append(b);
    return ret;
}

utils::sPOP operator*(const utils::sOP& a, const utils::sPOP& b)
{
    utils::sPOP ret(b);
    ret.prepend(a);
    return ret;
}


utils::sPOP operator*(const utils::sPOP& a, const utils::sPOP& b)
{
    utils::sPOP ret(a);
    for(const auto& op : b.ops())
    {
        ret.append(op);
    }
    return ret;
}




namespace utils
{

template <typename T> class sNBO;
template <typename T> std::ostream& operator<<(std::ostream& os, const sNBO<T>& op);

template <typename T> 
class sNBO
{
public:
    using iterator = typename sPOP::iterator;
    using const_iterator = typename sPOP::const_iterator;
    using reverse_iterator = typename sPOP::reverse_iterator;
    using const_reverse_iterator = typename sPOP::const_reverse_iterator;
public:
    sNBO(){}

    sNBO(const sPOP& p) : m_coeff(T(1.0)), m_ops(p) {}
    sNBO(sPOP&& p) : m_coeff(T(1.0)), m_ops(std::forward<sPOP>(p)) {}

    template <typename ... Args>
    sNBO(const sOP& p, Args&&... args) : m_coeff(T(1.0)), m_ops(p, std::forward<Args>(args)...) {}
    template <typename ... Args>
    sNBO(sOP&& p, Args&&... args) : m_coeff(T(1.0)), m_ops(std::forward<sOP>(p), std::forward<Args>(args)...) {}

    sNBO(const T& coeff, const sPOP& p) : m_coeff(coeff), m_ops(p) {}
    sNBO(const T& coeff, sPOP&& p) : m_coeff(coeff), m_ops(std::forward<sPOP>(p)) {}

    template <typename ... Args>
    sNBO(const T& coeff, const sOP& p, Args&&... args) : m_coeff(coeff), m_ops(p, std::forward<Args>(args)...) {}
    template <typename ... Args>
    sNBO(const T& coeff, sOP&& p, Args&&... args) : m_coeff(coeff), m_ops(std::forward<sOP>(p), std::forward<Args>(args)...) {}

    sNBO(const sNBO& o) = default;
    sNBO(sNBO&& o) = default;

    template <typename U>
    sNBO(const sNBO<U>& o) : m_coeff(o.coeff()), m_ops(o.ops()){}

    sNBO& operator=(const sNBO& o) = default;
    sNBO& operator=(sNBO&& o) = default;

    template <typename U>
    sNBO& operator=(const sNBO<U>& o)
    {
        m_coeff = o.coeff();
        m_ops = o.ops();
    }

    const T& coeff() const{return m_coeff;}
    T& coeff(){return m_coeff;}

    const std::list<sOP>& ops() const{return m_ops.ops();}
    std::list<sOP>& ops(){return m_ops.ops();}

    const sPOP& pop() const{return m_ops;}
    sPOP& pop(){return m_ops;}

    void clear(){m_coeff = T(1);    m_ops.clear();}
    void append(const sOP& o){m_ops.append(o);}

    sNBO<T>& operator*=(const sOP& b)
    {
        m_ops.append(b);
        return *this;
    }

    sNBO<T>& operator*=(const sPOP& b)
    {
        for(const auto& op : b.ops())
        {
            m_ops.append(op);
        }
        return *this;
    }

    template <typename U>
    sNBO<T>& operator*=(const sNBO<U>& b)
    {
        m_coeff *= b.coeff();
        for(const auto& op : b.ops())
        {
            m_ops.append(op);
        }
        return *this;
    }

    sNBO<T>& operator*=(const T& b)
    {
        m_coeff*=b;
        return *this;
    }

    friend std::ostream& operator<< <T>(std::ostream& os, const sNBO<T>& op);


public:
    iterator begin() {  return iterator(m_ops.begin());  }
    iterator end() {  return iterator(m_ops.end());  }
    const_iterator begin() const {  return const_iterator(m_ops.begin());  }
    const_iterator end() const {  return const_iterator(m_ops.end());  }

    reverse_iterator rbegin() {  return reverse_iterator(m_ops.rbegin());  }
    reverse_iterator rend() {  return reverse_iterator(m_ops.rend());  }
    const_reverse_iterator rbegin() const {  return const_reverse_iterator(m_ops.rbegin());  }
    const_reverse_iterator rend() const {  return const_reverse_iterator(m_ops.rend());  }

protected:
    T m_coeff;
    sPOP m_ops;
};

template <typename T> 
std::ostream& operator<<(std::ostream& os, const utils::sNBO<T>& op)
{
    os << op.coeff() << op.pop();
    return os;
}

}   //namespace utils

template <typename T>
utils::sNBO<T> operator*(const T& a, const utils::sOP& b)
{
    utils::sNBO<T> ret(a, {b});
    return ret;
}


template <typename T>
utils::sNBO<T> operator*(const utils::sOP& b, const T& a)
{
    utils::sNBO<T> ret(a, {b});
    return ret;
}


template <typename T>
utils::sNBO<T> operator*(const utils::sNBO<T>& a, const T& b)
{
    utils::sNBO<T> ret(a);
    ret.m_coeff*=b;
    return ret;
}

template <typename T>
utils::sNBO<T> operator*(const T& b, const utils::sNBO<T>& a)
{
    utils::sNBO<T> ret(a);
    ret.m_coeff*=b;
    return ret;
}



template <typename T>
utils::sNBO<T> operator*(const utils::sNBO<T>& a, const utils::sOP& b)
{
    utils::sNBO<T> ret(a);
    ret.pop().append(b);
    return ret;
}


template <typename T>
utils::sNBO<T> operator*(const utils::sOP& b, const utils::sNBO<T>& a)
{
    utils::sNBO<T> ret(a);
    ret.pop().prepend(b);
    return ret;
}

template <typename T>
utils::sNBO<T> operator*(const utils::sNBO<T>& a, const utils::sNBO<T>& b)
{
    utils::sNBO<T> ret(a);
    ret.m_coeff *= b.m_coeff;
    for(const auto& op : b.ops())
    {
        ret.pop().append(op);
    }
    return ret;
}

template <typename T>
utils::sNBO<T> operator*(const utils::sNBO<T>& a, const utils::sPOP& b)
{
    utils::sNBO<T> ret(a);
    for(const auto& op : b.ops())
    {
        ret.pop().append(op);
    }
    return ret;
}

template <typename T>
utils::sNBO<T> operator*(const utils::sPOP& a, const utils::sNBO<T>& b)
{
    utils::sNBO<T> ret(a);
    ret.coeff() = b.coeff();
    for(const auto& op : b.ops())
    {
        ret.pop().append(op);
    }
    return ret;
}

namespace utils
{

template <typename T> class sSOP;
template <typename T> std::ostream& operator<<(std::ostream& os, const sSOP<T>& op);

//the string sum of product operator class used for storing the representation of the Hamiltonian of interest.
template <typename T> 
class sSOP
{
public:
    using iterator = typename std::list<sNBO<T>>::iterator;
    using const_iterator = typename std::list<sNBO<T>>::const_iterator;
    using reverse_iterator = typename std::list<sNBO<T>>::reverse_iterator;
    using const_reverse_iterator = typename std::list<sNBO<T>>::const_reverse_iterator;

public:
    sSOP(){}
    sSOP(const sNBO<T>& str){m_string.push_back(str);}
    sSOP(const std::list<sNBO<T>>& str) : m_string(str){}
    template <typename ... Args>
    sSOP(const sNBO<T>& str){m_string.push_back(str);}

    sSOP(sNBO<T>&& str){m_string.push_back(std::forward(str));}
    sSOP(std::list<sNBO<T>>&& str) : m_string(std::forward(str)){}

    sSOP(const sSOP& o) = default;
    sSOP(sSOP&& o) = default;

    sSOP& operator=(const sSOP& o) = default;
    sSOP& operator=(sSOP&& o) = default;

    sSOP<T>& operator+=(const sSOP<T>& a)
    {
        for(auto& t : a.m_string)
        {
            m_string.push_back(t);
        }
        return *this;
    }

    sSOP<T>& operator+=(const sNBO<T>& a)
    {
        m_string.push_back(a);
        return *this;
    }

    sSOP<T>& operator*=(const T& a)
    {
        for(auto& op : m_string)
        {
            op *= a;
        }
        return *this;
    }

    friend std::ostream& operator<< <T>(std::ostream& os, const sSOP<T>& op);

    size_t nterms() const{return m_string.size();}
public:
    iterator begin() {  return iterator(m_string.begin());  }
    iterator end() {  return iterator(m_string.end());  }
    const_iterator begin() const {  return const_iterator(m_string.begin());  }
    const_iterator end() const {  return const_iterator(m_string.end());  }

    reverse_iterator rbegin() {  return reverse_iterator(m_string.rbegin());  }
    reverse_iterator rend() {  return reverse_iterator(m_string.rend());  }
    const_reverse_iterator rbegin() const {  return const_reverse_iterator(m_string.rbegin());  }
    const_reverse_iterator rend() const {  return const_reverse_iterator(m_string.rend());  }
protected:
    std::list<sNBO<T>> m_string;
};


template <typename T> 
std::ostream& operator<<(std::ostream& os, const utils::sSOP<T>& op)
{
    const auto separator = "";    const auto* sep = "";
    const auto plus = "+";
    for(const auto& t : op)
    {
        sep = std::real(t.coeff()) > 0 ? plus : separator;
        os << sep << t;
    }
    return os;
}
}




template <typename T>
utils::sSOP<T> operator+(const utils::sSOP<T>& a, const utils::sSOP<T>& b)
{
    utils::sSOP<T> ret(a);
    for(auto& t : b.m_string)
    {
        ret.m_string.push_back(t);
    }
    return  ret;
}

template <typename T>
utils::sSOP<T> operator+(const utils::sNBO<T>& a, const utils::sSOP<T>& b)
{
    utils::sSOP<T> ret(b);
    ret.m_string.push_back(a);
    return  ret;
}

template <typename T>
utils::sSOP<T> operator+(const utils::sSOP<T>& b, const utils::sNBO<T>& a)
{
    utils::sSOP<T> ret(b);
    ret.m_string.push_back(a);
    return  ret;
}

template <typename T>
utils::sSOP<T> operator+(const utils::sOP& a, const utils::sSOP<T>& b)
{
    utils::sSOP<T> ret(b);
    ret.m_string.push_back({1.0, a});
    return  ret;
}

template <typename T>
utils::sSOP<T> operator+(const utils::sSOP<T>& b, const utils::sOP& a)
{
    utils::sSOP<T> ret(b);
    ret.m_string.push_back({1.0, a});
    return  ret;
}



#endif

