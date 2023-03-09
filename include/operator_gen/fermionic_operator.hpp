#ifndef UTILS_OPERATOR_GEN_FERMION_OPERATORS_HPP
#define UTILS_OPERATOR_GEN_FERMION_OPERATORS_HPP

#include <linalg/linalg.hpp>
#include "one_body_operator.hpp"

namespace utils
{

namespace fermion
{

//form the one body fermion creation operator.  This operator does not correctly satisfy the fermion anticommutation operators  
template <typename T> 
class creation : public one_body_operator<T>
{
public:
    creation() {}
    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form fermion creation operator as a diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_raised_state(i, index))
                {   
                    ++nnz;
                }
            }

            mat.resize(nnz, op->nstates(), op->nstates());
            auto rowptr = mat.rowptr();    rowptr[0] = 0;
            auto colind = mat.colind();
            auto buffer = mat.buffer();

            size_t counter = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_raised_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = 1.0;
                    colind[counter] = op->get_raised_index(i, index);
                    ++counter;
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion creation operator as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
        ASSERT(index < op->nmodes(), "Index out of bounds.");
        mat.resize(op->nstates(), op->nstates());
        mat.fill_zeros();
    
        try
        {
        for(size_t i = 0; i < op->nstates(); ++i)
        {
            if(op->contains_raised_state(i, index))
            {
                mat(i, op->get_raised_index(i, index)) = 1.0;
            }
        }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion creation operator as dense.");
        }
    }
};

template <typename T> 
class annihilation : public one_body_operator<T>
{
public:
    annihilation() {}

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form fermion annihilation operator as a diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {   
                    ++nnz;
                }
            }

            mat.resize(nnz, op->nstates(), op->nstates());
            auto rowptr = mat.rowptr();    rowptr[0] = 0;
            auto colind = mat.colind();
            auto buffer = mat.buffer();

            size_t counter = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = 1.0;
                    colind[counter] = op->get_lowered_index(i, index);
                    ++counter;
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion annihilation operator as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {
                    mat(i, op->get_lowered_index(i, index)) = 1.0;
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion annihilation operator as dense.");
        }
    }
};

template <typename T> 
class number : public one_body_operator<T>
{
public:
    number() {}

    virtual bool is_sparse() const{return true;}
    virtual bool is_diagonal() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                mat(i, i) = (n == 0 ? 0.0 : 1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion number operator as dense.");
        }
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
    
            mat.resize(op->nstates(), op->nstates(), op->nstates());
            auto rowptr = mat.rowptr();    rowptr[0] = 0;
            auto colind = mat.colind();
            auto buffer = mat.buffer();

            size_t counter = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                buffer[counter] = (n == 0 ? 0.0 : 1.0);
                colind[counter] = i;
                rowptr[i+1] = counter;
            }
            
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion number operator as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                mat(i, i) = (n == 0 ? 0.0 : 1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct fermion number operator as dense.");
        }
    }
};


template <typename T> 
class jordan_wigner: public one_body_operator<T>
{
public:
    jordan_wigner() {}

    virtual bool is_sparse() const{return true;}
    virtual bool is_diagonal() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                std::cerr << i << " " << n << " " << index << std::endl;
                mat(i, i) = (n == 0 ? 1.0 : -1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct jordan_wigner as dense.");
        }
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
    
            mat.resize(op->nstates(), op->nstates(), op->nstates());
            auto rowptr = mat.rowptr();    rowptr[0] = 0;
            auto colind = mat.colind();
            auto buffer = mat.buffer();

            size_t counter = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                buffer[counter] = (n == 0 ? 1.0 : -1.0);
                colind[counter] = i;
                rowptr[i+1] = counter;
            }
            
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct jordan_wigner as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create fermion operator matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                std::cerr << i << " " << n << " " << index << std::endl;
                mat(i, i) = (n == 0 ? 1.0 : -1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct jordan_wigner as dense.");
        }
    }
};
}  //namespace fermion
}   //namespace utils

#endif

