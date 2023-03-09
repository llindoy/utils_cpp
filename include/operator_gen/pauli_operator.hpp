#ifndef UTILS_OPERATOR_GEN_PAULI_OPERATORS_HPP
#define UTILS_OPERATOR_GEN_PAULI_OPERATORS_HPP

#include <linalg/linalg.hpp>
#include "one_body_operator.hpp"

namespace utils{

namespace pauli
{


template <typename T> 
class sigma_p : public one_body_operator<T>
{
public:
    sigma_p() {}
    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form sigma_+ as a diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
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
            RAISE_EXCEPTION("Failed to construct sigma_+ as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
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
            RAISE_EXCEPTION("Failed to construct sigma_+ as dense.");
        }
    }
};

template <typename T> 
class sigma_m : public one_body_operator<T>
{
public:
    sigma_m() {}

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form sigma_+ as a diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
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
            RAISE_EXCEPTION("Failed to construct sigma_- as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
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
            RAISE_EXCEPTION("Failed to construct sigma_- as dense.");
        }
    }
};

template <typename T> 
class sigma_x : public one_body_operator<T>
{
public:
    sigma_x() {}

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form sigma_+ as a diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {   
                    ++nnz;
                }
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
                if(op->contains_lowered_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = 1.0;
                    colind[counter] = op->get_lowered_index(i, index);
                    ++counter;
                }
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
            RAISE_EXCEPTION("Failed to construct sigma_x as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {
                    mat(i, op->get_lowered_index(i, index)) = 1.0;
                }
                if(op->contains_raised_state(i, index))
                {
                    mat(i, op->get_raised_index(i, index)) = 1.0;
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sigma_x as dense.");
        }
    }
};

template <typename T> 
class sigma_y : public one_body_operator<T>
{
public:
    static_assert(linalg::is_complex<T>::value, "Cannot create sigma_y object in a real valued matrix.");
public:
    sigma_y() {}

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form sigma_+ as a diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {   
                    ++nnz;
                }
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
                if(op->contains_lowered_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = T(0, 1.0);
                    colind[counter] = op->get_lowered_index(i, index);
                    ++counter;
                }
                if(op->contains_raised_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = T(0.0,-1.0);
                    colind[counter] = op->get_raised_index(i, index);
                    ++counter;
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sigma_y as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {
                    mat(i, op->get_lowered_index(i, index)) = T(0, 1);
                }
                if(op->contains_raised_state(i, index))
                {
                    mat(i, op->get_raised_index(i, index)) = T(0, -1);
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sigma_y as dense.");
        }
    }
};

template <typename T> 
class sigma_z : public one_body_operator<T>
{
public:
    sigma_z() {}

    virtual bool is_sparse() const{return true;}
    virtual bool is_diagonal() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                mat(i, i) = (n == 0 ? 1.0 : -1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sigma_z as dense.");
        }
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
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
            RAISE_EXCEPTION("Failed to construct sigma_z as csr.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(op->dim(index) == 2, "Unable to create pauli matrix for mode with dimension greater than 2.");
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                mat(i, i) = (n == 0 ? 1.0 : -1.0);
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct sigma_z as dense.");
        }
    }
};

}//namespace pauli

}//namespace utils

#endif

