#ifndef UTILS_OPERATOR_GEN_BOSONIC_OPERATORS_HPP
#define UTILS_OPERATOR_GEN_BOSONIC_OPERATORS_HPP

#include "one_body_operator.hpp"

namespace utils
{

namespace bosonic
{

template <typename T> 
class one_body_operator
{
public:
    one_body_operator() {}
    virtual ~one_body_operator(){}

    virtual bool is_diagonal() const{return false;}
    virtual bool is_sparse() const{return false;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const = 0;
    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const = 0;
    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const = 0;
};

/*
 * Class for handling creation operators and displaced and squeezed creation operators
 */
template <typename T> 
class creation : public one_body_operator<T>
{
public:
    creation() {}
    
    creation& displace(const T& disp)
    {
        ASSERT(!m_has_disp, "Cannot handle double displacement.");
        m_disp = disp;
        m_has_disp = true;
        return *this;
    }

    creation& squeeze(const T& squeeze)
    {
        ASSERT(!m_has_squeeze, "Cannot handle double displacement.");
        if(m_has_squeeze && !m_has_disp){m_squeeze_first = true;}
        m_squeeze = squeeze;
        m_has_squeeze = true;
        return *this;
    }

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form creation operator as diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            using RT = typename linalg::get_real_type<T>::type;
            ASSERT(index < op->nmodes(), "Index out of bounds.");
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {   
                    ++nnz;
                }
                if(m_has_disp)
                {
                    ++nnz;
                }       
                if(m_has_squeeze)
                {
                    if(op->contains_raised_state(i, index))
                    {   
                        ++nnz;
                    }
                }
            }

            mat.resize(nnz, op->nstates(), op->nstates());
            auto rowptr = mat.rowptr();    rowptr[0] = 0;
            auto colind = mat.colind();
            auto buffer = mat.buffer();

            size_t counter = 0;


            RT coeff_a = 1.0;
            T coeff_b = 0.0;
            if(m_has_squeeze)
            {
                RT r = linalg::abs(m_squeeze);

                coeff_a = std::cosh(r);
                coeff_b = std::sinh(r);
                if(linalg::is_complex<T>::value)
                {
                    RT theta = linalg::arg(m_squeeze);
                    coeff_b *= std::exp(T(0, -1)*theta);
                }
            }

            for(size_t i = 0; i < op->nstates(); ++i)
            {
                //add in the creation operator contribution.  This is scaled by coeff_a which depends on the squeeze operator
                if(op->contains_lowered_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = coeff_a*std::sqrt((1.0*n));
                    colind[counter] = op->get_lowered_index(i, index);
                    ++counter;
                }

                //next add on the diagonal term associated with the displacement.  The form of this changes if we have applied the 
                //squeeze operator before the displacement operator.
                if(m_has_disp)
                {
                    if(m_has_squeeze && m_squeeze_first)
                    {
                        buffer[counter] = (linalg::conj(m_disp)*coeff_a-m_disp*coeff_b);
                    }   
                    else
                    {
                        buffer[counter] = linalg::conj(m_disp);
                    }
                    colind[counter] = i;
                    ++counter;
                }

                //finally if we have applied a squeeze operator we also need to add the contribution from the annihilation operator.
                if(m_has_squeeze)
                {
                    if(op->contains_raised_state(i, index))
                    {
                        size_t n = op->get_occupation(i, index);
                        buffer[counter] = -coeff_b*std::sqrt((n+1.0));
                        colind[counter] = op->get_raised_index(i, index);
                        ++counter;
                    }
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();

            using RT = typename linalg::get_real_type<T>::type;

            RT coeff_a = 1.0;
            T coeff_b = 0.0;
            if(m_has_squeeze)
            {
                RT r = linalg::abs(m_squeeze);

                coeff_a = std::cosh(r);
                coeff_b = std::sinh(r);
                if(linalg::is_complex<T>::value)
                {
                    RT theta = linalg::arg(m_squeeze);
                    coeff_b *= std::exp(T(0, -1)*theta);
                }
            }
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                if(op->contains_lowered_state(i, index))
                {
                    mat(i, op->get_lowered_index(i, index)) = coeff_a*std::sqrt((1.0*n));
                }

                if(m_has_disp)
                {
                    if(m_has_squeeze && m_squeeze_first)
                    {
                        mat(i, i) = (linalg::conj(m_disp)*coeff_a-m_disp*coeff_b);
                    }   
                    else
                    {
                        mat(i, i) = linalg::conj(m_disp);
                    }       
                }

                if(m_has_squeeze)
                {
                    if(op->contains_raised_state(i, index))
                    {
                        mat(i, op->get_raised_index(i, index)) = -coeff_b*std::sqrt((n+1.0));
                    }
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

protected:
    bool m_has_disp = false;
    T m_disp = T(0);

    bool m_has_squeeze = false;
    T m_squeeze = T(0);

    bool m_squeeze_first = false;
};

/*
 * Class for handling annihilation operators and displaced and squeezed annihilation operators
 */
template <typename T> 
class annihilation : public one_body_operator<T>
{
public:
    annihilation() {}
    
    annihilation& displace(const T& disp)
    {
        ASSERT(!m_has_disp, "Cannot handle double displacement.");
        m_disp = disp;
        m_has_disp = true;
        return *this;
    }

    annihilation& squeeze(const T& squeeze)
    {
        ASSERT(!m_has_squeeze, "Cannot handle double displacement.");
        if(m_has_squeeze && !m_has_disp){m_squeeze_first = true;}
        m_squeeze = squeeze;
        m_has_squeeze = true;
        return *this;
    }

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form creation operator as diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }


    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            using RT = typename linalg::get_real_type<T>::type;
            ASSERT(index < op->nmodes(), "Index out of bounds.");
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(m_has_squeeze)
                {
                    if(op->contains_lowered_state(i, index))
                    {   
                        ++nnz;
                    }
                }
                if(m_has_disp)
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


            RT coeff_a = 1.0;
            T coeff_b = 0.0;
            if(m_has_squeeze)
            {
                RT r = linalg::abs(m_squeeze);

                coeff_a = std::cosh(r);
                coeff_b = std::sinh(r);
                if(linalg::is_complex<T>::value)
                {
                    RT theta = linalg::arg(m_squeeze);
                    coeff_b *= std::exp(T(0, 1)*theta);
                }
            }

            for(size_t i = 0; i < op->nstates(); ++i)
            {
                //add in the annihilation operator contribution.  This is scaled by coeff_a which depends on the squeeze operator
                if(m_has_squeeze)
                {
                    if(op->contains_lowered_state(i, index))
                    {
                        size_t n = op->get_occupation(i, index);
                        buffer[counter] = -coeff_b*std::sqrt((1.0*n));
                        colind[counter] = op->get_lowered_index(i, index);
                        ++counter;
                    }
                }

                //next add on the diagonal term associated with the displacement.  The form of this changes if we have applied the 
                //squeeze operator before the displacement operator.
                if(m_has_disp)
                {
                    if(m_has_squeeze && m_squeeze_first)
                    {
                        buffer[counter] = (m_disp*coeff_a-linalg::conj(m_disp)*coeff_b);
                    }   
                    else
                    {
                        buffer[counter] = m_disp;
                    }
                    colind[counter] = i;
                    ++counter;
                }

                //finally if we have applied a squeeze operator we also need to add the contribution from the annihilation operator.
                if(op->contains_raised_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = coeff_a*std::sqrt((n+1.0));
                    colind[counter] = op->get_raised_index(i, index);
                    ++counter;
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();

            using RT = typename linalg::get_real_type<T>::type;

            RT coeff_a = 1.0;
            T coeff_b = 0.0;
            if(m_has_squeeze)
            {
                RT r = linalg::abs(m_squeeze);

                coeff_a = std::cosh(r);
                coeff_b = std::sinh(r);
                if(linalg::is_complex<T>::value)
                {
                    RT theta = linalg::arg(m_squeeze);
                    coeff_b *= std::exp(T(0, -1)*theta);
                }
            }
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                if(m_has_squeeze)
                {
                    if(op->contains_lowered_state(i, index))
                    {
                        mat(i, op->get_lowered_index(i, index)) = -coeff_b*std::sqrt((1.0*n));
                    }
                }

                if(m_has_disp)
                {
                    if(m_has_squeeze && m_squeeze_first)
                    {
                        mat(i, i) = (m_disp*coeff_a-linalg::conj(m_disp)*coeff_b);

                    }   
                    else
                    {
                        mat(i, i) = m_disp;
                    }       
                }

                if(op->contains_raised_state(i, index))
                {
                    mat(i, op->get_raised_index(i, index)) = coeff_a*std::sqrt((n+1.0));
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

protected:
    bool m_has_disp = false;
    T m_disp = T(0);

    bool m_has_squeeze = false;
    T m_squeeze = T(0);

    bool m_squeeze_first = false;
};


template <typename T> 
class number : public one_body_operator<T>
{
public:
    number() {}

    virtual bool is_diagonal() const{return true;}
    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const 
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                mat(i, i) = n;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
    
            mat.resize(op->nstates(), op->nstates(), op->nstates());
            auto rowptr = mat.rowptr();    rowptr[0] = 0;
            auto colind = mat.colind();
            auto buffer = mat.buffer();

            size_t counter = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                buffer[counter] = n;
                colind[counter] = i;
                ++counter;
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                mat(i, i) = n;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }
};

template <typename T> 
class displacement : public one_body_operator<T>
{
protected:
    T m_alpha = T(0);
public:
    displacement(const T& alpha) : m_alpha(alpha) {}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Diagonal displacement operator is invalid.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Sparse displacement operator is currently not supported.");
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            //if we are only treating a single mode here then we just form the dense matrix representation of the displacement operator.
            if(op->nmodes() == 1)
            {
                form_single_mode_displacement_operator(mat, m_alpha, op->nstates());
            }
            else
            {
                size_t maxdim = op->dim(index);
                linalg::matrix<T> D;
                form_single_mode_displacement_operator(D, m_alpha, maxdim);
        
                //now that we have formed the single mode representation of the displacement operator we can go through and attempt to construct
                //its form in the full occupation_number_basis object space.

                std::vector<size_t> state(op->nmodes());
                for(size_t i = 0; i < op->nstates(); ++i)
                {
                    op->get_state(i, state);
                    size_t si = state[index];
                    for(size_t j=0; j < maxdim; ++j)
                    {
                        state[index] = j;
                        if(op->contains_state(state))
                        {
                            mat(i, op->get_index(state)) = D(si, j);
                        }
                    }
                }
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }


public:
    static void form_single_mode_displacement_operator(linalg::matrix<T>& Dk, const T& a, size_t ni)
    {   
        using real_type = typename linalg::get_real_type<T>::type;
        //form the dense displacement operator associated with a single mode.  
        T alpha = a;
        T nalpha_conj = -conj(alpha);
        real_type abs_alpha = abs(alpha);
        real_type a2 = abs_alpha*abs_alpha;
        real_type expa2 = exp(-a2/2.0);

        Dk.resize(ni, ni);
        Dk(0, 0) = expa2;
        if(ni > 1)
        {
            Dk(1, 1) = expa2;

            for(size_t n=2; n < ni; ++n)
            {
                Dk(n, 1) = alpha*Dk(n-1, 1)/sqrt(static_cast<real_type>(n));
            }
            for(size_t m=2; m < ni; ++m)
            {
                Dk(1, m) = nalpha_conj*Dk(1, m-1)/sqrt(static_cast<real_type>(m));
            }

            Dk(1, 1) = expa2*(1.0-a2);
            for(size_t n=2; n<ni; ++n){Dk(n, 1) *= (n-a2);}
            for(size_t m=2; m<ni; ++m){Dk(1, m) *= (m-a2);}

            //Now we populate the diagonals
            Dk(0, 0) = expa2; Dk(1,1) = expa2*(1-a2);
            for(size_t i=2; i < ni; ++i)
            {
                Dk(i, i) = ((2.0*i-1.0-a2) * Dk(i-1, i-1) - (i-1.0)*Dk(i-2, i-2))/static_cast<real_type>(i);
            }
            
            //now populate the first column (all values with m=0)
            for(size_t n=1; n < ni; ++n)
            {
                Dk(n, 0) = alpha/sqrt(static_cast<real_type>(n))*Dk(n-1, 0);
            }

            //and first row (all values with n=0)
            for(size_t m=1; m < ni; ++m)
            {
                Dk(0, m) = nalpha_conj/sqrt(static_cast<real_type>(m))*Dk(0, m-1);
            }

            //now we can compute all values with n > m
            for(size_t d=1; d<ni; ++d)
            {
                for(size_t n=d+2; n < ni; ++n)
                {   
                    size_t m= n-d;
                    Dk(n, m) = (m+n-1.0-a2)/sqrt(static_cast<real_type>(m*n))*Dk(n-1, m-1) - sqrt((m-1.0)*(n-1.0)/(m*n))*Dk(n-2, m-2);
                }
            }
            for(size_t d=1; d<ni; ++d)
            {
                for(size_t m=d+2; m < ni; ++m)
                {   
                    size_t n = m-d;
                    Dk(n, m) = (m+n-1.0-a2)/sqrt(static_cast<real_type>(m*n))*Dk(n-1, m-1) - sqrt((m-1.0)*(n-1.0)/(m*n))*Dk(n-2, m-2);
                }
            }
        }
    }
};

}   //namespace bosonic


/*
 * Class for handling the position operator for a boson in nondimensional units
 */
template <typename T> 
class position : public one_body_operator<T>
{
public:
    using real_type = typename linalg::get_real_type<T>::value;
public:
    position() {}
    position& displace(const T& disp)
    {
        ASSERT(!m_has_disp, "Cannot handle double displacement.");
        m_disp = disp;
        m_has_disp = true;
        return *this;
    }

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form position operator as diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            using RT = typename linalg::get_real_type<T>::type;
            ASSERT(index < op->nmodes(), "Index out of bounds.");
    
            size_t nnz = 0;
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {   
                    ++nnz;
                }
                if(m_has_disp)
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


            RT coeff_a = 1.0/std::sqrt(2.0);

            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = coeff_a*std::sqrt((1.0*n));
                    colind[counter] = op->get_lowered_index(i, index);
                    ++counter;
                }

                if(m_has_disp)
                {
                    buffer[counter] = linalg::conj(m_disp);
                    colind[counter] = i;
                    ++counter;
                }

                if(op->contains_raised_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = coeff_a*std::sqrt((n+1.0));
                    colind[counter] = op->get_raised_index(i, index);
                    ++counter;
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();

            using RT = typename linalg::get_real_type<T>::type;

            RT coeff_a = 1.0/std::sqrt(2.0);
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                if(op->contains_lowered_state(i, index))
                {
                    mat(i, op->get_lowered_index(i, index)) = coeff_a*std::sqrt((1.0*n));
                }

                if(m_has_disp)
                {
                    mat(i, i) = linalg::conj(m_disp);
                }

                if(op->contains_raised_state(i, index))
                {
                    mat(i, op->get_raised_index(i, index)) = coeff_a*std::sqrt((n+1.0));
                }
            }

        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

protected:
    bool m_has_disp = false;
    T m_disp = T(0);
};

/*
 * Class for handling the position operator for a boson in momentum units
 */
template <typename T> 
class momentum : public one_body_operator<T>
{
public:
    using real_type = typename linalg::get_real_type<T>::value;
public:
    momentum() {}

    virtual bool is_sparse() const{return true;}

    virtual void as_diagonal(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::diagonal_matrix<T>& mat) const
    {
        RAISE_EXCEPTION("Cannot form momentum operator as diagonal operator.  It contains off-diagonal terms in the occupation number basis.");
    }

    virtual void as_csr(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::csr_matrix<T>& mat) const
    {
        try
        {
            using RT = typename linalg::get_real_type<T>::type;
            ASSERT(index < op->nmodes(), "Index out of bounds.");
    
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


            RT coeff_a = 1.0/std::sqrt(2.0);

            for(size_t i = 0; i < op->nstates(); ++i)
            {
                if(op->contains_lowered_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = coeff_a*std::sqrt((1.0*n));
                    colind[counter] = op->get_lowered_index(i, index);
                    ++counter;
                }

                if(op->contains_raised_state(i, index))
                {
                    size_t n = op->get_occupation(i, index);
                    buffer[counter] = coeff_a*std::sqrt((n+1.0));
                    colind[counter] = op->get_raised_index(i, index);
                    ++counter;
                }
                rowptr[i+1] = counter;
            }
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

    virtual void as_dense(const std::shared_ptr<occupation_number_basis>& op, size_t index, linalg::matrix<T>& mat) const
    {
        try
        {
            ASSERT(index < op->nmodes(), "Index out of bounds.");
            mat.resize(op->nstates(), op->nstates());
            mat.fill_zeros();

            using RT = typename linalg::get_real_type<T>::type;

            RT coeff_a = 1.0/std::sqrt(2.0);
    
            for(size_t i = 0; i < op->nstates(); ++i)
            {
                size_t n = op->get_occupation(i, index);
                if(op->contains_lowered_state(i, index))
                {
                    mat(i, op->get_lowered_index(i, index)) = coeff_a*std::sqrt((1.0*n));
                }

                if(op->contains_raised_state(i, index))
                {
                    mat(i, op->get_raised_index(i, index)) = coeff_a*std::sqrt((n+1.0));
                }
            }

        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to construct bosonic operator.");
        }
    }

};

}   //utils

#endif

