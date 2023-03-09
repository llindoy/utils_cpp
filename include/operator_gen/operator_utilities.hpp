#ifndef OPERATOR_UTILITIES_HPP
#define OPERATOR_UTILITIES_HPP

#include <memory>
#include <linalg/linalg.hpp>

/*  
template <typename T> 
static inline void csr_csr_product(const linalg::csr_matrix<T>& A, const linalg::csr_matrix<T>& B, linalg::csr_matrix<T>& C)
{
    auto Arp = A.rowptr();
    auto Aci= A.colind();
    auto Ab = A.buffer();

    auto Brp = B.rowptr();
    auto Bci= B.colind();
    auto Bb = B.buffer();

    //first we determine the number of non-zeros in the matrix
    size_t nnz = 0;
    for(size_t i = 0; i < A.size(0); ++i

    for(size_t i = 0; i < S.size(0); ++i)
    {
        for(size_t j = static_cast<size_t>(rowptr[i]); j < static_cast<size_t>(rowptr[i+1]); ++j)
        {
            D(i, colind[j]) = buffer[j];
        }
    }
}

static inline void csr_csr_sum
*/

template <typename T> 
static inline void csr_to_dense(const linalg::csr_matrix<T>& S, linalg::matrix<T>& D)
{
    CALL_AND_HANDLE(D.resize(S.size(0), S.size(1)), "Failed to resize dense matrix.");   D.fill_zeros();

    auto rowptr = S.rowptr();
    auto colind = S.colind();
    auto buffer = S.buffer();
    for(size_t i = 0; i < S.size(0); ++i)
    {
        for(size_t j = static_cast<size_t>(rowptr[i]); j < static_cast<size_t>(rowptr[i+1]); ++j)
        {
            D(i, colind[j]) = buffer[j];
        }
    }
}

template <typename T, typename RT>
static inline void dense_to_csr(const linalg::matrix<T>& D, linalg::csr_matrix<T>& S, RT prune_tol = -1) 
{
    CALL_AND_HANDLE(S.resize(D.size(0), D.size(1)), "Failed to resize sparse matrix.");
    auto rowptr = S.rowptr();   rowptr[0] = 0;
    size_t nnz = 0;
    for(size_t i=0; i < D.size(0); ++i)
    {
        for(size_t j=0; j < D.size(1); ++j)
        {
            if(prune_tol <= linalg::abs(D(i, j)) )
            {
                ++nnz;
            }
        }
        rowptr[i+1] = nnz;
    }
    CALL_AND_HANDLE(S.resize(nnz), "Failed to set nnz in sparse matrix.");
    
    auto buffer = S.buffer();
    auto colind = S.colind();

    size_t counter = 0;
    for(size_t i=0; i < D.size(0); ++i)
    {
        for(size_t j=0; j < D.size(1); ++j)
        {
            if(prune_tol <= linalg::abs(D(i, j)) )
            {
                buffer[counter] = D(i, j);
                colind[counter] = j;
                ++counter;
            }
        }
    }
}

#endif

