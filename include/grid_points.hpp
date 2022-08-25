#ifndef UTILS_GRID_POINTS_HPP
#define UTILS_GRID_POINTS_HPP

#include <linalg/linalg.hpp>

namespace utils
{

template <typename T, typename ret_type = linalg::vector<T>>
ret_type linspace(T start, T stop, size_t N, bool endpoint = true)
{
    ASSERT(N > 1, "Cannot compute linspace with one point.");
    ret_type ret(N);
    T dx;
    if(endpoint)
    {
        dx = (stop-start)/(N-1);
    }
    else
    {
        dx = (stop-start)/N;
    }

    for(size_t i = 0; i < N; ++i)
    {
        ret[i] = start + dx*T(i);
    }
    return ret;
}


template <typename T, typename ret_type = linalg::vector<T>>
ret_type logspace(T start, T stop, size_t N, bool endpoint = true, T base = T(10.0))
{
    ASSERT(N > 1, "Cannot compute logspace with one point.");
    ASSERT(start > 0 && stop > 0, "Cannot compute logspace from negative exponents");
    start = std::log(start)/std::log(base);
    stop = std::log(stop)/std::log(base);
    ret_type ret(N);
    T dx;
    if(endpoint)
    {
        dx = (stop-start)/(N-1);
    }
    else
    {
        dx = (stop-start)/N;
    }

    for(size_t i = 0; i < N; ++i)
    {
        ret[i] = std::pow(base, start + dx*T(i));
    }
    return ret;
}

template <typename T, typename ret_type = linalg::vector<T>>
ret_type softmspace(T start, T stop, size_t N, T beta = 1, bool endpoint = true)
{
    ASSERT(N > 1, "Cannot compute logspace with one point.");
    ASSERT(start > 0 && stop > 0, "Cannot compute logspace from negative exponents");

    ret_type ret(N);
    T dx;

    start = std::log(std::exp(beta*start)-1)/beta;
    stop = std::log(std::exp(beta*stop)-1)/beta;
    if(endpoint)
    {
        dx = (stop-start)/(N-1);
    }
    else
    {
        dx = (start-stop)/N;
    }

    for(size_t i = 0; i < N; ++i)
    {
        ret[i] = std::log(std::exp(beta*(start + dx*i))+1)/beta;

    }
    return ret;
}

};

#endif

