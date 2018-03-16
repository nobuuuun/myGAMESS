/*
 * hermite.hpp
 *
 *  Created on: Sep 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_HERMITE_HERMITE_HPP_
#define LIBCCHEM_RYSQ_SRC_HERMITE_HERMITE_HPP_

#ifdef __CUDACC__
#define constant__ __device__ const
#else
#define constant__ static const
#endif

namespace hermite {
namespace detail {

template<size_t N, typename T>
struct polynomial;

#define TYPE__ double
#include "hermite/hermite_type.hpp"
#undef TYPE__
}}

namespace hermite {

template<size_t N, typename T>
typename boost::enable_if_c<(N == 0)>::type
polynomial(T (&t2)[1], T (&W)[1]) {
	return detail::polynomial<N,T>::evaluate(t2, W);
}

template<size_t N, typename T>
void
polynomial(T (&t2)[N], T (&W)[N]) {
	return detail::polynomial<N,T>::evaluate(t2, W);
}

}//namespace hermite



namespace hermite {

namespace detail {

template<size_t N, typename T>
struct polynomial{
	static void evaluate(T *t2, T *W){}

};//struct hermite



} //namespace detail

} //namespace hermite


#endif /* LIBCCHEM_RYSQ_SRC_HERMITE_HERMITE_HPP_ */
