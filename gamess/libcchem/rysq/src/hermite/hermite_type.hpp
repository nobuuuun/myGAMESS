/*
 * hermite_type.hpp
 *
 *  Created on: Sep 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_HERMITE_HERMITE_TYPE_HPP_
#define LIBCCHEM_RYSQ_SRC_HERMITE_HERMITE_TYPE_HPP_

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#ifndef TYPE__
#error "type is not defined"
#endif

namespace BOOST_PP_CAT(TYPE__ , __) {
#include "hermite/polynomials.hpp"
}

//template<>
//struct polynomial<0,TYPE__> {
//    static void evaluate(TYPE__ (&t2)[1], TYPE__ (&W)[1]) {
//	BOOST_PP_CAT(TYPE__ , __)::hermite0::evaluate(t2, W);
//    }
//};

#define HERMITE_N(R, N, DATA)					\
    template<>								\
    struct polynomial<N,TYPE__> {						\
    static void evaluate(TYPE__ (&t2)[N], TYPE__ (&W)[N]) {	\
	BOOST_PP_CAT(TYPE__ , __)::					\
	    BOOST_PP_CAT(hermite, N)::evaluate(t2, W);			\
    }									\
};

//BOOST_PP_SEQ_FOR_EACH(HERMITE_N, (), (1)(2)(3)(4))

BOOST_PP_REPEAT_FROM_TO(1, 16, HERMITE_N, ())

#endif /* LIBCCHEM_RYSQ_SRC_HERMITE_HERMITE_TYPE_HPP_ */
