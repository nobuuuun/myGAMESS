/*
 * twoc-new.hpp
 *
 *  Created on: Apr 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT_NEW_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT_NEW_HPP_

#include "foreach.hpp"
#include <sstream>
#include <stdexcept>

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/elem.hpp>

#include "rysq-core.hpp"
#include "rysq-types.hpp"
#include "kernel/ri-int.hpp"



namespace rysq {

    namespace kernel {

	struct invalid_doublet : std::runtime_error {
	    static std::string str(const Doublet<Shell> &doublet) {
		std::ostringstream os;
		os << "invalid doublet " << doublet;
		return os.str();
	    }
	    invalid_doublet(const Doublet<Shell> &doublet)
		: std::runtime_error(str(doublet)) {}
	};








	template< template<class> class T >
	kernel::TwoCenterInt<>* twoc_kernel_eri(const Doublet<Shell> &doublet) {
	    type a = type(doublet[0]);
	    type b = type(doublet[1]);

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types)) {		\
		typedef typename					\
		    kernel::TwoCenterEriFind<BOOST_PP_SEQ_ENUM(types)>::type kernel; \
		return new kernel(doublet, new T<kernel::bra>());	\
	    } else

	    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (RI_RYSQ_TYPES)(RI_RYSQ_TYPES)) {
		throw invalid_doublet(doublet);
	    }

#undef ERI

	}//kernel::TwoCenterInt<>* twoc_kernel_eri



	template< template<class,class> class T >
	kernel::TwoCenterInt<>* twoc_kernel_derivative_eri(const Doublet<Shell> &doublet) {
	    type a = type(doublet[0]);
	    type b = type(doublet[1]);

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types)) {		\
		typedef typename					\
		    kernel::TwoCenterDerivativeEriFind<BOOST_PP_SEQ_ENUM(types)>::type kernel; \
		return new kernel(doublet, new T<kernel::bra,kernel::ket>()); \
	    } else

	    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (RI_RYSQ_TYPES)(RI_RYSQ_TYPES)) {
		throw invalid_doublet(doublet);
	    }

#undef ERI

	}//kernel::TwoCenterInt<>* twoc_kernel_derivative_eri




	template< template<class,class> class T >
	kernel::TwoCenterInt<>* twoc_kernel_derivative_overlap(const Doublet<Shell> &doublet) {
	    type a = type(doublet[0]);
	    type b = type(doublet[1]);

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types)) {		\
		typedef typename					\
		    kernel::TwoCenterDerivativeOverlapFind<BOOST_PP_SEQ_ENUM(types)>::type kernel; \
		return new kernel(doublet, new T<kernel::bra,kernel::ket>()); \
	    } else

	    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (RYSQ_TYPES)(RYSQ_TYPES)) {
		throw invalid_doublet(doublet);
	    }

#undef ERI

	}//kernel::TwoCenterInt<>* twoc_kernel_derivative_overlap












	struct invalid_triplet : std::runtime_error {
	    static std::string str(const Triplet<Shell> &triplet) {
		std::ostringstream os;
		os << "invalid triplet " << triplet;
		return os.str();
	    }
	    invalid_triplet(const Triplet<Shell> &triplet)
		: std::runtime_error(str(triplet)) {}
	};


	template< template<class> class T >
	kernel::ThreeCenterInt<>* threec_kernel_eri(const Triplet<Shell> &triplet) {
	    type a = type(triplet[0]);
	    type b = type(triplet[1]);
	    type c = type(triplet[2]);

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types) &&		\
			  c == BOOST_PP_SEQ_ELEM(2, types)) {		\
		typedef typename					\
		    kernel::ThreeCenterEriFind<BOOST_PP_SEQ_ENUM(types)>::type kernel; \
		return new kernel(triplet, new T<kernel::bra>());	\
	    } else

	    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (RYSQ_TYPES)(RYSQ_TYPES)(RI_RYSQ_TYPES)) {
		throw invalid_triplet(triplet);
	    }

#undef ERI

	}//kernel::ThreeCenterInt<>* threec_kernel_eri






	template< template<class, class> class T >
	kernel::ThreeCenterInt<>* threec_kernel_derivative_eri(const Triplet<Shell> &triplet) {
	    type a = type(triplet[0]);
	    type b = type(triplet[1]);
	    type c = type(triplet[2]);

#define ERI(r, types) if (a == BOOST_PP_SEQ_ELEM(0, types) &&		\
			  b == BOOST_PP_SEQ_ELEM(1, types) &&		\
			  c == BOOST_PP_SEQ_ELEM(2, types)) {		\
		typedef typename					\
		    kernel::ThreeCenterDerivativeEriFind<BOOST_PP_SEQ_ENUM(types)>::type kernel; \
		return new kernel(triplet, new T<kernel::bra, kernel::ket>()); \
	    } else

	    BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (RYSQ_TYPES)(RYSQ_TYPES)(RI_RYSQ_TYPES)) {
		throw invalid_triplet(triplet);
	    }

#undef ERI

	}//kernel::ThreeCenterInt<>* threec_kernel_eri







    } // namespace kernel

} // namespace rysq

#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT_NEW_HPP_ */
