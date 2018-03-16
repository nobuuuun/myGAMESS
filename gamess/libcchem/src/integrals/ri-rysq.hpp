/*
 * ri-rysq.hpp
 *
 *  Created on: Aug 14, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_INTEGRALS_RI_RYSQ_HPP_
#define LIBCCHEM_SRC_INTEGRALS_RI_RYSQ_HPP_


#include "basis/basis.hpp"
#include "integrals/quartet.hpp"
#include "exception.hpp"


#include <ri-rysq.hpp>

#include "adapter/rysq.hpp"

#include <boost/array.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include "boost/utility/profiler.hpp"



namespace integrals {
namespace rysq {


template<typename T>
struct TwoCenterInt {

public:
	TwoCenterInt(Doublet<const Basis::Shell::Data&> doublet) :
		A_(doublet.q), B_(doublet.s),
		twoc_int_(::rysq::Doublet< ::rysq::Shell >(A_, B_)) {}

	void operator()(const std::vector<Basis::Center> &centers,
			const std::vector< boost::array<int,2> > &doublet,
			double *G, double cutoff) {
		twoc_int_(centers, doublet, G, cutoff);
	}

private:
	adapter::rysq::Shell A_, B_;
	T twoc_int_; //::rysq::TwoCenterEri or ::rysq::TwoCenterDerivativeEri
};//TwoCenterInt






template<typename T>
struct ThreeCenterInt {
	ThreeCenterInt(Triplet<const Basis::Shell::Data&> triplet) :
		A_(triplet.q), B_(triplet.s),C_(triplet.l),
		threec_eri_(::rysq::Triplet< ::rysq::Shell >(A_, B_, C_)) {}

	void operator()(const std::vector<Basis::Center> &obscenters,
			const std::vector<Basis::Center> &auxcenters,
			const std::vector< boost::array<int,3> > &triplet,
			double *G, double cutoff) {
		threec_eri_(obscenters, auxcenters, triplet, G, cutoff);
	}

private:
	adapter::rysq::Shell A_, B_, C_;
	T threec_eri_;

};//ThreeCenterInt



} // rysq
} // integrals



#endif /* LIBCCHEM_SRC_INTEGRALS_RI_RYSQ_HPP_ */
