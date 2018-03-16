/*
 * twoc-eri.hpp
 *
 *  Created on: Apr 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT_HPP_




#include "ri-rysq-int.hpp"
#include "meta.hpp"


namespace rysq {

namespace kernel {

template<class bra = void> //primary template, defaults parameter set to void
struct TwoCenterTransform; //forward declaration, empty base class

template<> //enclosing struct template
struct TwoCenterTransform<> //explicit specialization (no template parameters)
{
	struct Data {};
};

// abstract base class (ABC) is a class with one or more pure virtual functions (ABCs cannot be instantiated)
//  a derived class must define ALL the pure virtual members (provide a body), otherwise the derived class is an ABC too.
template<class bra> // TwoCenterTransform definition when template argument is NOT void
//struct TwoCenterTransform<bra> {
	struct TwoCenterTransform { //primary template
	typedef typename TwoCenterTransform<>::Data Data;
	virtual ~TwoCenterTransform() {}
	virtual TwoCenterTransform& operator()(Data &data) = 0;
	virtual void operator()(int kl,	const double *Q, double scale) = 0;
};
// since all functions are virtual -> interface

template<class A = void, class B = void, class enable = void>
struct TwoCenterInt {
	static const bool value = false;
};

template<>
struct TwoCenterInt<> {
	typedef TwoCenterTransform<>::Data Data;
	TwoCenterInt(const Doublet<Shell> &doublet) : doublet_(doublet) {}
	virtual ~TwoCenterInt() {}
	const Doublet<Shell> &doublet() const { return doublet_; }
	virtual void operator()(const Doublet<Center> &r, Data &data,
			const rysq::TwoCenterInt::Parameters &parameters) = 0;
protected:
	const Doublet<Shell> doublet_;
};



} // namespace kernel
} // namespace rysq







//
//for three-center integrals
//
namespace rysq {
namespace kernel {

template<class bra = void, class ket = void>
struct ThreeCenterTransform;

template<>
struct ThreeCenterTransform<> {
	struct Data {};
};

template<class bra>
struct ThreeCenterTransform<bra> {
    typedef typename ThreeCenterTransform<>::Data Data;
    virtual ~ThreeCenterTransform() {}
    virtual ThreeCenterTransform& operator()(Data &data) = 0;
    virtual void operator()(int k,
			    const double *Q, double scale) = 0;
    virtual void operator()() = 0;
};

template<class A = void, class B = void, class enable = void>
struct ThreeCenterInt {
	static const bool value = false;
};

template<>
struct ThreeCenterInt<> {
	typedef ThreeCenterTransform<>::Data Data;
	ThreeCenterInt(const Triplet<Shell> &triplet) : triplet_(triplet) {}
	virtual ~ThreeCenterInt() {}
	const Triplet<Shell> &triplet() const { return triplet_; }
	virtual void operator()(const Triplet<Center> &r, Data &data,
			const rysq::ThreeCenterInt::Parameters &parameters) = 0;
protected:
	const Triplet<Shell> triplet_;
};


} // namespace kernel
} // namespace rysq












#include "kernel/ri-int1.hpp"


namespace rysq {
namespace kernel {

template<type T0, type T1>
struct TwoCenterEriFind {
	typedef typename meta::OneCenterState<T0> onec_bra; //one-center state
	typedef typename meta::OneCenterState<T1> onec_ket; //one-center state
	typedef boost::mpl::int_<(onec_bra::L + onec_ket::L)/2 + 1> roots;
	typedef TwoCenterEri<onec_bra,roots> type;
};//TwoCenterEriFinD

template<type T0, type T1>
struct TwoCenterDerivativeEriFind {
	typedef typename meta::OneCenterState<T0> onec_bra; //one-center state
	typedef typename meta::OneCenterState<T1> onec_ket; //one-center state
//	typedef boost::mpl::int_<(onec_bra::L + onec_ket::L)/2 + 1> roots;
	// needs higher moments ERIs, bump up the bra moment
	typedef boost::mpl::int_<(onec_bra::L +1 + onec_ket::L)/2 + 1> roots;
	typedef TwoCenterDerivativeEri<onec_bra,onec_ket,roots> type;
};//TwoCenterDerivativeEriFind

template<type T0, type T1>
struct TwoCenterDerivativeOverlapFind {
	typedef typename meta::OneCenterState<T0> onec_bra; //one-center state
	typedef typename meta::OneCenterState<T1> onec_ket; //one-center state
	typedef boost::mpl::int_<(onec_bra::L +1 + onec_ket::L)/2 + 1> npoints;
	typedef TwoCenterDerivativeOverlap<onec_bra,onec_ket,npoints> type;
};//TwoCenterDerivativeEriFind



template<type T0, type T1, type T2>
struct ThreeCenterEriFind {
	typedef typename meta::state<T0,T1> bra; //two-center state
	typedef typename meta::OneCenterState<T2> onec_ket; //one-center state
	typedef boost::mpl::int_<(bra::L + onec_ket::L)/2 + 1> roots;
	typedef ThreeCenterEri<bra,roots> type;
};//ThreeCenterEriFind

template<type T0, type T1, type T2>
struct ThreeCenterDerivativeEriFind {
	typedef typename meta::state<T0,T1> bra; //two-center state
	typedef typename meta::OneCenterState<T2> onec_ket; //one-center state
	// needs higher moments ERIs, bump up the bra moment by one (total)
	typedef boost::mpl::int_<(bra::L + onec_ket::L +1)/2 + 1> roots;

	typedef ThreeCenterDerivativeEri<bra, onec_ket, roots> type;
};//ThreeCenterDerivativeEriFind


} // namespace kernel
} // namespace rysq

#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT_HPP_ */
