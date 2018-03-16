#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT1_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT1_HPP_


#include <math.h>
#include <boost/mpl/assert.hpp>

#include "rysq-core.hpp"
#include "ri-rysq-core.hpp"
#include "vector.hpp"

#include "kernel/ri-quadrature1.hpp"
#include "kernel/ri-derivative-quadrature1.hpp"
#include "kernel/derivative-overlap-quadrature.hpp"

namespace rysq {
namespace kernel {

static const double CUTOFF_SCALE = 1e-2;

template<typename T, size_t N>
struct align {

#define ALIGN_(n) ((A/sizeof(T) - (n)%(A/sizeof(T)))%(A/sizeof(T)))
    static const size_t A = 16;
    static const size_t value = ALIGN_(N);
    static size_t get(size_t M) { return ALIGN_(M); }
#undef ALIGN_

};//align





template<class T, class N>
struct TwoCenterEri;

template<class bra_, int N>
struct TwoCenterEri<bra_, boost::mpl::int_<N> > : public TwoCenterInt<>
{
	typedef bra_ bra;
	typedef void ket;
	typedef kernel::TwoCenterTransform<bra> Transform;
	typedef TwoCenterInt<>::Data Data;

	TwoCenterEri(const Doublet<Shell> &doublet, Transform *transform)
	: TwoCenterInt<>(doublet), transform_(transform)
	  {
		primitives_.allocate<align>(doublet);
	  }

	~TwoCenterEri() { delete transform_; }

	virtual void operator()(const Doublet<Center> &r, Data &data,
			const rysq::TwoCenterInt::Parameters &parameters) {

		double scale = rysq::SQRT_4PI5;
		double cutoff = parameters.cutoff/(scale*this->doublet_.K());
		cutoff *= CUTOFF_SCALE;
		typedef kernel::vector<3> vector;

		twoc_quadrature::apply<bra,N,align>(this->doublet_,
				vector(r[0]), vector(r[1]),
				scale, cutoff, primitives_,
				(*transform_)(data));
	}

private:
	Transform *transform_;
	twoc_quadrature::Primitives<double, double> primitives_;

};//TwoCenterEri<bra_, boost::mpl::int_<N> >









template<class T, class U, class N>
struct TwoCenterDerivativeEri;

template<class bra_, class ket_, int N>
struct TwoCenterDerivativeEri<bra_, ket_, boost::mpl::int_<N> > : public TwoCenterInt<>
{
	typedef bra_ bra;
//	typedef void ket;
	typedef ket_ ket;
	typedef kernel::TwoCenterTransform<bra> Transform;
	typedef TwoCenterInt<>::Data Data;

	TwoCenterDerivativeEri(const Doublet<Shell> &doublet, Transform *transform)
	: TwoCenterInt<>(doublet), transform_(transform)
	  {
		primitives_.allocate<align>(doublet);
	  }

	~TwoCenterDerivativeEri() { delete transform_; }

	virtual void operator()(const Doublet<Center> &r, Data &data,
			const rysq::TwoCenterInt::Parameters &parameters) {

		double scale = rysq::SQRT_4PI5;
		double cutoff = parameters.cutoff/(scale*this->doublet_.K());
		cutoff *= CUTOFF_SCALE;
		typedef kernel::vector<3> vector;
		twoc_derivative_quadrature::apply<bra,N,align>(this->doublet_,
				vector(r[0]), vector(r[1]),
				scale, cutoff, primitives_,
				(*transform_)(data));
	}

private:
	Transform *transform_;
	twoc_derivative_quadrature::Primitives<double, double> primitives_;

};//TwoCenterDerivativeEri<bra_, boost::mpl::int_<N> >



template<class T, class U, class N>
struct TwoCenterDerivativeOverlap;

template<class bra_, class ket_, int N>
struct TwoCenterDerivativeOverlap<bra_, ket_, boost::mpl::int_<N> > : public TwoCenterInt<>
{
	typedef bra_ bra;
//	typedef void ket;
	typedef ket_ ket;
	typedef kernel::TwoCenterTransform<bra> Transform;
	typedef TwoCenterInt<>::Data Data;

	TwoCenterDerivativeOverlap(const Doublet<Shell> &doublet, Transform *transform)
	: TwoCenterInt<>(doublet), transform_(transform)
	  {
		primitives_.allocate<align>(doublet);
	  }

	~TwoCenterDerivativeOverlap() { delete transform_; }

	virtual void operator()(const Doublet<Center> &r, Data &data,
			const rysq::TwoCenterInt::Parameters &parameters) {

		double scale = rysq::SQRT_4PI5;
		double cutoff = parameters.cutoff/(scale*this->doublet_.K());
		cutoff *= CUTOFF_SCALE;
		typedef kernel::vector<3> vector;
//		std::cout << "what what" << std::endl;
				twoc_derivative_overlap_quadrature::apply<bra,ket,N,align>(this->doublet_,
				vector(r[0]), vector(r[1]),
				scale, cutoff, primitives_,
				(*transform_)(data));
	}

private:
	Transform *transform_;
	twoc_derivative_overlap_quadrature::Primitives<double, double> primitives_;

};//TwoCenterDerivativeOverlap<bra_, boost::mpl::int_<N> >













template<class T, class N>
struct ThreeCenterEri;

template<class bra_, int N>
struct ThreeCenterEri<bra_, boost::mpl::int_<N> > : public ThreeCenterInt<>
{
	typedef bra_ bra;
	typedef void ket;
	typedef kernel::ThreeCenterTransform<bra> Transform;
	typedef ThreeCenterInt<>::Data Data;

	ThreeCenterEri(const Triplet<Shell> &triplet, Transform *transform)
	: ThreeCenterInt<>(triplet), transform_(transform)
	  {
		primitives_.allocate<align>(triplet);
	  }

	~ThreeCenterEri() { delete transform_; }

	virtual void operator()(const Triplet<Center> &r, Data &data,
			const rysq::ThreeCenterInt::Parameters &parameters) {
		double scale = rysq::SQRT_4PI5;
		double cutoff = parameters.cutoff/(scale*this->triplet_.K());
		cutoff *= CUTOFF_SCALE;
		typedef kernel::vector<3> vector;
		threec_quadrature::apply<bra,N,align>(this->triplet_,
				vector(r[0]), vector(r[1]), vector(r[2]),
				scale, cutoff, primitives_,
				(*transform_)(data));
	}

private:
	Transform *transform_;
	threec_quadrature::Primitives<double, double> primitives_;

};//ThreeCenterEri<bra_, boost::mpl::int_<N> >








template<class T, class U, class N>
struct ThreeCenterDerivativeEri;

template<class bra_, class ket_, int N>
struct ThreeCenterDerivativeEri<bra_, ket_, boost::mpl::int_<N> > : public ThreeCenterInt<>
{
	typedef bra_ bra;
//	typedef void ket;
	typedef ket_ ket;
	typedef kernel::ThreeCenterTransform<bra> Transform;
	typedef ThreeCenterInt<>::Data Data;

	ThreeCenterDerivativeEri(const Triplet<Shell> &triplet, Transform *transform)
	: ThreeCenterInt<>(triplet), transform_(transform)
	  {
		primitives_.allocate<align>(triplet);
	  }

	~ThreeCenterDerivativeEri() { delete transform_; }

	virtual void operator()(const Triplet<Center> &r, Data &data,
			const rysq::ThreeCenterInt::Parameters &parameters) {
		double scale = rysq::SQRT_4PI5;
		double cutoff = parameters.cutoff/(scale*this->triplet_.K());
		cutoff *= CUTOFF_SCALE;
		typedef kernel::vector<3> vector;

		threec_derivative_quadrature::apply<bra,N,align>(this->triplet_,
				vector(r[0]), vector(r[1]), vector(r[2]),
				scale, cutoff, primitives_,
				(*transform_)(data));
	}

private:
	Transform *transform_;
	threec_derivative_quadrature::Primitives<double, double> primitives_;

};//ThreeCenterDerivativeEri<bra_, boost::mpl::int_<N> >







} // namespace kernel
} // namespace rysq




#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_RI_INT1_HPP_ */
