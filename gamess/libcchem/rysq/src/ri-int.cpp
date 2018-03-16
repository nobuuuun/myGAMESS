/*
 * twoc-eri.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: luke
 */


#include "ri-rysq-int.hpp"
#include "kernel/ri-int-new.hpp"
#include "ri-int-transform.hpp"

#include <boost/typeof/typeof.hpp>

using namespace rysq;

typedef kernel::TwoCenterInt<> twoc_kernel_type;

TwoCenterInt::~TwoCenterInt() { delete static_cast<twoc_kernel_type*>(kernel_); }

void TwoCenterInt::operator()(const Doublet<Center> &centers,
		double *Q, const Parameters &parameters) {
	twoc_transform::TwoCenterTransform<>::Data data(Q);
	Parameters p = parameters;
	(*static_cast<twoc_kernel_type*>(kernel_))(centers, data, p);
}

void TwoCenterInt::operator()(const std::vector<Center> &centers,
		const std::vector<Int2> &doublets,
		double *Q,const Parameters &parameters) {

	foreach (const Int2 &index, doublets) {
		(*this)(centers, index, Q, parameters);
		Q += size_;
	}
}

//two-center eri type
TwoCenterEri::TwoCenterEri(const Doublet<Shell> &doublet){
	kernel_ = kernel::twoc_kernel_eri<twoc_transform::TwoCenterEriTransform>(doublet);
	size_ = doublet.size();
}

//two-center derivative eri type
TwoCenterDerivativeEri::TwoCenterDerivativeEri(const Doublet<Shell> &doublet){
	kernel_ = kernel::twoc_kernel_derivative_eri<twoc_transform::TwoCenterDerivativeEriTransform>(doublet);
	size_ = doublet.size();
}

//two-center derivative eri type
TwoCenterDerivativeOverlap::TwoCenterDerivativeOverlap(const Doublet<Shell> &doublet){
	kernel_ = kernel::twoc_kernel_derivative_overlap<twoc_transform::TwoCenterDerivativeEriTransform>(doublet);
	size_ = doublet.size();
}








typedef kernel::ThreeCenterInt<> threec_kernel_type;

ThreeCenterInt::~ThreeCenterInt() { delete static_cast<threec_kernel_type*>(kernel_); }

void ThreeCenterInt::operator()(const Triplet<Center> &centers,
		double *Q, const Parameters &parameters) {
	threec_transform::ThreeCenterTransform<>::Data data(Q);
	Parameters p = parameters;
	(*static_cast<threec_kernel_type*>(kernel_))(centers, data, p);
}

void ThreeCenterInt::operator()(const std::vector<Center> &obscenters,
		const std::vector<Center> &auxcenters,
		const std::vector<Int3> &triplets,
		double *Q,const Parameters &parameters) {

	foreach (const Int3 &index, triplets) {
		(*this)(obscenters, auxcenters, index, Q, parameters);
		Q += size_;
	}
}


//three-center eri type
ThreeCenterEri::ThreeCenterEri(const Triplet<Shell> &triplet){
	kernel_ = kernel::threec_kernel_eri<threec_transform::ThreeCenterEriTransform>(triplet); //this one, i think transpose is not needed
	size_ = triplet.size();
}


//three-center derivative eri type
ThreeCenterDerivativeEri::ThreeCenterDerivativeEri(const Triplet<Shell> &triplet){
	kernel_ = kernel::threec_kernel_derivative_eri<threec_transform::ThreeCenterDerivativeEriTransform>(triplet); //this one, i think transpose is not needed
	size_ = triplet.size();
}




