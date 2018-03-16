/*
 * derivative-overlap-qaudrature.hpp
 *
 *  Created on: Sep 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_OVERLAP_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_OVERLAP_HPP_

#include "kernel/derivative-overlap-quadrature-impl.hpp"
//#include "kernel/derivative-twoc-quadrature-impl.hpp"

#include <math.h>
#include <cmath>
#include "rysq-core.hpp"
#include "vector.hpp"

#include "hermite/hermite.hpp"
#include "kernel/overlap-recurrence.hpp"

namespace rysq {
    namespace kernel {


	namespace twoc_derivative_overlap_quadrature {

	    //declaration for second function below
	    template<class bra, size_t N, size_t NT,
		     class Transform, typename T, typename U>
	    static void apply(rysq::type c,
			      const Primitives<T,U> &primitives,
			      Transform &transform, double scale, double alpha, bool debug, double *vfac);

	    template<class bra, class ket, size_t N,template<typename, size_t> class align,
		     class Transform, typename T, typename U>
	    static void apply(const Doublet<Shell> &doublet,
			      const vector<3> &ri, const vector<3> &rj,
			      double scale, double cutoff,
			      Primitives<T,U> &primitives,
			      Transform &transform) {

		static const size_t NT = N + align<T,N>::value;

		const Shell &a = doublet[0];
		const Shell &b = doublet[1];

		bool debug = 0;

		vector<3> rr = ri - rj;
		double RR = rr.dot();

		double vfac[primitives.Kmax] __attribute__ ((aligned(16)));

		int K = 0;
		for(int Kj = 0; Kj < b.K; ++Kj) { //loop over primitives  : normally 1 for RI basis
		    double aj = b(Kj);
		    double arrj = aj*RR;
		    double axj = aj*rj[0];
		    double ayj = aj*rj[1];
		    double azj = aj*rj[2];


		    double Cj[2] __attribute__ ((aligned(16)));
		    for(int j = 0; j < b.nc; ++j){  // (nc = 2 --> for sp-hybrid)
			Cj[j] = b(Kj,j);
		    }

		    for(int Ki = 0; Ki < a.K; ++Ki){ //primitives
			double ai = a(Ki);
			double aa = ai+aj;
			double aa1 = 1/aa;
			double dum = ai*arrj*aa1;
			double fac = exp(-dum);
			double ax = (axj + ai*ri[0])*aa1;
			double ay = (ayj + ai*ri[1])*aa1;
			double az = (azj + ai*ri[2])*aa1;

			vfac[K] = fac;

			double Ci[bra::nc] __attribute__ ((aligned(16)));
			for(int i = 0; i < bra::nc; ++i){ // (nc = 2 --> for sp-hybrid)
			    Ci[i] = a(Ki,i);
			}

			for(int j = 0; j < b.nc; ++j){
			    for(int i = 0; i < a.nc; ++i){
				primitives.C[j][i +K*bra::nc] = Ci[i]*Cj[j];
			    }
			}

			double aa12 = sqrt(aa1);

			T *Ix = primitives.Ix(K);
			T *Iy = primitives.Iy(K);
			T *Iz = primitives.Iz(K);

			//set exponents for centers A and B
			T *AExp = primitives.EAlpha(K);
			*AExp = 2*a(Ki);

			for(int lj1 = 0; lj1 < ket::L+1; lj1++){
			    for(int li1 = 0; li1 < bra::L+2; li1++){

				int npts = (li1+lj1)/2+1;
				double x0 = ax;
				double y0 = ay;
				double z0 = az;

#define HERMITE_CASE(R,CASE,DATA)					\
				if(npts == CASE){			\
				    vector<CASE,U> Wp, t2p;		\
				    hermite::polynomial<CASE>(t2p.elems, Wp.elems); \
				    typedef overlap_recurrence<CASE,bra::L+2> s_recurrence; \
				    s_recurrence::apply( li1, lj1, aa12, x0, y0, z0, \
							 ri.data(), rj.data(), t2p.data(), Wp.data(), \
							 Ix,Iy,Iz);	\
				};
				
				BOOST_PP_REPEAT_FROM_TO(1, 16, HERMITE_CASE, ())
				    
				    }//lj1
			}//li1

			++K;

		    // primitives.K = K;
		    // apply<bra,1,NT>(b, primitives, transform, scale, debug, &vfac[0] );
		    // K = 0;

		    }//ki
		}//kj

		// contract primitives
		//		if (K == primitives.Kmax) {
		    primitives.K = K;
		    apply<bra,1,NT>(b, primitives, transform, scale, debug, &vfac[0] );
		    K = 0;
		    //		}

	    } //apply






	template<class bra, size_t N, size_t NT,
		 class Transform, typename T, typename U>
	static void apply(rysq::type b,
			  const Primitives<T,U> &primitives,
			  Transform &transform, double scale, bool debug, double *vfac ) {

	    const int K = primitives.K;
	    static const int Ni = bra::A::L+2;


	    const int Lb = abs(b);
	    const int dim2d = Ni*(Lb+1);

	    const int b_first = shell(b).begin();
	    const int b_last = shell(b).end() - 1;
	    const double* const *C = primitives.C;
	    const T *Ix = primitives.Ix(0);
	    const T *Iy = primitives.Iy(0);
	    const T *Iz = primitives.Iz(0);

	    //A+B center exponents
	    const T *AExp = primitives.EAlpha(0);

	    //is b an sp function?
	    int spb = (b < 0);

	    //loop moment b (eg. dxx,dyy,dzz,dxy,dxz,dyz)
	    for(int b = b_first, kl = 0; b <= b_last; ++b, ++kl) { 


		const int klx = Ni*LX[b];
		const int kly = Ni*LY[b];
		const int klz = Ni*LZ[b];
		const double *Ckl = C[(spb && b)];

		int flags = 0;
		double screen = 0.0;

		double I[bra::A::size*3] __attribute__((aligned(16))) = { 0.0 };


		typedef twoc_derivative_overlap_quadrature::impl<bra> impl;
		impl::template apply<1,T,NT>(flags, screen, K, Ckl, AExp, dim2d,
		 			     &Ix[klx], &Iy[kly], &Iz[klz], 1.0, I, &vfac[0] );


		//typedef twoc_derivative_quadrature::impl<bra> impl;
		//bool normalize=0;
		//impl::template apply<1,T,NT>(normalize, screen, K, Ckl, AExp, dim2d,
		//&Ix[klx], &Iy[kly], &Iz[klz], 1.0, I, &vfac[0] );

		transform(kl, I, (double)1);

	    } //b_first

	} //apply


    } //namespace twoc_derivative_quadrature
} //namespace kernel

    }//namespace hermite


#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_OVERLAP_HPP_ */
