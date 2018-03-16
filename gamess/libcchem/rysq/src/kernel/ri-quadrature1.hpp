/*
 * twoc-quadrature1.hpp
 *
 *  Created on: Apr 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_TWOC_QUADRATURE1_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_TWOC_QUADRATURE1_HPP_

#include "kernel/ri-eri-quadrature.hpp"
#include "kernel/twoc-quadrature-impl.hpp"

#include <math.h>
#include "rysq-core.hpp"
#include "vector.hpp"

#include "roots/roots.hpp"
#include "kernel/ri-twoc-recurrence.hpp"

namespace rysq {
    namespace kernel {

	namespace twoc_quadrature {

	    //declaration for second function below
	    template<class bra, size_t N, size_t NT,
		     class Transform, typename T, typename U>
	    static void apply(rysq::type c,
			      const Primitives<T,U> &primitives,
			      Transform &transform, double scale, bool debug);

	    template<class bra, size_t N,template<typename, size_t> class align,
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

		double ei[a.K];
		for (int Ki = 0; Ki < a.K; ++Ki) {
		    ei[Ki] = exp(-a(Ki));
		} //Ki


		int K = 0;
		for(int Kj = 0; Kj < b.K; ++Kj) { //loop over primitives  : normally 1 for RI basis
		    double aj = b(Kj);

		    double B = aj;

		    vector<3> rB = rj;

		    double Cj[2] __attribute__ ((aligned(16)));
		    for(int j = 0; j < b.nc; ++j){
			Cj[j] = b(Kj,j)/B;
		    }

		    for(int Ki = 0; Ki < a.K; ++Ki){ //primitives
			double ai = a(Ki);
			double A = ai;

			double Ci[bra::nc] __attribute__ ((aligned(16)));

			for(int i = 0; i < bra::nc; ++i){
			    Ci[i] = a(Ki,i)/A;
			}

			for(int j = 0; j < b.nc; ++j){
			    for(int i = 0; i < a.nc; ++i){
				primitives.C[j][i +K*bra::nc] = Ci[i]*Cj[j];
			    }
			}

			vector<3> rA = ri;
			vector<3> rAB = rB - rA;

			double rho = A*B/(A+B);
			U X = rho*rAB.dot();
			vector<N,U> W, t2;
			rysq::roots<N>(X, t2.elems, W.elems);
			t2 /= (A + B);

			for(int i = 0; i<N; i++){
			    W.elems[i] /= sqrt(A+B);
			}

			T *Ix = primitives.Ix(K);
			T *Iy = primitives.Iy(K);
			T *Iz = primitives.Iz(K);

			if( doublet.L() ) {
			    typedef kernel::twoc_recurrence<bra::L,N,align> twoc_recurrence;
			    twoc_recurrence::apply(b.L, A, B, rAB, rA, rB,
						   t2.data(), W.data(), Ix, Iy, Iz);
			}else{
			    for(int a = 0; a < N; a++){
				Ix[a] = Iy[a] = 1.0;
				Iz[a] = W[a];
			    }
			}

		    }//Ki

		} //Kj

		primitives.K=1;
		apply<bra,N,NT>(b, primitives, transform, scale,debug);

	    }

	template<class bra, size_t N, size_t NT,
		 class Transform, typename T, typename U>
	static void apply(rysq::type b,
			  const Primitives<T,U> &primitives,
			  Transform &transform, double scale,bool debug) {

	    static const int Ni = NT*(bra::A::L+1);

	    const int Lb = abs(b);
	    const int dim2d = Ni*(Lb+1);

	    const int K = primitives.K;
	    const double* const *C = primitives.C;
	    const T *Ix = primitives.Ix(0);
	    const T *Iy = primitives.Iy(0);
	    const T *Iz = primitives.Iz(0);

	    const int b_first = shell(b).begin();
	    const int b_last = shell(b).end() - 1;

	    for(int b = b_first, kl = 0; b <= b_last; ++b, ++kl) {

		const double *Ckl = C[0];

		const int klx = Ni*LX[b];
		const int kly = Ni*LY[b];
		const int klz = Ni*LZ[b];

		static const int ni = bra::A::size;
		double I[ni] __attribute__((aligned(16))) = { 0.0 };

		int flags = 0;
		double screen = 0.0;

		typedef twoc_quadrature::impl<bra> impl;
		impl::template apply<N,T,NT>(flags, screen, K, Ckl, dim2d,
					     &Ix[klx], &Iy[kly], &Iz[klz], 1.0, I);

		transform(kl, I, scale*NORMALIZE[b]); //kernel/twoc-eri-transform.hpp

	    }

	}

    } // namespace twoc_quadrature


}//kernel
}//rysq


#include "kernel/ri-threec-recurrence.hpp"
#include "kernel/transfer.hpp"
#include "kernel/quadrature1-impl.hpp"

#include "normalize.h"
namespace rysq {
namespace kernel {


namespace threec_quadrature {


//const double ANORMALIZE[] = {
//	RI_RYSQ_NORMAL_000,
//	RI_RYSQ_NORMAL_100, RI_RYSQ_NORMAL_100, RI_RYSQ_NORMAL_100,
//	RI_RYSQ_NORMAL_200, RI_RYSQ_NORMAL_200, RI_RYSQ_NORMAL_200,
//	RI_RYSQ_NORMAL_110, RI_RYSQ_NORMAL_110, RI_RYSQ_NORMAL_110,
//	RI_RYSQ_NORMAL_300, RI_RYSQ_NORMAL_300, RI_RYSQ_NORMAL_300,
//	RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210,
//	RI_RYSQ_NORMAL_111,
//	RI_RYSQ_NORMAL_400, RI_RYSQ_NORMAL_400, RI_RYSQ_NORMAL_400,
//	RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310,
//	RI_RYSQ_NORMAL_220, RI_RYSQ_NORMAL_220, RI_RYSQ_NORMAL_220,
//	RI_RYSQ_NORMAL_211, RI_RYSQ_NORMAL_211, RI_RYSQ_NORMAL_211,
//	RI_RYSQ_NORMAL_500, RI_RYSQ_NORMAL_500, RI_RYSQ_NORMAL_500,
//	RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410,
//	RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320,
//	RI_RYSQ_NORMAL_311, RI_RYSQ_NORMAL_311, RI_RYSQ_NORMAL_311,
//	RI_RYSQ_NORMAL_221, RI_RYSQ_NORMAL_221, RI_RYSQ_NORMAL_221};

const double *NORMALIZE = &(::rysq::kernel::twoc_quadrature::NORMALIZE[0]);

//declaration for second function
template<class bra, size_t N, size_t NT,
class Transform, typename T, typename U>
static void apply(rysq::type c,
		const Primitives<T,U> &primitives,
		Transform &transform, double scale, bool debug);

template<class bra, size_t N,template<typename, size_t> class align,
class Transform, typename T, typename U>
static void apply(const Triplet<Shell> &triplet,
		const vector<3> &ri, const vector<3> &rj, const vector<3> &rk,
		double scale, double cutoff,
		Primitives<T,U> &primitives,
		Transform &transform) {

	static const size_t NT = N + align<T,N>::value;

	const Shell &a = triplet[0];
	const Shell &b = triplet[1];
	const Shell &c = triplet[2];

	bool debug = 0;

	vector<3> rij = ri - rj;

	const double rij2 = rij.dot();

	double eij[a.K*b.K];
	for (int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
	    for (int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
		const double A = a(Ki) + b(Kj);
		const double A1 = 1.0/A;
		eij[Kij] = exp(-a(Ki)*b(Kj)*A1*rij2);
	    }
	}

	int K = 0;

	for (int Kk = 0; Kk < c.K; ++Kk) {
	    const double ak = c(Kk);
	    const double B = ak;
	    
	    vector<3> rB = rk;
	    
		// double Ck[4] __attribute__ ((aligned(16)));
		// double Ck_max = 0.0;
		// for(int k = 0; k < c.nc; ++k) {
		// 	Ck[k] = c(Kk,k);
		// 	Ck_max = std::max(Ck_max, fabs(Ck[k]));
		// }
	    double Ck[1] __attribute__ ((aligned(16)));
		//		double Ck_max = 0.0;
		//for(int k = 0; k < c.nc; ++k) {
	    Ck[0] = c(Kk,0);
	    const double Ck_max = fabs(Ck[0]); //std::max(Ck_max, fabs(Ck[k]));
			//}

	    for(int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
		for(int Ki = 0; Ki < a.K; ++Ki, ++Kij) {

		    const double ai = a(Ki);
		    const double aj = b(Kj);
		    const double A = ai + aj;

		    const double test = ai*aj/A*rij2;
		    if ( test > (double)23.0258)continue;

		    const double e = eij[Ki+Kj*a.K];;
		    const double eAB = e/(A*B*sqrt(A+B));

		    double Cij[bra::nc] __attribute__ ((aligned(16)));
		    double Cij_max = 0.0;

		    for(int j = 0, ij = 0; j < bra::B::nc; ++j) {
			for(int i = 0; i < bra::A::nc; ++i, ++ij) {
			    Cij[ij] = eAB*a(Ki,i)*b(Kj,j);
			    Cij_max = std::max(fabs(Cij[ij]), Cij_max);
			}
		    }

		    //for(int k = 0; k < c.nc; ++k) {
		    for(int ij = 0; ij < bra::nc; ++ij) {
			primitives.C[0][ij + K*bra::nc] = Cij[ij]*Ck[0];
		    }
		    //}

		    vector<3> rA = vector<3>::center(ai, ri, aj, rj);
		    vector<3> rAi = rA - ri;
		    vector<3> rAB = rA - rB;

		    const double rho = A*B/(A + B);
		    U X =  rho*rAB.dot();
		    
		    vector<N,U> W, t2;
		    rysq::roots<N>(X, t2.elems, W.elems);
		    t2 /= (A + B);

		    static const int LA = bra::A::L;
		    static const int LB = bra::B::L;
		    
		    T *Ix = primitives.Ix(K);
		    T *Iy = primitives.Iy(K);
		    T *Iz = primitives.Iz(K);

		    typedef kernel::threec_recurrence<bra::L,N,align> recurrence;
		    typedef kernel::transfer<LA,LB,N,align> transfer;

		    if (LB == 0) {  //nothing to transfer -> if (LB == 0) {
			recurrence::apply(LA, c.L, A, B, rAB, rAi, rk,
					  t2.data(), W.data(), Ix, Iy, Iz);
			
		    }else{
			U *Gx = primitives.Gx(K);
			U *Gy = primitives.Gy(K);
			U *Gz = primitives.Gz(K);
			
			double * tmp;
			
			recurrence::apply(LA, c.L, A, B, rAB, rAi, rk,
					  t2.data(), W.data(), Gx, Gy, Gz);
			
			transfer::apply(c.L, 0, U(rij[0]),U(rij[0]), Gx, Ix,tmp);
			transfer::apply(c.L, 0, U(rij[1]),U(rij[1]), Gy, Iy,tmp);
			transfer::apply(c.L, 0, U(rij[2]),U(rij[2]), Gz, Iz,tmp);
		    }

		    ++K;
		    
		    if (K == primitives.Kmax) {
			primitives.K = K;
			apply<bra,N,NT>(c, primitives, transform, scale,debug);
			K = 0;
		    }
		    
		}//Ki :A
	    }//Kj :B
	}//Kk :C
	
	// Final contraction
	primitives.K = K;
	if (primitives.K) apply<bra,N,NT>(c, primitives, transform, scale,debug);
	
}// apply



    template<class bra, size_t N, size_t NT,
	     class Transform, typename T, typename U>
    static void apply(rysq::type c,
		      const Primitives<T,U> &primitives,
		      Transform &transform, double scale, bool debug) {
	
	static const int Nij = NT*(bra::A::L+1)*(bra::B::L+1);
	
	const int Lc = abs(c);
	
	const int dim2d = Nij*(Lc+1);
	
	const int K = primitives.K;
	const double* const *C = primitives.C;
	const T *Ix = primitives.Ix(0);
	const T *Iy = primitives.Iy(0);
	const T *Iz = primitives.Iz(0);
	
	const int c_first = shell(c).begin();
	const int c_last = shell(c).end() - 1;
	
	for(int k = c_first; k <= c_last; ++k) {
	    const double *Ckl = C[0];
	    
	    const int kx = Nij*(LX[k]);
	    const int ky = Nij*(LY[k]);
	    const int kz = Nij*(LZ[k]);
	    
	    static const int ni = bra::A::size;
	    static const int nj = bra::B::size;
	    double I[ni*nj] __attribute__((aligned(16))) = { 0.0 };
	    
	    int flags = 0;
	    double screen = 0.0;
	    
	    typedef quadrature::impl<bra> impl;
	    impl::template apply<N,T,NT>(flags, screen, K, Ckl, dim2d,
					 &Ix[kx], &Iy[ky], &Iz[kz], 1.0, I);
	    
	    transform(k - c_first, I, scale*NORMALIZE[k]);
	    
	}//k
	
    }//apply
    
    
    
} // namespace threec_quadrature





}//kernel
}//rysq

#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_TWOC_QUADRATURE1_HPP_ */
