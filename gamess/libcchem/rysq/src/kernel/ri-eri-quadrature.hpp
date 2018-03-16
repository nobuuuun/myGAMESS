/*
 * twoc-quadrature.hpp
 *
 *  Created on: Apr 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_TWOC_QUADRATURE_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_TWOC_QUADRATURE_HPP_


#include "rysq-core.hpp"
#include "vector.hpp"


#include <memory>
#include <assert.h>
#include <stdlib.h>
#include <math.h>


namespace rysq {
    namespace kernel {




	namespace twoc_quadrature {


	    template<int N>
	    struct roots {
		static const size_t value = N;
	    };


	    template<class braket>
	    void apply(const Doublet<Shell> &doublet,
		       const vector<3> &ri, const vector<3> &rj,
		       double (&Q)[braket::size], int transpose,
		       double cutoff);

	    template<typename T, typename U>
	    struct Primitives {
		int Kmax, K;
		int size2d;
		void *memory16_;

		T *Ints2d;
		double *C[2];  //coefficients

		template<template<typename, size_t> class align>
		void allocate(const Doublet<Shell> &doublet);

		~Primitives() {
		    if (memory16_) free(memory16_); }

		T* Ix(int K = 0) { return get<T,0>(K); }
		T* Iy(int K = 0) { return get<T,1>(K); }
		T* Iz(int K = 0) { return get<T,2>(K); }

		const T* Ix(int K = 0) const { return get<T,0>(K); }
		const T* Iy(int K = 0) const { return get<T,1>(K); }
		const T* Iz(int K = 0) const { return get<T,2>(K); }

	    private:
		template<typename R, size_t q> R* get(int i) {
		    return reinterpret_cast<R*>(&Ints2d[(3*i + q)*size2d]);
		}
		template<typename R, size_t q>  const R* get(int i)const {
		    return (const_cast<Primitives*>(this))->get<R,q>(i);
		}
	    };//Primitives


	    template<class bra, size_t N, template<typename, size_t> class align,
		     class Transform, typename T, typename U>
	    void apply(const Doublet<Shell> &doublet,
		       const vector<3> &ri, const vector<3> &rj,
		       double scale, double cutoff,
		       Primitives<T,U> &primitives,
		       Transform &transform);

	template<typename T, typename U>
	template< template<typename, size_t> class align>
	void Primitives<T,U>::allocate(const Doublet<Shell> &doublet) {
	    memory16_ = NULL;
	    C[0] = C[1] = NULL;

	    int N = doublet.L()/2 + 1;

	    const Shell &a = doublet[0], &b = doublet[1];

	    size_t NU = N + align<U,0>::get(N);
	    size2d = NU*(a.L + 1)*(b.L + 1);

	    size_t NT = N + align<T,0>::get(N);

	    //int nc = a->nc*b->nc;  --- for hybrid sp shells  a->nc => 1:non-hybrid or 2:hybrid
	    //int K = a->K*b->K;  --- a->K : number of primitives in shell a

	    int nc = doublet.nc();
	    int K = doublet.K();  //should always be 1 for fully un-contracted RI basis

	    Kmax = K;//(RYSQ_WITH_INT2D_BLOCK*1024)/((3*size2d+nc)*sizeof(double));
	    Kmax = std::max(Kmax, 1);
	    Kmax = std::min(Kmax, K);

	    size_t size = 0;
	    size += (3*size2d)*Kmax;
	    size_t bytes = size*sizeof(U) + nc*Kmax*sizeof(double);

	    int status = posix_memalign(&memory16_, 16, bytes);

	    assert(status == 0);
	    if (status != 0) throw std::bad_alloc();

	    U *ptr = (U*)memory16_;

	    T *ptr64 = (T*)ptr;
	    Ints2d = ptr64;
	    ptr64 += 3*size2d*Kmax;

	    for (int i = 0; i < (b.nc); ++i) {
		C[i] = (double*)ptr64 + i*(a.nc)*Kmax;
	    }
	}//Primitives<T,U>::allocate


    } // namespace twoc_quadrature





















	namespace threec_quadrature {


	    template<int N>
	    struct roots {
		static const size_t value = N;
	    };


	    template<class braket>
	    void apply(const Triplet<Shell> &triplet,
		       const vector<3> &ri, const vector<3> &rj,
		       const vector<3> &rk,
		       double (&Q)[braket::size], int transpose,
		       double cutoff);

	    template<typename T, typename U>
	    struct Primitives {
		int Kmax, K;
		int size2d;
		void *memory16_;
		U *Gx_, *Gy_, *Gz_;
		T *Ints2d;
		double *C[2];  //coefficients

		template<template<typename, size_t> class align>
		void allocate(const Triplet<Shell> &triplet);

		~Primitives() {	if (memory16_) free(memory16_); }

		U* Gx(int K) { return (Gx_) ? Gx_ : get<U,0>(K); }
		U* Gy(int K) { return (Gy_) ? Gy_ : get<U,1>(K); }
		U* Gz(int K) { return (Gz_) ? Gz_ : get<U,2>(K); }



		T* Ix(int K = 0) { return get<T,0>(K); }
		T* Iy(int K = 0) { return get<T,1>(K); }
		T* Iz(int K = 0) { return get<T,2>(K); }

		const T* Ix(int K = 0) const { return get<T,0>(K); }
		const T* Iy(int K = 0) const { return get<T,1>(K); }
		const T* Iz(int K = 0) const { return get<T,2>(K); }

	    private:
		template<typename R, size_t q> R* get(int i) {
		    return reinterpret_cast<R*>(&Ints2d[(3*i + q)*size2d]);
		}
		template<typename R, size_t q>  const R* get(int i)const { //q = 0,1,2  i = primitive
		    return (const_cast<Primitives*>(this))->get<R,q>(i);
		}
	    };


	    template<class bra, size_t N, template<typename, size_t> class align,
		     class Transform, typename T, typename U>
	    void apply(const Triplet<Shell> &triplet,
		       const vector<3> &ri, const vector<3> &rj,
		       const vector<3> &rk,
		       double scale, double cutoff,
		       Primitives<T,U> &primitives,
		       Transform &transform);

    template<typename T, typename U>
    template< template<typename, size_t> class align>
    void Primitives<T,U>::allocate(const Triplet<Shell> &triplet) {
	memory16_ = NULL;
	Gx_ = Gy_ = Gz_ = NULL;
	C[0] = C[1] = NULL;

	int N = triplet.L()/2 + 1;

	const Shell &a = triplet[0], &b = triplet[1], &c = triplet[2];

	int transfer = (a.L && b.L);

	size_t NU = N + align<U,0>::get(N);
	int size2d0 = NU*(a.L + b.L + 1)*(c.L + 1);


	size_t NT = N + align<T,0>::get(N);
	size2d = NT*(a.L + 1)*(b.L + 1)*(c.L + 1);

	//int nc = a->nc*b->nc*c->nc;  --- for hybrid sp shells  a->nc => 1:non-hybrid or 2:hybrid
	//int K = a->K*b->K*c->K;  --- a->K : number of primitives in shell a

	int nc = triplet.nc();  //nc=8,4,2,1
	int K = triplet.K();    //a->K*b->K*c->K  -> c->K is always "1"


	Kmax = K;//(RYSQ_WITH_INT2D_BLOCK*1024)/((3*size2d+nc)*sizeof(double));
	Kmax = std::max(Kmax, 1);
	Kmax = std::min(Kmax, K);

	size_t size_tmp = 0;
	size_tmp += (transfer) ? 3*size2d0 : 0;  //space for un-transfered 2d integrals

	//allocate space for 2D integrals
	size_t size = 0;
	size += (3*size2d)*Kmax; //space for 2D integrals

	size_t bytes = size_tmp*sizeof(U) + size*sizeof(T) + nc*Kmax*sizeof(double);


	int status = posix_memalign(&memory16_, 16, bytes);

	assert(status == 0);
	if (status != 0) throw std::bad_alloc();

	U *ptr = (U*)memory16_;
	if (transfer) {
	    Gx_ = ptr + 0*size2d0;
	    Gy_ = ptr + 1*size2d0;
	    Gz_ = ptr + 2*size2d0;
	    ptr += 3*size2d0;        //un-transfered 2d integrals
	}

	T *ptr64 = (T*)ptr;
	Ints2d = ptr64;
	ptr64 += 3*size2d*Kmax;  //transfered 2d integrals integrals
	//
	for (int i = 0; i < (c.nc); ++i) {
	    C[i] = (double*)ptr64 + i*(a.nc*b.nc)*Kmax; //finished 6d integrals
	} //i
    }//void Primitives<T,U>::allocate

} // namespace threec_quadrature





    } // namespace kernel
} // namespace rysq


#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_TWOC_QUADRATURE_HPP_ */
