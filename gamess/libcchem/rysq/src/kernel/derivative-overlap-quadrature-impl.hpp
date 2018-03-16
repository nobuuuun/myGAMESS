/**  
 @file 
 @warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/



#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_OVERLAP_QUADRATURE_IMPL_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_OVERLAP_QUADRATURE_IMPL_HPP_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include "meta.hpp"
#include "kernel/ri-forward.hpp"


namespace rysq {
namespace kernel {
namespace twoc_derivative_overlap_quadrature {
//unrolled bras
#define Ix(a,i) (Ix[(a) + (i)])
#define Iy(a,i) (Iy[(a) + (i)])
#define Iz(a,i) (Iz[(a) + (i)])


/** 
    @brief <s| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param alpha exponent of bra gaussian
    @param[out] I integral batch
    @return number of screened integrals
*/

template<>
struct impl< meta::OneCenterState<rysq::S> > {

    template<typename T>
    struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
                        const double *__restrict AExp,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I, double *__restrict vfac);
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
                             const double *__restrict AExp,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I, double *__restrict vfac) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q1 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q2 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[0] += q0*C_[0];
        I[1] += q1*C_[0];
        I[2] += q2*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
    return num;
}//apply



/** 
    @brief <p| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param alpha exponent of bra gaussian
    @param[out] I integral batch
    @return number of screened integrals
*/

template<>
struct impl< meta::OneCenterState<rysq::P> > {

    template<typename T>
    struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
                        const double *__restrict AExp,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I, double *__restrict vfac);
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
                             const double *__restrict AExp,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I, double *__restrict vfac) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q1 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q2 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q3 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q4 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q5 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q6 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q7 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q8 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
        I[0] += q0*C_[0];
        I[1] += q1*C_[0];
        I[2] += q2*C_[0];
        I[3] += q3*C_[0];
        I[4] += q4*C_[0];
        I[5] += q5*C_[0];
        I[6] += q6*C_[0];
        I[7] += q7*C_[0];
        I[8] += q8*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
    return num;
}//apply



/** 
    @brief <d| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param alpha exponent of bra gaussian
    @param[out] I integral batch
    @return number of screened integrals
*/

template<>
struct impl< meta::OneCenterState<rysq::D> > {

    template<typename T>
    struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
                        const double *__restrict AExp,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I, double *__restrict vfac);
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
                             const double *__restrict AExp,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I, double *__restrict vfac) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += alpha*Ix(a,3)*Iy(a,0)*Iz(a,0) - 2*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q1 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	    q2 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	    q3 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,0) - 1*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q4 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,1) - 1*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q5 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q6 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q7 += alpha*Ix(a,0)*Iy(a,3)*Iz(a,0) - 2*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q8 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q9 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,0) - 1*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
        I[0] += q0*C_[0];
        I[1] += q1*C_[0];
        I[2] += q2*C_[0];
        I[3] += q3*C_[0];
        I[4] += q4*C_[0];
        I[5] += q5*C_[0];
        I[6] += q6*C_[0];
        I[7] += q7*C_[0];
        I[8] += q8*C_[0];
        I[9] += q9*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q11 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,1) - 1*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q12 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q13 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q14 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,3) - 2*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q15 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q16 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,2) - 1*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q17 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,2) - 1*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
        I[10] += q10*C_[0];
        I[11] += q11*C_[0];
        I[12] += q12*C_[0];
        I[13] += q13*C_[0];
        I[14] += q14*C_[0];
        I[15] += q15*C_[0];
        I[16] += q16*C_[0];
        I[17] += q17*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
    return num;
}//apply



/** 
    @brief <f| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param alpha exponent of bra gaussian
    @param[out] I integral batch
    @return number of screened integrals
*/

template<>
struct impl< meta::OneCenterState<rysq::F> > {

    template<typename T>
    struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
                        const double *__restrict AExp,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I, double *__restrict vfac);
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
                             const double *__restrict AExp,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I, double *__restrict vfac) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += alpha*Ix(a,4)*Iy(a,0)*Iz(a,0) - 3*Ix(a,2)*Iy(a,0)*Iz(a,0); 
	    q1 += alpha*Ix(a,1)*Iy(a,3)*Iz(a,0); 
	    q2 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,3); 
	    q3 += alpha*Ix(a,3)*Iy(a,1)*Iz(a,0) - 2*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q4 += alpha*Ix(a,3)*Iy(a,0)*Iz(a,1) - 2*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q5 += alpha*Ix(a,2)*Iy(a,2)*Iz(a,0) - 1*Ix(a,0)*Iy(a,2)*Iz(a,0); 
	    q6 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	    q7 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,2); 
	    q8 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	    q9 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,1) - 1*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[0] += q0*C_[0];
        I[1] += q1*C_[0];
        I[2] += q2*C_[0];
        I[3] += q3*C_[0];
        I[4] += q4*C_[0];
        I[5] += q5*C_[0];
        I[6] += q6*C_[0];
        I[7] += q7*C_[0];
        I[8] += q8*C_[0];
        I[9] += q9*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += alpha*Ix(a,3)*Iy(a,1)*Iz(a,0); 
	    q11 += alpha*Ix(a,0)*Iy(a,4)*Iz(a,0) - 3*Ix(a,0)*Iy(a,2)*Iz(a,0); 
	    q12 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,3); 
	    q13 += alpha*Ix(a,2)*Iy(a,2)*Iz(a,0) - 1*Ix(a,2)*Iy(a,0)*Iz(a,0); 
	    q14 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	    q15 += alpha*Ix(a,1)*Iy(a,3)*Iz(a,0) - 2*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q16 += alpha*Ix(a,0)*Iy(a,3)*Iz(a,1) - 2*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q17 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	    q18 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,2); 
	    q19 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,1) - 1*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[10] += q10*C_[0];
        I[11] += q11*C_[0];
        I[12] += q12*C_[0];
        I[13] += q13*C_[0];
        I[14] += q14*C_[0];
        I[15] += q15*C_[0];
        I[16] += q16*C_[0];
        I[17] += q17*C_[0];
        I[18] += q18*C_[0];
        I[19] += q19*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += alpha*Ix(a,3)*Iy(a,0)*Iz(a,1); 
	    q21 += alpha*Ix(a,0)*Iy(a,3)*Iz(a,1); 
	    q22 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,4) - 3*Ix(a,0)*Iy(a,0)*Iz(a,2); 
	    q23 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	    q24 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,2) - 1*Ix(a,2)*Iy(a,0)*Iz(a,0); 
	    q25 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	    q26 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,2) - 1*Ix(a,0)*Iy(a,2)*Iz(a,0); 
	    q27 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,3) - 2*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q28 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,3) - 2*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q29 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,2) - 1*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
        I[20] += q20*C_[0];
        I[21] += q21*C_[0];
        I[22] += q22*C_[0];
        I[23] += q23*C_[0];
        I[24] += q24*C_[0];
        I[25] += q25*C_[0];
        I[26] += q26*C_[0];
        I[27] += q27*C_[0];
        I[28] += q28*C_[0];
        I[29] += q29*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
    return num;
}//apply



/** 
    @brief <g| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param alpha exponent of bra gaussian
    @param[out] I integral batch
    @return number of screened integrals
*/

template<>
struct impl< meta::OneCenterState<rysq::G> > {

    template<typename T>
    struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
                        const double *__restrict AExp,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I, double *__restrict vfac);
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::G> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
                             const double *__restrict AExp,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I, double *__restrict vfac) {

    const int Li1 = 5;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += alpha*Ix(a,5)*Iy(a,0)*Iz(a,0) - 4*Ix(a,3)*Iy(a,0)*Iz(a,0); 
	    q1 += alpha*Ix(a,1)*Iy(a,4)*Iz(a,0); 
	    q2 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,4); 
	    q3 += alpha*Ix(a,4)*Iy(a,1)*Iz(a,0) - 3*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q4 += alpha*Ix(a,4)*Iy(a,0)*Iz(a,1) - 3*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q5 += alpha*Ix(a,2)*Iy(a,3)*Iz(a,0) - 1*Ix(a,0)*Iy(a,3)*Iz(a,0); 
	    q6 += alpha*Ix(a,1)*Iy(a,3)*Iz(a,1); 
	    q7 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,3) - 1*Ix(a,0)*Iy(a,0)*Iz(a,3); 
	    q8 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,3); 
	    q9 += alpha*Ix(a,3)*Iy(a,2)*Iz(a,0) - 2*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
        I[0] += q0*C_[0];
        I[1] += q1*C_[0];
        I[2] += q2*C_[0];
        I[3] += q3*C_[0];
        I[4] += q4*C_[0];
        I[5] += q5*C_[0];
        I[6] += q6*C_[0];
        I[7] += q7*C_[0];
        I[8] += q8*C_[0];
        I[9] += q9*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q10 = 0.0;
	T q11 = 0.0;
	T q12 = 0.0;
	T q13 = 0.0;
	T q14 = 0.0;
	T q15 = 0.0;
	T q16 = 0.0;
	T q17 = 0.0;
	T q18 = 0.0;
	T q19 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += alpha*Ix(a,3)*Iy(a,0)*Iz(a,2) - 2*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	    q11 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,2); 
	    q12 += alpha*Ix(a,3)*Iy(a,1)*Iz(a,1) - 2*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q13 += alpha*Ix(a,2)*Iy(a,2)*Iz(a,1) - 1*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q14 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,2) - 1*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q15 += alpha*Ix(a,4)*Iy(a,1)*Iz(a,0); 
	    q16 += alpha*Ix(a,0)*Iy(a,5)*Iz(a,0) - 4*Ix(a,0)*Iy(a,3)*Iz(a,0); 
	    q17 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,4); 
	    q18 += alpha*Ix(a,3)*Iy(a,2)*Iz(a,0) - 1*Ix(a,3)*Iy(a,0)*Iz(a,0); 
	    q19 += alpha*Ix(a,3)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[10] += q10*C_[0];
        I[11] += q11*C_[0];
        I[12] += q12*C_[0];
        I[13] += q13*C_[0];
        I[14] += q14*C_[0];
        I[15] += q15*C_[0];
        I[16] += q16*C_[0];
        I[17] += q17*C_[0];
        I[18] += q18*C_[0];
        I[19] += q19*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q20 = 0.0;
	T q21 = 0.0;
	T q22 = 0.0;
	T q23 = 0.0;
	T q24 = 0.0;
	T q25 = 0.0;
	T q26 = 0.0;
	T q27 = 0.0;
	T q28 = 0.0;
	T q29 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q20 += alpha*Ix(a,1)*Iy(a,4)*Iz(a,0) - 3*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	    q21 += alpha*Ix(a,0)*Iy(a,4)*Iz(a,1) - 3*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q22 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,3); 
	    q23 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,3) - 1*Ix(a,0)*Iy(a,0)*Iz(a,3); 
	    q24 += alpha*Ix(a,2)*Iy(a,3)*Iz(a,0) - 2*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q25 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,2); 
	    q26 += alpha*Ix(a,0)*Iy(a,3)*Iz(a,2) - 2*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q27 += alpha*Ix(a,2)*Iy(a,2)*Iz(a,1) - 1*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q28 += alpha*Ix(a,1)*Iy(a,3)*Iz(a,1) - 2*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q29 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,2) - 1*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	}//a
	    
	//contraction coefficients and normalize
        I[20] += q20*C_[0];
        I[21] += q21*C_[0];
        I[22] += q22*C_[0];
        I[23] += q23*C_[0];
        I[24] += q24*C_[0];
        I[25] += q25*C_[0];
        I[26] += q26*C_[0];
        I[27] += q27*C_[0];
        I[28] += q28*C_[0];
        I[29] += q29*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q30 = 0.0;
	T q31 = 0.0;
	T q32 = 0.0;
	T q33 = 0.0;
	T q34 = 0.0;
	T q35 = 0.0;
	T q36 = 0.0;
	T q37 = 0.0;
	T q38 = 0.0;
	T q39 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q30 += alpha*Ix(a,4)*Iy(a,0)*Iz(a,1); 
	    q31 += alpha*Ix(a,0)*Iy(a,4)*Iz(a,1); 
	    q32 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,5) - 4*Ix(a,0)*Iy(a,0)*Iz(a,3); 
	    q33 += alpha*Ix(a,3)*Iy(a,1)*Iz(a,1); 
	    q34 += alpha*Ix(a,3)*Iy(a,0)*Iz(a,2) - 1*Ix(a,3)*Iy(a,0)*Iz(a,0); 
	    q35 += alpha*Ix(a,1)*Iy(a,3)*Iz(a,1); 
	    q36 += alpha*Ix(a,0)*Iy(a,3)*Iz(a,2) - 1*Ix(a,0)*Iy(a,3)*Iz(a,0); 
	    q37 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,4) - 3*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	    q38 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,4) - 3*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q39 += alpha*Ix(a,2)*Iy(a,2)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[30] += q30*C_[0];
        I[31] += q31*C_[0];
        I[32] += q32*C_[0];
        I[33] += q33*C_[0];
        I[34] += q34*C_[0];
        I[35] += q35*C_[0];
        I[36] += q36*C_[0];
        I[37] += q37*C_[0];
        I[38] += q38*C_[0];
        I[39] += q39*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[1];
        const double alpha = AExp[k];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k*1 + i]*vfac[k];


        // function registers
	T q40 = 0.0;
	T q41 = 0.0;
	T q42 = 0.0;
	T q43 = 0.0;
	T q44 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q40 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,3) - 2*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q41 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,3) - 2*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q42 += alpha*Ix(a,2)*Iy(a,1)*Iz(a,2) - 1*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q43 += alpha*Ix(a,1)*Iy(a,2)*Iz(a,2) - 1*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	    q44 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,3) - 2*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[40] += q40*C_[0];
        I[41] += q41*C_[0];
        I[42] += q42*C_[0];
        I[43] += q43*C_[0];
        I[44] += q44*C_[0];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
    return num;
}//apply



/** 
    @brief <sp| shell quadrature
    @param normalize
    @param tol tolerance
    @param K number of contractions
    @param C contraction coefficients
    @param dim2d 2-D integrals dimensions
    @param Ix 2-D integral, Ix(N,Li,Lj,K,Lk,Ll)
    @param Iy 2-D integral, Iy(N,Li,Lj,K,Lk,Ll)
    @param Iz 2-D integral, Iz(N,Li,Lj,K,Lk,Ll)
    @param scale scale factor
    @param alpha exponent of bra gaussian
    @param[out] I integral batch
    @return number of screened integrals
*/

template<>
struct impl< meta::OneCenterState<rysq::SP> > {

    template<typename T>
    struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	typedef T type;
#else
	typedef T type __attribute__((aligned (16)));
#endif
    };

    template<size_t N, typename T, size_t NT>
    static size_t apply(bool normalize,
			double tol, int K,
			const double *__restrict C,
                        const double *__restrict AExp,
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
			double *__restrict I, double *__restrict vfac);
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
                             const double *__restrict AExp,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     double *__restrict I, double *__restrict vfac) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[2];
        const double alpha = AExp[k];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k*2 + i]*vfac[k];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
	T q3 = 0.0;
	T q4 = 0.0;
	T q5 = 0.0;
	T q6 = 0.0;
	T q7 = 0.0;
	T q8 = 0.0;
	T q9 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q1 += alpha*Ix(a,2)*Iy(a,0)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q2 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q3 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q4 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q5 += alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q6 += alpha*Ix(a,0)*Iy(a,2)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q7 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q8 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q9 += alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
        I[0] += q0*C_[0];
        I[1] += q1*C_[1];
        I[2] += q2*C_[1];
        I[3] += q3*C_[1];
        I[4] += q4*C_[0];
        I[5] += q5*C_[1];
        I[6] += q6*C_[1];
        I[7] += q7*C_[1];
        I[8] += q8*C_[0];
        I[9] += q9*C_[1];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K; k += 1) {
	double C_[2];
        const double alpha = AExp[k];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k*2 + i]*vfac[k];


        // function registers
	T q10 = 0.0;
	T q11 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q11 += alpha*Ix(a,0)*Iy(a,0)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
        I[10] += q10*C_[1];
        I[11] += q11*C_[1];

	Ix += 3*dim2d;
	Iy += 3*dim2d;
	Iz += 3*dim2d;
	
    }//k
    Ix = Ix - 3*dim2d*K;
    Iy = Iy - 3*dim2d*K;
    Iz = Iz - 3*dim2d*K;
    
    
    return num;
}//apply





#undef Ix
#undef Iy
#undef Iz


} // namespace rysq
} // namespace kernel
} // namespace twoc_derivative_quadrature



#endif // LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_TWOC_QUADRATURE_IMPL_HPP_

