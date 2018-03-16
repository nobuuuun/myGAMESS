/**  
 @file 
 @warning Automatically Generated
*/
/**  
 @warning AUTOMATICALLY GENERATED
*/



#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_TWOC_QUADRATURE_IMPL_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_DERIVATIVE_TWOC_QUADRATURE_IMPL_HPP_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include "meta.hpp"
#include "kernel/forward.hpp"


namespace rysq {
namespace kernel {
namespace twoc_derivative_quadrature {


#define RI_RYSQ_NORMAL_000 1.0000000000000000 /**< @brief (0 0 0) normalization constant */ 
#define RI_RYSQ_NORMAL_100 1.0000000000000000 /**< @brief (1 0 0) normalization constant */ 
#define RI_RYSQ_NORMAL_200 1.0000000000000000 /**< @brief (2 0 0) normalization constant */ 
#define RI_RYSQ_NORMAL_110 1.7320508075688772 /**< @brief (1 1 0) normalization constant */ 
#define RI_RYSQ_NORMAL_300 1.0000000000000000 /**< @brief (3 0 0) normalization constant */ 
#define RI_RYSQ_NORMAL_210 2.2360679774997898 /**< @brief (2 1 0) normalization constant */ 
#define RI_RYSQ_NORMAL_111 3.8729833462074170 /**< @brief (1 1 1) normalization constant */ 
#define RI_RYSQ_NORMAL_400 1.0000000000000000 /**< @brief (4 0 0) normalization constant */ 
#define RI_RYSQ_NORMAL_310 2.6457513110645907 /**< @brief (3 1 0) normalization constant */ 
#define RI_RYSQ_NORMAL_220 3.4156502553198660 /**< @brief (2 2 0) normalization constant */ 
#define RI_RYSQ_NORMAL_211 5.9160797830996161 /**< @brief (2 1 1) normalization constant */ 
#define RI_RYSQ_NORMAL_500 1.0000000000000000 /**< @brief (5 0 0) normalization constant */ 
#define RI_RYSQ_NORMAL_410 3.0000000000000000 /**< @brief (4 1 0) normalization constant */ 
#define RI_RYSQ_NORMAL_320 4.5825756949558398 /**< @brief (3 2 0) normalization constant */ 
#define RI_RYSQ_NORMAL_311 7.9372539331937721 /**< @brief (3 1 1) normalization constant */ 
#define RI_RYSQ_NORMAL_221 10.2469507659595980 /**< @brief (2 2 1) normalization constant */ 



const double NORMALIZE[] = {
	RI_RYSQ_NORMAL_000,
	RI_RYSQ_NORMAL_100, RI_RYSQ_NORMAL_100, RI_RYSQ_NORMAL_100,
	RI_RYSQ_NORMAL_200, RI_RYSQ_NORMAL_200, RI_RYSQ_NORMAL_200,
	RI_RYSQ_NORMAL_110, RI_RYSQ_NORMAL_110, RI_RYSQ_NORMAL_110,
	RI_RYSQ_NORMAL_300, RI_RYSQ_NORMAL_300, RI_RYSQ_NORMAL_300,
	RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210, RI_RYSQ_NORMAL_210,
	RI_RYSQ_NORMAL_111,
	RI_RYSQ_NORMAL_400, RI_RYSQ_NORMAL_400, RI_RYSQ_NORMAL_400,
	RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310, RI_RYSQ_NORMAL_310,
	RI_RYSQ_NORMAL_220, RI_RYSQ_NORMAL_220, RI_RYSQ_NORMAL_220,
	RI_RYSQ_NORMAL_211, RI_RYSQ_NORMAL_211, RI_RYSQ_NORMAL_211,
	RI_RYSQ_NORMAL_500, RI_RYSQ_NORMAL_500, RI_RYSQ_NORMAL_500,
	RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410, RI_RYSQ_NORMAL_410,
	RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320, RI_RYSQ_NORMAL_320,
	RI_RYSQ_NORMAL_311, RI_RYSQ_NORMAL_311, RI_RYSQ_NORMAL_311,
	RI_RYSQ_NORMAL_221, RI_RYSQ_NORMAL_221, RI_RYSQ_NORMAL_221};




//unrolled bras
#define Ix(a,i) (Ix[(a) + (i)*NT])
#define Iy(a,i) (Iy[(a) + (i)*NT])
#define Iz(a,i) (Iz[(a) + (i)*NT])


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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::S> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 1;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


        // function registers
	T q0 = 0.0;
	T q1 = 0.0;
	T q2 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q0 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[0];
	I[1] += q1*C_[0]*NORMALIZE[0];
	I[2] += q2*C_[0]*NORMALIZE[0];

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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::P> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q0 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q3 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q4 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q5 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q6 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q7 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q8 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[1];
	I[1] += q1*C_[0]*NORMALIZE[2];
	I[2] += q2*C_[0]*NORMALIZE[3];
	I[3] += q3*C_[0]*NORMALIZE[1];
	I[4] += q4*C_[0]*NORMALIZE[2];
	I[5] += q5*C_[0]*NORMALIZE[3];
	I[6] += q6*C_[0]*NORMALIZE[1];
	I[7] += q7*C_[0]*NORMALIZE[2];
	I[8] += q8*C_[0]*NORMALIZE[3];

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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::D> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 3;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q0 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,0) - 2*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	    q3 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,0) - 1*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q4 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,1) - 1*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q5 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q6 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q7 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,0) - 2*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q8 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q9 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,0) - 1*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[4];
	I[1] += q1*C_[0]*NORMALIZE[5];
	I[2] += q2*C_[0]*NORMALIZE[6];
	I[3] += q3*C_[0]*NORMALIZE[7];
	I[4] += q4*C_[0]*NORMALIZE[8];
	I[5] += q5*C_[0]*NORMALIZE[9];
	I[6] += q6*C_[0]*NORMALIZE[4];
	I[7] += q7*C_[0]*NORMALIZE[5];
	I[8] += q8*C_[0]*NORMALIZE[6];
	I[9] += q9*C_[0]*NORMALIZE[7];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q10 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q11 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,1) - 1*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q12 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q13 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q14 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,3) - 2*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q15 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q16 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,2) - 1*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q17 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,2) - 1*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[10] += q10*C_[0]*NORMALIZE[8];
	I[11] += q11*C_[0]*NORMALIZE[9];
	I[12] += q12*C_[0]*NORMALIZE[4];
	I[13] += q13*C_[0]*NORMALIZE[5];
	I[14] += q14*C_[0]*NORMALIZE[6];
	I[15] += q15*C_[0]*NORMALIZE[7];
	I[16] += q16*C_[0]*NORMALIZE[8];
	I[17] += q17*C_[0]*NORMALIZE[9];

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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::F> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 4;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q0 += 2*alpha*Ix(a,4)*Iy(a,0)*Iz(a,0) - 3*Ix(a,2)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,3); 
	    q3 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,0) - 2*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q4 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,1) - 2*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q5 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,0) - 1*Ix(a,0)*Iy(a,2)*Iz(a,0); 
	    q6 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	    q7 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,2); 
	    q8 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	    q9 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,1) - 1*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[10];
	I[1] += q1*C_[0]*NORMALIZE[11];
	I[2] += q2*C_[0]*NORMALIZE[12];
	I[3] += q3*C_[0]*NORMALIZE[13];
	I[4] += q4*C_[0]*NORMALIZE[14];
	I[5] += q5*C_[0]*NORMALIZE[15];
	I[6] += q6*C_[0]*NORMALIZE[16];
	I[7] += q7*C_[0]*NORMALIZE[17];
	I[8] += q8*C_[0]*NORMALIZE[18];
	I[9] += q9*C_[0]*NORMALIZE[19];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q10 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,0); 
	    q11 += 2*alpha*Ix(a,0)*Iy(a,4)*Iz(a,0) - 3*Ix(a,0)*Iy(a,2)*Iz(a,0); 
	    q12 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,3); 
	    q13 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,0) - 1*Ix(a,2)*Iy(a,0)*Iz(a,0); 
	    q14 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	    q15 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,0) - 2*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q16 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,1) - 2*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q17 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	    q18 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,2); 
	    q19 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,1) - 1*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[10] += q10*C_[0]*NORMALIZE[10];
	I[11] += q11*C_[0]*NORMALIZE[11];
	I[12] += q12*C_[0]*NORMALIZE[12];
	I[13] += q13*C_[0]*NORMALIZE[13];
	I[14] += q14*C_[0]*NORMALIZE[14];
	I[15] += q15*C_[0]*NORMALIZE[15];
	I[16] += q16*C_[0]*NORMALIZE[16];
	I[17] += q17*C_[0]*NORMALIZE[17];
	I[18] += q18*C_[0]*NORMALIZE[18];
	I[19] += q19*C_[0]*NORMALIZE[19];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q20 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,1); 
	    q21 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,1); 
	    q22 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,4) - 3*Ix(a,0)*Iy(a,0)*Iz(a,2); 
	    q23 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	    q24 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,2) - 1*Ix(a,2)*Iy(a,0)*Iz(a,0); 
	    q25 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	    q26 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,2) - 1*Ix(a,0)*Iy(a,2)*Iz(a,0); 
	    q27 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,3) - 2*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q28 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,3) - 2*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q29 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,2) - 1*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[20] += q20*C_[0]*NORMALIZE[10];
	I[21] += q21*C_[0]*NORMALIZE[11];
	I[22] += q22*C_[0]*NORMALIZE[12];
	I[23] += q23*C_[0]*NORMALIZE[13];
	I[24] += q24*C_[0]*NORMALIZE[14];
	I[25] += q25*C_[0]*NORMALIZE[15];
	I[26] += q26*C_[0]*NORMALIZE[16];
	I[27] += q27*C_[0]*NORMALIZE[17];
	I[28] += q28*C_[0]*NORMALIZE[18];
	I[29] += q29*C_[0]*NORMALIZE[19];

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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::G> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 5;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q0 += 2*alpha*Ix(a,5)*Iy(a,0)*Iz(a,0) - 4*Ix(a,3)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,1)*Iy(a,4)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,4); 
	    q3 += 2*alpha*Ix(a,4)*Iy(a,1)*Iz(a,0) - 3*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q4 += 2*alpha*Ix(a,4)*Iy(a,0)*Iz(a,1) - 3*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q5 += 2*alpha*Ix(a,2)*Iy(a,3)*Iz(a,0) - 1*Ix(a,0)*Iy(a,3)*Iz(a,0); 
	    q6 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,1); 
	    q7 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,3) - 1*Ix(a,0)*Iy(a,0)*Iz(a,3); 
	    q8 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,3); 
	    q9 += 2*alpha*Ix(a,3)*Iy(a,2)*Iz(a,0) - 2*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[20];
	I[1] += q1*C_[0]*NORMALIZE[21];
	I[2] += q2*C_[0]*NORMALIZE[22];
	I[3] += q3*C_[0]*NORMALIZE[23];
	I[4] += q4*C_[0]*NORMALIZE[24];
	I[5] += q5*C_[0]*NORMALIZE[25];
	I[6] += q6*C_[0]*NORMALIZE[26];
	I[7] += q7*C_[0]*NORMALIZE[27];
	I[8] += q8*C_[0]*NORMALIZE[28];
	I[9] += q9*C_[0]*NORMALIZE[29];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q10 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,2) - 2*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	    q11 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,2); 
	    q12 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,1) - 2*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q13 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,1) - 1*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q14 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,2) - 1*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q15 += 2*alpha*Ix(a,4)*Iy(a,1)*Iz(a,0); 
	    q16 += 2*alpha*Ix(a,0)*Iy(a,5)*Iz(a,0) - 4*Ix(a,0)*Iy(a,3)*Iz(a,0); 
	    q17 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,4); 
	    q18 += 2*alpha*Ix(a,3)*Iy(a,2)*Iz(a,0) - 1*Ix(a,3)*Iy(a,0)*Iz(a,0); 
	    q19 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[10] += q10*C_[0]*NORMALIZE[30];
	I[11] += q11*C_[0]*NORMALIZE[31];
	I[12] += q12*C_[0]*NORMALIZE[32];
	I[13] += q13*C_[0]*NORMALIZE[33];
	I[14] += q14*C_[0]*NORMALIZE[34];
	I[15] += q15*C_[0]*NORMALIZE[20];
	I[16] += q16*C_[0]*NORMALIZE[21];
	I[17] += q17*C_[0]*NORMALIZE[22];
	I[18] += q18*C_[0]*NORMALIZE[23];
	I[19] += q19*C_[0]*NORMALIZE[24];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q20 += 2*alpha*Ix(a,1)*Iy(a,4)*Iz(a,0) - 3*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	    q21 += 2*alpha*Ix(a,0)*Iy(a,4)*Iz(a,1) - 3*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q22 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,3); 
	    q23 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,3) - 1*Ix(a,0)*Iy(a,0)*Iz(a,3); 
	    q24 += 2*alpha*Ix(a,2)*Iy(a,3)*Iz(a,0) - 2*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q25 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,2); 
	    q26 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,2) - 2*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q27 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,1) - 1*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q28 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,1) - 2*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	    q29 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,2) - 1*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	}//a
	    
	//contraction coefficients and normalize
	I[20] += q20*C_[0]*NORMALIZE[25];
	I[21] += q21*C_[0]*NORMALIZE[26];
	I[22] += q22*C_[0]*NORMALIZE[27];
	I[23] += q23*C_[0]*NORMALIZE[28];
	I[24] += q24*C_[0]*NORMALIZE[29];
	I[25] += q25*C_[0]*NORMALIZE[30];
	I[26] += q26*C_[0]*NORMALIZE[31];
	I[27] += q27*C_[0]*NORMALIZE[32];
	I[28] += q28*C_[0]*NORMALIZE[33];
	I[29] += q29*C_[0]*NORMALIZE[34];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q30 += 2*alpha*Ix(a,4)*Iy(a,0)*Iz(a,1); 
	    q31 += 2*alpha*Ix(a,0)*Iy(a,4)*Iz(a,1); 
	    q32 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,5) - 4*Ix(a,0)*Iy(a,0)*Iz(a,3); 
	    q33 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,1); 
	    q34 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,2) - 1*Ix(a,3)*Iy(a,0)*Iz(a,0); 
	    q35 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,1); 
	    q36 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,2) - 1*Ix(a,0)*Iy(a,3)*Iz(a,0); 
	    q37 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,4) - 3*Ix(a,1)*Iy(a,0)*Iz(a,2); 
	    q38 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,4) - 3*Ix(a,0)*Iy(a,1)*Iz(a,2); 
	    q39 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[30] += q30*C_[0]*NORMALIZE[20];
	I[31] += q31*C_[0]*NORMALIZE[21];
	I[32] += q32*C_[0]*NORMALIZE[22];
	I[33] += q33*C_[0]*NORMALIZE[23];
	I[34] += q34*C_[0]*NORMALIZE[24];
	I[35] += q35*C_[0]*NORMALIZE[25];
	I[36] += q36*C_[0]*NORMALIZE[26];
	I[37] += q37*C_[0]*NORMALIZE[27];
	I[38] += q38*C_[0]*NORMALIZE[28];
	I[39] += q39*C_[0]*NORMALIZE[29];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q40 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,3) - 2*Ix(a,2)*Iy(a,0)*Iz(a,1); 
	    q41 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,3) - 2*Ix(a,0)*Iy(a,2)*Iz(a,1); 
	    q42 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,2) - 1*Ix(a,2)*Iy(a,1)*Iz(a,0); 
	    q43 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,2) - 1*Ix(a,1)*Iy(a,2)*Iz(a,0); 
	    q44 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,3) - 2*Ix(a,1)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[40] += q40*C_[0]*NORMALIZE[30];
	I[41] += q41*C_[0]*NORMALIZE[31];
	I[42] += q42*C_[0]*NORMALIZE[32];
	I[43] += q43*C_[0]*NORMALIZE[33];
	I[44] += q44*C_[0]*NORMALIZE[34];

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
    @brief <h| shell quadrature
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
struct impl< meta::OneCenterState<rysq::H> > {

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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::H> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 6;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q0 += 2*alpha*Ix(a,6)*Iy(a,0)*Iz(a,0) - 5*Ix(a,4)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,1)*Iy(a,5)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,5); 
	    q3 += 2*alpha*Ix(a,5)*Iy(a,1)*Iz(a,0) - 4*Ix(a,3)*Iy(a,1)*Iz(a,0); 
	    q4 += 2*alpha*Ix(a,5)*Iy(a,0)*Iz(a,1) - 4*Ix(a,3)*Iy(a,0)*Iz(a,1); 
	    q5 += 2*alpha*Ix(a,2)*Iy(a,4)*Iz(a,0) - 1*Ix(a,0)*Iy(a,4)*Iz(a,0); 
	    q6 += 2*alpha*Ix(a,1)*Iy(a,4)*Iz(a,1); 
	    q7 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,4) - 1*Ix(a,0)*Iy(a,0)*Iz(a,4); 
	    q8 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,4); 
	    q9 += 2*alpha*Ix(a,4)*Iy(a,2)*Iz(a,0) - 3*Ix(a,2)*Iy(a,2)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[35];
	I[1] += q1*C_[0]*NORMALIZE[36];
	I[2] += q2*C_[0]*NORMALIZE[37];
	I[3] += q3*C_[0]*NORMALIZE[38];
	I[4] += q4*C_[0]*NORMALIZE[39];
	I[5] += q5*C_[0]*NORMALIZE[40];
	I[6] += q6*C_[0]*NORMALIZE[41];
	I[7] += q7*C_[0]*NORMALIZE[42];
	I[8] += q8*C_[0]*NORMALIZE[43];
	I[9] += q9*C_[0]*NORMALIZE[44];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q10 += 2*alpha*Ix(a,4)*Iy(a,0)*Iz(a,2) - 3*Ix(a,2)*Iy(a,0)*Iz(a,2); 
	    q11 += 2*alpha*Ix(a,3)*Iy(a,3)*Iz(a,0) - 2*Ix(a,1)*Iy(a,3)*Iz(a,0); 
	    q12 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,2); 
	    q13 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,3) - 2*Ix(a,1)*Iy(a,0)*Iz(a,3); 
	    q14 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,3); 
	    q15 += 2*alpha*Ix(a,4)*Iy(a,1)*Iz(a,1) - 3*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	    q16 += 2*alpha*Ix(a,2)*Iy(a,3)*Iz(a,1) - 1*Ix(a,0)*Iy(a,3)*Iz(a,1); 
	    q17 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,3) - 1*Ix(a,0)*Iy(a,1)*Iz(a,3); 
	    q18 += 2*alpha*Ix(a,3)*Iy(a,2)*Iz(a,1) - 2*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	    q19 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,2) - 2*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	}//a
	    
	//contraction coefficients and normalize
	I[10] += q10*C_[0]*NORMALIZE[45];
	I[11] += q11*C_[0]*NORMALIZE[46];
	I[12] += q12*C_[0]*NORMALIZE[47];
	I[13] += q13*C_[0]*NORMALIZE[48];
	I[14] += q14*C_[0]*NORMALIZE[49];
	I[15] += q15*C_[0]*NORMALIZE[50];
	I[16] += q16*C_[0]*NORMALIZE[51];
	I[17] += q17*C_[0]*NORMALIZE[52];
	I[18] += q18*C_[0]*NORMALIZE[53];
	I[19] += q19*C_[0]*NORMALIZE[54];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q20 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,2) - 1*Ix(a,0)*Iy(a,2)*Iz(a,2); 
	    q21 += 2*alpha*Ix(a,5)*Iy(a,1)*Iz(a,0); 
	    q22 += 2*alpha*Ix(a,0)*Iy(a,6)*Iz(a,0) - 5*Ix(a,0)*Iy(a,4)*Iz(a,0); 
	    q23 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,5); 
	    q24 += 2*alpha*Ix(a,4)*Iy(a,2)*Iz(a,0) - 1*Ix(a,4)*Iy(a,0)*Iz(a,0); 
	    q25 += 2*alpha*Ix(a,4)*Iy(a,1)*Iz(a,1); 
	    q26 += 2*alpha*Ix(a,1)*Iy(a,5)*Iz(a,0) - 4*Ix(a,1)*Iy(a,3)*Iz(a,0); 
	    q27 += 2*alpha*Ix(a,0)*Iy(a,5)*Iz(a,1) - 4*Ix(a,0)*Iy(a,3)*Iz(a,1); 
	    q28 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,4); 
	    q29 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,4) - 1*Ix(a,0)*Iy(a,0)*Iz(a,4); 
	}//a
	    
	//contraction coefficients and normalize
	I[20] += q20*C_[0]*NORMALIZE[55];
	I[21] += q21*C_[0]*NORMALIZE[35];
	I[22] += q22*C_[0]*NORMALIZE[36];
	I[23] += q23*C_[0]*NORMALIZE[37];
	I[24] += q24*C_[0]*NORMALIZE[38];
	I[25] += q25*C_[0]*NORMALIZE[39];
	I[26] += q26*C_[0]*NORMALIZE[40];
	I[27] += q27*C_[0]*NORMALIZE[41];
	I[28] += q28*C_[0]*NORMALIZE[42];
	I[29] += q29*C_[0]*NORMALIZE[43];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


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
	    q30 += 2*alpha*Ix(a,3)*Iy(a,3)*Iz(a,0) - 2*Ix(a,3)*Iy(a,1)*Iz(a,0); 
	    q31 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,2); 
	    q32 += 2*alpha*Ix(a,2)*Iy(a,4)*Iz(a,0) - 3*Ix(a,2)*Iy(a,2)*Iz(a,0); 
	    q33 += 2*alpha*Ix(a,0)*Iy(a,4)*Iz(a,2) - 3*Ix(a,0)*Iy(a,2)*Iz(a,2); 
	    q34 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,3); 
	    q35 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,3) - 2*Ix(a,0)*Iy(a,1)*Iz(a,3); 
	    q36 += 2*alpha*Ix(a,3)*Iy(a,2)*Iz(a,1) - 1*Ix(a,3)*Iy(a,0)*Iz(a,1); 
	    q37 += 2*alpha*Ix(a,1)*Iy(a,4)*Iz(a,1) - 3*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	    q38 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,3) - 1*Ix(a,1)*Iy(a,0)*Iz(a,3); 
	    q39 += 2*alpha*Ix(a,2)*Iy(a,3)*Iz(a,1) - 2*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[30] += q30*C_[0]*NORMALIZE[44];
	I[31] += q31*C_[0]*NORMALIZE[45];
	I[32] += q32*C_[0]*NORMALIZE[46];
	I[33] += q33*C_[0]*NORMALIZE[47];
	I[34] += q34*C_[0]*NORMALIZE[48];
	I[35] += q35*C_[0]*NORMALIZE[49];
	I[36] += q36*C_[0]*NORMALIZE[50];
	I[37] += q37*C_[0]*NORMALIZE[51];
	I[38] += q38*C_[0]*NORMALIZE[52];
	I[39] += q39*C_[0]*NORMALIZE[53];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


        // function registers
	T q40 = 0.0;
	T q41 = 0.0;
	T q42 = 0.0;
	T q43 = 0.0;
	T q44 = 0.0;
	T q45 = 0.0;
	T q46 = 0.0;
	T q47 = 0.0;
	T q48 = 0.0;
	T q49 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q40 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,2) - 1*Ix(a,2)*Iy(a,0)*Iz(a,2); 
	    q41 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,2) - 2*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	    q42 += 2*alpha*Ix(a,5)*Iy(a,0)*Iz(a,1); 
	    q43 += 2*alpha*Ix(a,0)*Iy(a,5)*Iz(a,1); 
	    q44 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,6) - 5*Ix(a,0)*Iy(a,0)*Iz(a,4); 
	    q45 += 2*alpha*Ix(a,4)*Iy(a,1)*Iz(a,1); 
	    q46 += 2*alpha*Ix(a,4)*Iy(a,0)*Iz(a,2) - 1*Ix(a,4)*Iy(a,0)*Iz(a,0); 
	    q47 += 2*alpha*Ix(a,1)*Iy(a,4)*Iz(a,1); 
	    q48 += 2*alpha*Ix(a,0)*Iy(a,4)*Iz(a,2) - 1*Ix(a,0)*Iy(a,4)*Iz(a,0); 
	    q49 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,5) - 4*Ix(a,1)*Iy(a,0)*Iz(a,3); 
	}//a
	    
	//contraction coefficients and normalize
	I[40] += q40*C_[0]*NORMALIZE[54];
	I[41] += q41*C_[0]*NORMALIZE[55];
	I[42] += q42*C_[0]*NORMALIZE[35];
	I[43] += q43*C_[0]*NORMALIZE[36];
	I[44] += q44*C_[0]*NORMALIZE[37];
	I[45] += q45*C_[0]*NORMALIZE[38];
	I[46] += q46*C_[0]*NORMALIZE[39];
	I[47] += q47*C_[0]*NORMALIZE[40];
	I[48] += q48*C_[0]*NORMALIZE[41];
	I[49] += q49*C_[0]*NORMALIZE[42];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


        // function registers
	T q50 = 0.0;
	T q51 = 0.0;
	T q52 = 0.0;
	T q53 = 0.0;
	T q54 = 0.0;
	T q55 = 0.0;
	T q56 = 0.0;
	T q57 = 0.0;
	T q58 = 0.0;
	T q59 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q50 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,5) - 4*Ix(a,0)*Iy(a,1)*Iz(a,3); 
	    q51 += 2*alpha*Ix(a,3)*Iy(a,2)*Iz(a,1); 
	    q52 += 2*alpha*Ix(a,3)*Iy(a,0)*Iz(a,3) - 2*Ix(a,3)*Iy(a,0)*Iz(a,1); 
	    q53 += 2*alpha*Ix(a,2)*Iy(a,3)*Iz(a,1); 
	    q54 += 2*alpha*Ix(a,0)*Iy(a,3)*Iz(a,3) - 2*Ix(a,0)*Iy(a,3)*Iz(a,1); 
	    q55 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,4) - 3*Ix(a,2)*Iy(a,0)*Iz(a,2); 
	    q56 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,4) - 3*Ix(a,0)*Iy(a,2)*Iz(a,2); 
	    q57 += 2*alpha*Ix(a,3)*Iy(a,1)*Iz(a,2) - 1*Ix(a,3)*Iy(a,1)*Iz(a,0); 
	    q58 += 2*alpha*Ix(a,1)*Iy(a,3)*Iz(a,2) - 1*Ix(a,1)*Iy(a,3)*Iz(a,0); 
	    q59 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,4) - 3*Ix(a,1)*Iy(a,1)*Iz(a,2); 
	}//a
	    
	//contraction coefficients and normalize
	I[50] += q50*C_[0]*NORMALIZE[43];
	I[51] += q51*C_[0]*NORMALIZE[44];
	I[52] += q52*C_[0]*NORMALIZE[45];
	I[53] += q53*C_[0]*NORMALIZE[46];
	I[54] += q54*C_[0]*NORMALIZE[47];
	I[55] += q55*C_[0]*NORMALIZE[48];
	I[56] += q56*C_[0]*NORMALIZE[49];
	I[57] += q57*C_[0]*NORMALIZE[50];
	I[58] += q58*C_[0]*NORMALIZE[51];
	I[59] += q59*C_[0]*NORMALIZE[52];

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

    for(int k = 0; k < K*1; k += 1) {
	double C_[1];
	for (int i = 0; i < 1; ++ i) C_[i] = C[k + i];


        // function registers
	T q60 = 0.0;
	T q61 = 0.0;
	T q62 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q60 += 2*alpha*Ix(a,2)*Iy(a,2)*Iz(a,2) - 1*Ix(a,2)*Iy(a,2)*Iz(a,0); 
	    q61 += 2*alpha*Ix(a,2)*Iy(a,1)*Iz(a,3) - 2*Ix(a,2)*Iy(a,1)*Iz(a,1); 
	    q62 += 2*alpha*Ix(a,1)*Iy(a,2)*Iz(a,3) - 2*Ix(a,1)*Iy(a,2)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[60] += q60*C_[0]*NORMALIZE[53];
	I[61] += q61*C_[0]*NORMALIZE[54];
	I[62] += q62*C_[0]*NORMALIZE[55];

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
			int dim2d,
			const typename aligned<T>::type *__restrict Ix,
			const typename aligned<T>::type *__restrict Iy, 
			const typename aligned<T>::type *__restrict Iz, 
			double scale,
                        const double alpha,
			double *__restrict I);// __attribute__((pure));
};

template<size_t N, typename T, size_t NT>
size_t impl< meta::OneCenterState<rysq::SP> >::apply(bool normalize,
			     double tol, int K,
			     const double *__restrict C,
			     int dim2d,
			     const typename aligned<T>::type *__restrict Ix,
			     const typename aligned<T>::type *__restrict Iy, 
			     const typename aligned<T>::type *__restrict Iz, 
			     double scale,
			     const double alpha,
			     double *__restrict I) {

    const int Li1 = 2;
    
    int num = 0;
    
    
    
#ifdef __GNUG__
     asm("#begin vector loop");
#endif

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];


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
	    q0 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,0); 
	    q1 += 2*alpha*Ix(a,2)*Iy(a,0)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q2 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q3 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	    q4 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,0); 
	    q5 += 2*alpha*Ix(a,1)*Iy(a,1)*Iz(a,0); 
	    q6 += 2*alpha*Ix(a,0)*Iy(a,2)*Iz(a,0) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	    q7 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q8 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,1); 
	    q9 += 2*alpha*Ix(a,1)*Iy(a,0)*Iz(a,1); 
	}//a
	    
	//contraction coefficients and normalize
	I[0] += q0*C_[0]*NORMALIZE[1];
	I[1] += q1*C_[1]*NORMALIZE[2];
	I[2] += q2*C_[1]*NORMALIZE[3];
	I[3] += q3*C_[1]*NORMALIZE[4];
	I[4] += q4*C_[0]*NORMALIZE[1];
	I[5] += q5*C_[1]*NORMALIZE[2];
	I[6] += q6*C_[1]*NORMALIZE[3];
	I[7] += q7*C_[1]*NORMALIZE[4];
	I[8] += q8*C_[0]*NORMALIZE[1];
	I[9] += q9*C_[1]*NORMALIZE[2];

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

    for(int k = 0; k < K*2; k += 2) {
	double C_[2];
	for (int i = 0; i < 2; ++ i) C_[i] = C[k + i];


        // function registers
	T q10 = 0.0;
	T q11 = 0.0;
    
#if defined (__INTEL_COMPILER) 
#pragma vector aligned
#endif // alignment attribute

	for (int a = 0; a < int(N); ++a) {
	    q10 += 2*alpha*Ix(a,0)*Iy(a,1)*Iz(a,1); 
	    q11 += 2*alpha*Ix(a,0)*Iy(a,0)*Iz(a,2) - 1*Ix(a,0)*Iy(a,0)*Iz(a,0); 
	}//a
	    
	//contraction coefficients and normalize
	I[10] += q10*C_[1]*NORMALIZE[3];
	I[11] += q11*C_[1]*NORMALIZE[4];

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

