

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>

/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>




#define index(a,b,offset) (a + b*offset)

__global__ void invokeGradientPvv(const int nv,
				  const int iocc,
				  const int jocc,
				  const double *d_ea,
				  const double *d_ev,
				  //double *d_pAC,
				  //double *d_pBC,
				  double *d_pEri,
				  double *d_pEnergy,
				  double *d_pTijLong,
				  double *d_pEriLong)
{


    const double denom_ij = d_ea[iocc] + d_ea[jocc];

    for(unsigned int a = blockIdx.x*blockDim.x + threadIdx.x;
    	a < nv; a += blockDim.x * gridDim.x ){

	const double denom_ij_a = denom_ij - d_ev[a];

	//for(unsigned int c = 0; c < nv; c++ ){
	for(unsigned int c = 0; c < a; c++ ){

	    const double denom_ij_ac = denom_ij_a - d_ev[c];
	    const double denominator = 1.0/denom_ij_ac;

	    const double iajc = d_pEri[ index(a,c,nv) ];
	    const double icja = d_pEri[ index(c,a,nv) ];
	    
	    const double E1 = (2*iajc - icja) * denominator;
	    const double E2 = (2*icja - iajc) * denominator;

	    // (2*(ia|jc)-(ic|ja))/(e[iocc]+e[jocc]-e[a]-e[c])
	    //	    d_pAC[ index(a,c,nv) ] = (2*iajc - icja) * denominator;
	    // d_pAC[ index(a,c,nv) ] = E1;
	    // d_pAC[ index(c,a,nv) ] = E2;
		
	    // (2*(ic|ja)-(ia|jc))/(e[iocc]+e[jocc]-e[a]-e[c])
	    //	    d_pAC[ index(c,a,nv) ] = (2*icja - iajc) * denominator;
		
	    // (ib|jc)/(e[iocc]+e[jocc]-e[b]-e[c])  (nb. b = a)
	    //d_pBC [ index(a,c,nv) ] = iajc * denominator;
	    //d_pBC [ index(c,a,nv) ] = icja * denominator;
		
	    // (ic|jb)/(e[iocc]+e[jocc]-e[b]-e[c])  (nb. b = a)
	    //	    d_pBC [ index(c,a,nv) ] = icja * denominator;

	    //compute energy here
	    d_pEnergy[blockIdx.x*blockDim.x + threadIdx.x] += E1*iajc;
	    d_pEnergy[blockIdx.x*blockDim.x + threadIdx.x] += E2*icja;
	    
	    //build TijLong here
	    d_pTijLong[ index(a,c,nv) ] = E1;
	    d_pTijLong[ index(c,a,nv) ] = E2;

	    //build EriLong here
	    d_pEriLong[ index(a,c,nv) ] = iajc*denominator;
	    d_pEriLong[ index(c,a,nv) ] = icja*denominator;

	}//b
	
	const double denom_ij_aa = denom_ij_a - d_ev[a];
	const double denominator = 1.0/denom_ij_aa;
	const double iaja = d_pEri[ index(a,a,nv) ];
	const double E1 = iaja * denominator;
	//d_pAC[ index(a,a,nv) ] = E1;
	//d_pBC[ index(a,a,nv) ] = iaja * denominator;
	d_pEnergy[blockIdx.x*blockDim.x + threadIdx.x] += E1*iaja;
	d_pTijLong[ index(a,a,nv) ] = E1;
	d_pEriLong[ index(a,a,nv) ] = E1;

    }//a
    
};//invokeGradientPvv


extern "C" void deviceGradientPvv(const int nv,
				  const int iocc,
				  const int jocc,
				  const double *ea,
				  const double *ev,
				  //double *pAC,
				  //double *pBC,
				  double *pEri,
				  double *pEnergy,
				  double *pTijLong,
				  double *pEriLong,
				  const int cudaBlocks, 
				  const int cudaThreadsPB, 
				  const cudaStream_t stream)

{

    //    invokeGradientPvv<<<cudaBlocks,cudaThreadsPB,0,stream>>>
    invokeGradientPvv<<<cudaBlocks,cudaThreadsPB,0,stream>>>
	//(nv,iocc,jocc,ea,ev,pAC,pBC,pEri,pEnergy,pTijLong,pEriLong);
	//(nv,iocc,jocc,ea,ev,pAC,pEri,pEnergy,pTijLong,pEriLong);
	(nv,iocc,jocc,ea,ev,pEri,pEnergy,pTijLong,pEriLong);

	return;


} //deviceGradientPvv



























__global__ void invokeJKCopy(int mu,
			     size_t N,
			     size_t naux,
			     size_t mnSize,
			     double *pQmnSym,
			     double *pAOBatch)
{
    
cublasHandle_t cnpHandle;
cublasStatus_t status = cublasCreate(&cnpHandle);

for(int nu = blockIdx.x*blockDim.x + threadIdx.x;
    nu < N; 
    nu += blockDim.x * gridDim.x ){
    
    const unsigned int start = (mu >= nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu );
    
    cublasStatus_t error = 
	cublasDcopy(cnpHandle,
		    (int)naux,
		    &pQmnSym[start], (int)mnSize,
		    &pAOBatch[nu], (int)N);
    
 }//nu

cublasDestroy(cnpHandle);

}//invokeJKCopy

extern "C" void deviceJKCopy(int mu,
			     size_t N,
			     size_t naux,
			     size_t mnSize,
			     double *pQmnSym,
			     double *pAOBatch,
			     int cudaBlocks, 
			     int cudaThreadsPB, 
			     cudaStream_t stream)


{
    //    invokeJKCopy<<< cudaBlocks, cudaThreadsPB, 0, stream >>>
    //    invokeJKCopy<<< 1, cudaThreadsPB, 0, stream >>>
    //    invokeJKCopy<<< cudaBlocks, 1, 0, stream >>>
    //    invokeJKCopy<<< cudaBlocks, 32, 0, stream >>>
    //        invokeJKCopy<<< 8, 1, 0, stream >>>
    //            invokeJKCopy<<< 1, 8, 0, stream >>>
    //        invokeJKCopy<<< 1, 16, 0, stream >>>
    //    invokeJKCopy<<< 1, 2, 0, stream >>>
    //invokeJKCopy<<< 1, 1, 0, stream >>>
    invokeJKCopy<<< 1, 2, 0, stream >>>
    //invokeJKCopy<<< 1, 3, 0, stream >>>
	(mu, N, naux, mnSize, pQmnSym, pAOBatch); //, nocc, pLeft, pRight, pElp, pErp);


    return;

}//deviceJKCopy

