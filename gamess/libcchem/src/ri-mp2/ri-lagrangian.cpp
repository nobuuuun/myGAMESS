/*
 * ri-lagrangian.cpp
 *
 *  Created on: Nov 11, 2015
 *      Author: luke
 */

#include <ri-lagrangian.hpp>
#include <math.hpp>

#include "god.hpp"


#define index(a,b,offset) (a + b*offset)

void ::cchem::rimp2_gradient::detail::CPUPmatGammaFunctor::operator()
    (const double *ptr_a1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
     size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
     const double *ptr_a2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
     size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols){

    const size_t rows1 = (A1_BRStop-A1_BRStart)*NA1BRows;
    const size_t rows2 = (A2_BRStop-A2_BRStart)*NA2BRows;

    //threads contribution to energy
    double &energy = threadEnergy_[omp_get_thread_num()];

    //nv_*nv_
    double *pEri = pVecEri_[omp_get_thread_num()];

    //no_*no_
    double *pPmatOcc = pVecPmatOcc_[omp_get_thread_num()];

    //nv_*nv_
    double *pPmatVirt = pVecPmatVirt_[omp_get_thread_num()];

    //nv_*nl_
    double *pGamma = pVecGamma_[omp_get_thread_num()];

    //nl_*nl_
    double *pGammaRS = gamma_rs_;

    //nv_*nv_
    //double *AC = pVecAC_[omp_get_thread_num()];

    //nv_*nv_
    //double *BC = pVecBC_[omp_get_thread_num()];

    for(size_t iocc = A1_BRStart; iocc < A1_BRStop; iocc++){
	const size_t a1off = iocc-A1_BRStart;

	//if(A2_BRStart == nf_) memset(&pGamma[0], 0, nv_*nl_*sizeof(double) );
	if(A2_BRStart == nf_) memset(pGamma, 0, nv_*nl_*sizeof(double) );

	//for this to work correctly (for now), you need the entire set of jocc
	//#pragma omp for schedule(dynamic) //nowait
#pragma omp for schedule(dynamic) nowait
	for(size_t jocc = A2_BRStart; jocc < A2_BRStop; jocc++){
	    const size_t a2off = jocc-A2_BRStart;

	    //build ERIs : ( iocc a | P ) * (P|Q)^-1 * (Q | jocc nv)
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
			NA1BRows,NA2BRows,NA1BCols,
			1.0, &ptr_a1[a1off*NA1BRows], rows1,
			&ptr_a2[a2off*NA2BRows], rows2,
			0.0, pEri, NA1BRows);

	    const double denom = ea_[iocc] + ea_[jocc];

	    //////////////////////////////////
	    //start pPmatVirt contribution
	    //////////////////////////////////
	    //pPmatVirt[ index[ (a,b,nv_) ] += ( 2*pEri[ index(a,c,nv_) ] - pEri[ index(c,a,nv_) ]) *
	    //                pEri[ index(b,c,nv_) ]/(denom_m_a - ev_[c]) * (denom_m_b - ev_[c]) );
	    //
	    //  since the denomenator is separable use a DGEMM, pPmatVirt +=  AC*BC		
	    //
	    // build AC and BC

	    double *pEriLong = pVecEriLong_[jocc-nf_];
	    double *pTijLong = pVecTijLong_[jocc-nf_];

	    for (int a=0; a<nv_; a++){
		const double denom_m_a = denom - ev_[a];
		
		for(int c = 0; c <a ;c++){
		    const double denominator = 1.0/(denom_m_a - ev_[c]);
		    // const double iajc = pEri[ index(a,c,nv_) ];
		    // const double icja = pEri[ index(c,a,nv_) ];
		    
		    // // (2*(ia|jc)-(ic|ja))/(e[iocc]+e[jocc]-e[a]-e[c])
		    // AC[ index(a,c,nv_) ] = (2*iajc - icja) * denominator;
		    
		    // // (2*(ic|ja)-(ia|jc))/(e[iocc]+e[jocc]-e[a]-e[c])
		    // AC[ index(c,a,nv_) ] = (2*icja - iajc) * denominator;
		    
		    // // (ib|jc)/(e[iocc]+e[jocc]-e[b]-e[c])  (nb. b = a)
		    // BC [ index(a,c,nv_) ] = iajc * denominator;
		    
		    // // (ic|jb)/(e[iocc]+e[jocc]-e[b]-e[c])  (nb. b = a)
		    // BC [ index(c,a,nv_) ] = icja * denominator;


		    const double iajc = pEri[ index(a,c,nv_) ];
		    const double icja = pEri[ index(c,a,nv_) ];

		    const double E1 = (2*iajc - icja) * denominator;
		    const double E2 = (2*icja - iajc) * denominator;

		    // (2*(ia|jc)-(ic|ja))/(e[iocc]+e[jocc]-e[a]-e[c])
		    //AC[ index(a,c,nv_) ] = 2*iajc - icja;
		    //AC[ index(a,c,nv_) ] = E1;
		    
		    // (2*(ic|ja)-(ia|jc))/(e[iocc]+e[jocc]-e[a]-e[c])
		    //AC[ index(c,a,nv_) ] = 2*icja - iajc;
		    //AC[ index(c,a,nv_) ] = E2;
		    
		    // (ib|jc)/(e[iocc]+e[jocc]-e[b]-e[c])  (nb. b = a)
		    //BC [ index(a,c,nv_) ] = iajc * denominator;
		    
		    // (ic|jb)/(e[iocc]+e[jocc]-e[b]-e[c])  (nb. b = a)
		    //BC [ index(c,a,nv_) ] = icja * denominator;
		    
		    energy += E1*iajc
			+E2*icja;

		    pTijLong[ index(a,c,nv_) ] = E1;
		    pTijLong[ index(c,a,nv_) ] = E2;

		    pEriLong[ index(a,c,nv_) ] = iajc * denominator;
		    pEriLong[ index(c,a,nv_) ] = icja * denominator;
		}//c

		const double denominator = 1.0/(denom_m_a - ev_[a]);
		const double iaja = pEri[ index(a,a,nv_) ];
		//AC[ index(a,a,nv_) ] = iaja * denominator;
		//BC[ index(a,a,nv_) ] = iaja * denominator;
		energy += iaja*iaja* denominator;
		pTijLong[ index(a,a,nv_) ] = iaja* denominator;
		pEriLong[ index(a,a,nv_) ] = iaja* denominator;

	    }//a
	    //memcpy(pEriLong, BC, nv_*nv_*sizeof(double));
	    //cblas_daxpy(nv_*nv_, 2.0, BC, 1, 1.0, AC, 1);

	    



	    //(2*(ia|jc)-(ic|ja))*(ib|jc)/( (e[iocc]+e[jocc]-e[b]-e[c])*(e[iocc]+e[jocc]-e[a]-e[c]) )
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
			nv_,nv_,nv_,
			//1.0, AC, nv_,
			1.0, pTijLong, nv_,
			//BC, nv_,
			pEriLong, nv_,
			1.0, pPmatVirt, nv_);
	    //////////////////////////////////
	    //done with pPmatVirt contribution
	    //////////////////////////////////

	    // double *pEriLong = pVecEriLong_[jocc-nf_];
	    // double *pTijLong = pVecTijLong_[jocc-nf_];

	    // for(int a=0; a<nv_; a++){
	    // 	const double denom_m_a = denom - ev_[a];
		
	    // 	for(int b=0; b<nv_; b++){
	    // 	    const double denom_m_ab = 1.0/(denom_m_a - ev_[b]);
	    // 	    const double ab = pEri[ index(a,b,nv_) ];
	    // 	    const double ba = pEri[ index(b,a,nv_) ];
	    // 	    const double element = ( 2*ab - ba )*denom_m_ab;

	    // 	    energy += ab*element;
	    // 	    pTijLong[ index(a,b,nv_) ] = element;
	    // 	    //modify ERIs to compute paa block pf P(2)
	    // 	    // (ia|jb)/(ea[iocc]+ea[jocc]-ev[a]-ev[b]) - amplitudes
	    // 	    pEriLong[ index(a,b,nv_) ] = ab*denom_m_ab;

	    // 	}//b
	    // }//a
	    
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			nv_,NA2BCols,nv_,
			1.0, pTijLong, nv_,
			&ptr_a2[a2off*NA2BRows], rows2,
			1.0, pGamma, nv_);
	    
	}//jocc

	//the next only happens when we reached the last of array 2 (CIAQ)
	if(A2_BRStop == no_){
	    
	    //accumulate pGamma to pGammaLong (one at a time)
#pragma omp critical
	    cblas_daxpy(nv_*nl_,1.0,pGamma,1,pGammaLong_[iocc-nf_],1);
	    
	    //make sure all threads are done with pVecTijLong[], pVecEriLong_[], &
	    //                                    pGammaLong_[]
#pragma omp barrier
	    
	    // contribution to occupied density (i.e. pPmatOcc)
#pragma omp single nowait
	    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
			no_-nf_, no_-nf_, nv_*nv_,
			-1.0, pVecTijLong_[0], nv_*nv_,
			pVecEriLong_[0], nv_*nv_,
			1.0, pPmatOcc, no_-nf_);
	    
#pragma omp single
	    {
		const double *pGammaLong = pGammaLong_[iocc-nf_];
		cblas_dsymm(CblasColMajor,CblasRight,CblasLower,
			    nv_,nl_,
			    1.0, data_pqm1_, nl_,
			    &ptr_a1[a1off*NA1BRows], rows1,
			    0.0, pTrans_, nv_);
		cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
			    nl_,nl_,nv_,
			    2.0, pGammaLong, nv_,
			    pTrans_, nv_,
			    1.0, pGammaRS, nl_);
	    }//omp single
	    
	}//(A2_BRStop == no_)

    }//iocc

#pragma omp barrier


}//CPUPmatGammaFunctor









void ::cchem::rimp2_gradient::detail::matrix_factor::
pseudoSquareRoot(const double &cutoff){
    
    Eigen::MatrixXd scratchDensity(no_+nv_,no_+nv_);
    
    scratchDensity = densityMatrix_.block(0,0,no_+nv_,no_+nv_);
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(scratchDensity);
    
    Eigen::VectorXd evalues(no_+nv_);
    Eigen::MatrixXd evectors(no_+nv_,no_+nv_);
    
    evalues  = eigensolver.eigenvalues();
    evectors = eigensolver.eigenvectors();
   
    Npi_ = 0;
    Nni_ = 0;
    
    for (int i = 0; i < evalues.size(); i++){
	if(fabs(evalues(i)) > cutoff){
	    if (evalues(i) >= 0.0){
		Npi_++;
	    }else{
		Nni_++;
	    }//logic
	}//logic
    }//i

    Eigen::MatrixXd Pos(no_+nv_,Npi_);
    Eigen::MatrixXd Neg(no_+nv_,Nni_);
    
    int Pcounter = 0;
    int Ncounter =0;
    
    for(int i = 0; i < no_+nv_; i++){
	if(fabs(evalues(i)) > cutoff){
	    if(evalues(i) >= 0.0){
		double val = sqrt(fabs(evalues(i)));
		Pos.block(0,Pcounter,no_+nv_,1) = val*evectors.block(0,i,no_+nv_,1);
		Pcounter++;
	    }else{
		double val = sqrt(fabs(evalues(i)));
		Neg.block(0,Ncounter,no_+nv_,1) = -val*evectors.block(0,i,no_+nv_,1);
		Ncounter++;
	    }//logic
	}//logic
    }//i
    
    //transform occupied block
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		Npi_, N_, no_,
		1.0, Pos.data(), no_+nv_,
		pCa_, no_,
		0.0, pLeft_[0], Npi_);
    //transform virtual block
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		Npi_, N_, nv_,
		1.0, &Pos.data()[no_], no_+nv_,
		pCv_, nv_,
		1.0, pLeft_[0], Npi_);
    
    //transform occupied block
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		Nni_, N_, no_,
		1.0, Neg.data(), no_+nv_,
		pCa_, no_,
		0.0, pLeft_[1], Nni_);
    //transform virtual block
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		Nni_, N_, nv_,
		1.0, &Neg.data()[no_], no_+nv_,
		pCv_, nv_,
		1.0, pLeft_[1], Nni_);
    
    memcpy(pRight_[0], pLeft_[0], Npi_*N_*sizeof(double) );
    memcpy(pRight_[1], pLeft_[1], Nni_*N_*sizeof(double) );
    
}//pseudoSquareRoot


void ::cchem::rimp2_gradient::detail::matrix_factor::buildJK()
{

    //utility::timer timerTest;

    std::vector<int> dims;
    dims.clear();
    dims.push_back(Npi_);
    dims.push_back(Nni_);
    
    //construct density from +/- factorization
    for(int iF = 0; iF < 2; iF++){
	cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		    N_, N_, dims[iF],
		    1.0, pLeft_[iF], dims[iF],
		    pRight_[iF], dims[iF],
		    0.0, pD_[iF], N_);
    }//iF

    memset( pJ_[0],0,N_*N_*sizeof(double) );
    memset( pJ_[1],0,N_*N_*sizeof(double) );
    
    memset( pK_[0],0,N_*N_*sizeof(double) );
    memset( pK_[1],0,N_*N_*sizeof(double) );
    
    const size_t mnSize = N_*(N_+1)/2;
    size_t NMB = 2000;

#if HAVE_CUBLAS
    size_t GPUMem;
    if(pGPU_ != NULL )
	GPUMem = ::cchem
	    ::rimp2_gradient
	    ::detail
	    ::GPUBuildJKMatrixFunctor
	    ::GPUMemLimit(nl_, dims, N_, mnSize, pGPU_->getNumStreams(),2 );
    
    size_t NMBDeviceTemp = GPUMem;
    for(size_t ip = 0; ip < pe_.size(); ip++){
	
	if(ip == pe_.rank()){
	    pe_.broadcast(&GPUMem, 1, ip);
	}else{
	    pe_.broadcast(&NMBDeviceTemp, 1, ip);
	}
	if(NMBDeviceTemp < GPUMem) GPUMem = NMBDeviceTemp;
	
    } //ip
    NMB = GPUMem;
#endif
    
    bool ReadRow = 0; bool debug = 1; bool printArray = 0;
    ::cchem::rimp2_gradient::detail::DFAsync 
	  JMatrix(NMB,
		  LM1_MN_SYM_, 1, 0, nl_, 0,
		  ReadRow, debug, pe_, printArray);
    
    //JK on CPU
    if(pGPU_ == NULL)
	{
	    ::cchem::rimp2_gradient::detail::buildJKMatrixFunctor 
		Functor(mnSize, N_, pJTemp_, pD_, pJ_, pDmn_, pDAux_, 
			pK_, pLeft_, pRight_, dims,JMatrix.getMaxNAB(0) );
	    JMatrix.DoOpAsync_R_NODE_PARALLEL( Functor );
	}
    
#if HAVE_CUBLAS
    //JK with GPU
    if(pGPU_ != NULL)
	{//scope
	    //timerTest.reset();
	    //int nThreads = omp_get_max_threads();
	    //omp_set_num_threads(1);
	    {
		int myRank = pe_.rank();
		::cchem::rimp2_gradient::detail::GPUBuildJKMatrixFunctor
		      Functor(mnSize, N_, no_, pJTemp_, pD_, pJ_, pDmn_, pDAux_,
			      pK_, pLeft_, pRight_, dims, JMatrix.getMaxNAB(0), pGPU_,
			      myRank);
		JMatrix.DoOpAsync_R_NODE_PARALLEL( Functor );
	    }
	    //omp_set_num_threads(nThreads);
	    //std::cout  << "  JK matrix with I/O: "<< timerTest << std::endl;
	}//scope
#endif
    
    pe_.reduce("+",pJ_[0], (size_t)N_*N_); 
    pe_.reduce("+",pJ_[1], (size_t)N_*N_); 
    pe_.reduce("+",pK_[0], (size_t)N_*N_); 
    pe_.reduce("+",pK_[1], (size_t)N_*N_); 
    
}//buildJK


void ::cchem::rimp2_gradient::detail::matrix_factor::buildWmat3(MapMatrixXd &wmat, MapMatrixXd &FAO)
{

    cblas_daxpy(N_*N_, -1.0, pJ_[1], 1, pJ_[0], 1);
    cblas_daxpy(N_*N_, 2.0, pJ_[0], 1, FAO.data(), 1);
    
    cblas_daxpy(N_*N_, -1.0, pK_[1], 1, pK_[0], 1);
    cblas_daxpy(N_*N_, -1.0, pK_[0], 1, FAO.data(), 1);
    
    //	wmat +=occ_coeff_mat*FAO*occ_coeff_mat.transpose();
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
		N_,no_,N_,
		1.0, FAO.data(), N_,
		pCa_, no_,
		0.0, pScrNno_, N_);
    
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		no_,no_,N_,
		(-0.5), pCa_, no_,
		pScrNno_, N_,
		1.0, wmat.data(), N_);
    
}//buildWmat3



void ::cchem::rimp2_gradient::detail::matrix_factor::buildLagMO(MapMatrixXd &lag_mo, MapMatrixXd &FAO)
{

    cblas_daxpy(N_*N_, -1.0, pJ_[1], 1, pJ_[0], 1);
    cblas_daxpy(N_*N_, 2.0, pJ_[0], 1, FAO.data(), 1);
    
    cblas_daxpy(N_*N_, -1.0, pK_[1], 1, pK_[0], 1);
    cblas_daxpy(N_*N_, -1.0, pK_[0], 1, FAO.data(), 1);
    
    //	lag_mo +=virt_coeff_mat*FAO*occ_coeff_mat.transpose();
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
		N_,no_,N_,
		1.0, FAO.data(), N_,
		pCa_, no_,
		0.0, pScrNno_, N_);
    
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		nv_,no_,N_,
		1.0, pCv_, nv_,
		pScrNno_, N_,
		1.0, lag_mo.data(), nv_);
    
}//buildLagMO

template<typename T>
void ::cchem::rimp2_gradient::detail::matrix_factor::print_matrix(T & mat){

	std::cout << std::right << std::setprecision(10);
	for (int i = 0; i < mat.cols() ; i+=5){

		int ij_max = std::min( (int)mat.cols(), (int)(i+5) );
		std::cout << std::endl<< std::endl<< std::endl<< std::endl;
		for (int j = i; j < ij_max; j++){
			if(j == i)std::cout << std::setw(8) << j ;
			if(j > i)std::cout << std::setw(14) << j ;
		}
		std::cout<<std::endl << std::endl;
		std::cout << mat.block(0,i,mat.rows(),ij_max-i);
		std::cout<<std::endl << std::endl;
	}//i

}//print_matrix







//using ::cchem::rimp2_gradient::detail::rimp2_lagrangian;

//inline
//void ::cchem::rimp2_gradient::detail::rimp2_lagrangian::vvvo
//(const MapMatrixXd &eri_vvvo, const MapMatrixXd &pab, MapMatrixXd &lag, const int i, const int b){
//
//	for (int a = 0; a < nv_; a++){
//		for (int c = 0; c < nv_; c++){
//
//			lag(a,i) += pab(b+no_,c+no_)*2.0*eri_vvvo(a,c);
//			lag(a,i) -= pab(b+no_,c+no_)*0.5*eri_vvvo(c,a);
//			lag(b,i) -= pab(a+no_,c+no_)*0.5*eri_vvvo(a,c);
//
//		}//c
//	}//b
//}//vvvo
//
//
//inline
//void ::cchem::rimp2_gradient::detail::rimp2_lagrangian::vooo
//(const MapMatrixXd &eri_vooo, const MapMatrixXd &pij, MapMatrixXd &lag, const int i, const int j){
//
//	for (int a = 0; a < nv_; a++){
//		for(int k = 0; k < no_; k++){
//
//			lag(a,i) += pij(j,k)*2.0*eri_vooo(a,k);
//			lag(a,j) -= pij(i,k)*0.5*eri_vooo(a,k);
//			lag(a,k) -= pij(j,i)*0.5*eri_vooo(a,k);
//
//		}//k
//	}//a
//
//}//vooo







// //this is used int RI-MP2 electronic hessian (for zvector)
// void ::cchem::rimp2_gradient::detail::lag_vovo_functor::operator()
// (MapMatrixXd &eri,
// 		MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 		MapMatrixXd &bjk, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
// 		bool &restriction){


// //	void * aptr_unxt_thread = NULL;
// //	posix_memalign(&aptr_unxt_thread, 16, nv_*no_*sizeof(double) );
// //	double *ptr_unxt_thread = new(aptr_unxt_thread) double[nv_*no_];
// 	double *ptr_unxt_thread = new double[nv_*no_];
// 	MapMatrixXd unxt_thread(ptr_unxt_thread,nv_,no_);
// 	memset(ptr_unxt_thread, 0, nv_*no_*sizeof(double) );
// //	unxt_thread.setZero();

// #pragma omp for schedule(dynamic) nowait
// 	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
// 		size_t a1off = iocc-a1_start;

// 		int a2_start_prime = (restriction ? iocc : a2_start);

// 		for(size_t jocc = a2_start_prime; jocc < a2_stop; jocc++){
// 			size_t a2off = jocc-a2_start;

// 			eri = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
// 					bjk.block(a2off*NA2BRows,0,NA2BRows,NA2BCols).transpose();

// 			//			for(int a = 0; a < nv_; a++){
// 			//				for(int b = 0; b < nv_; b++){
// 			//					unxt_thread(a,iocc) += 4.0*eri(a,b)*u_(b,jocc);
// 			//					unxt_thread(a,iocc) -= eri(b,a)*u_(b,jocc);
// 			//				}//b
// 			//			}//a

// 			unxt_thread.col(iocc) += 4.0*eri*u_.col(jocc);
// 			unxt_thread.col(iocc) -= eri.transpose()*u_.col(jocc);


// 			//			if(iocc != jocc){
// 			//				for(int a = 0; a < nv_; a++){
// 			//					for(int b = 0; b < nv_; b++){
// 			//						unxt_thread(b,jocc) += 4.0*eri(a,b)*u_(a,iocc);
// 			//						unxt_thread(b,jocc) -= eri(b,a)*u_(a,iocc);
// 			//					}//b
// 			//				}//a
// 			//			}//(iocc != jocc)

// 			if(iocc != jocc){
// 				unxt_thread.col(jocc) += 4.0*eri.transpose()*u_.col(iocc);
// 				unxt_thread.col(jocc) -= eri*u_.col(iocc);
// 			}

// 		}//jocc
// 	}//iocc

// #pragma omp critical
// 	unxt_ += unxt_thread;

// 	delete [] ptr_unxt_thread;

// }//lag_vovo_functor








// //	inline --this is bad
// //this is used int RI-MP2 electronic hessian (for zvector)
// void ::cchem::rimp2_gradient::detail::lag_vvoo_functor::operator()
// 													(MapMatrixXd &eri,
// 															MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 															MapMatrixXd &bjk, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
// 															bool &loop_restriction){

// 	Eigen::MatrixXd matter(NA1BRows,NA1BCols);
// 	Eigen::MatrixXd matter2(NA2BRows,NA2BCols);

// #pragma omp for collapse(2) schedule(dynamic) nowait
// 	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
// 		for(size_t a = a2_start; a < a2_stop; a++){

// 			size_t a1off = iocc-a1_start;
// 			size_t a2off = a - a2_start;

// 			eri = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
// 					bjk.block(a2off*NA2BRows,0,NA2BRows,NA2BCols).transpose();


// 			//the following is temporary

// 			//			matter = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols);


// 			//			matter2 = bjk.block(a2off*NA2BRows,0,NA2BRows,NA2BCols);


// 			//			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
// 			//					(int)(NA1BRows),(int)(NA2BRows),(int)NA1BCols,
// 			//					1.0,matter.data(),(int)(NA1BRows),
// 			//					matter2.data(),(int)(NA2BRows),
// 			//					0.0,eri.data(),(int)NA1BRows);


// 			//			for(int jocc =0; jocc < no_; jocc++){
// 			//				for(int b = 0; b < nv_; b++){
// 			//					unxt_(a,iocc) -= eri(jocc,b)*u_(b,jocc);
// 			////					unxt_(a,iocc) -= eri(jocc,b)*u_(b,jocc);
// 			////					tij.cwiseProduct(eri).sum();
// 			//
// 			//				}//b
// 			//			}//jocc

// 			//			unxt_(a,iocc) -= eri.cwiseProduct(u_.transpose()).sum();
// 			unxt_(iocc,a) -= eri.cwiseProduct(u_).sum();

// 		}//a
// 	}//iocc

// }//lag_vvoo_functor









// //	inline --this is bad
// //used to compute vooo contributions to Lagrangian (not for zvector)
// void ::cchem::rimp2_gradient::detail::lag_vooo_functor::operator()
// 													(MapMatrixXd &eri_vooo,
// 															MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 															MapMatrixXd &bjk, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
// 															bool &loop_restriction){

// 	double *ptr_lag_thread = new double[nv_*no_];
// 	MapMatrixXd lag_thread(ptr_lag_thread, nv_ ,no_);
// 	memset(ptr_lag_thread, 0, nv_*no_*sizeof(double) );

// #pragma omp for collapse(2) schedule(dynamic)
// 	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
// 		for(size_t jocc = a2_start; jocc < a2_stop; jocc++){

// 			size_t a2off = jocc-a2_start;
// 			size_t a1off = iocc-a1_start;

// 			eri_vooo = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
// 					bjk.block(a2off*NA2BRows,0,NA2BRows,NA2BCols).transpose();

// 			//			for (int a = 0; a < nv_; a++){
// 			//				for(int k = 0; k < no_; k++){
// 			//
// 			//					lag_thread(a,iocc) += pmat_(jocc,k)*2.0*eri_vooo(a,k);
// 			//					lag_thread(a,jocc) -= pmat_(iocc,k)*0.5*eri_vooo(a,k);
// 			//					lag_thread(a,k) -= pmat_(jocc,iocc)*0.5*eri_vooo(a,k);
// 			//
// 			//				}//k
// 			//			}//a

// 			lag_thread.col(iocc) += 2.0*pmat_.row(jocc)*eri_vooo.transpose();
// 			lag_thread.col(jocc) -= 0.5*pmat_.row(iocc)*eri_vooo.transpose();
// 			lag_thread -= 0.5*pmat_(jocc,iocc)*eri_vooo;

// 		}//jocc
// 	}//iocc

// #pragma omp critical
// 	lag_mo_ += lag_thread;

// 	delete [] ptr_lag_thread;

// }//lag_vooo_functor




// //	inline --this is bad
// //used to compute vvvo contributions to Lagrangian (not for zvector)
// void ::cchem::rimp2_gradient::detail::lag_vvvo_functor::operator()
// 													(MapMatrixXd &eri_vvvo,
// 															MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 															MapMatrixXd &bba, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
// 															bool &loop_restriction){

// #pragma omp for  schedule(dynamic) //do not collapse(2)?
// 	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
// 		size_t a1off = iocc-a1_start;

// 		for(size_t bvir = a3_start; bvir < a3_stop; bvir++){
// 			size_t a3off = bvir-a3_start;

// 			eri_vvvo = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
// 					bba.block(a3off*NA3BRows,0,NA3BRows,NA3BCols).transpose();

// 			//			for (int a = 0; a < nv_; a++){
// 			//				for (int c = 0; c < nv_; c++){
// 			//
// 			//					lag_mo_(a,iocc) += pmat_(bvir+no_,c+no_)*2.0*eri_vvvo(a,c);
// 			//					lag_mo_(a,iocc) -= pmat_(bvir+no_,c+no_)*0.5*eri_vvvo(c,a);
// 			//					lag_mo_(bvir,iocc) -= pmat_(a+no_,c+no_)*0.5*eri_vvvo(a,c);
// 			//
// 			//				}//c
// 			//			}//a

// 			lag_mo_.col(iocc) += 2.0*pmat_.block(no_+bvir,no_,1,nv_)*eri_vvvo.transpose();
// 			lag_mo_.col(iocc) -= 0.5*pmat_.block(no_+bvir,no_,1,nv_)*eri_vvvo;
// 			lag_mo_(bvir,iocc) -= 0.5*pmat_.block(no_,no_,nv_,nv_).cwiseProduct(eri_vvvo).sum();


// 		}//bvir
// 	}//iocc

// }//lag_vvvo_functor


