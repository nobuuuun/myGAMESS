/*
 * JK.cpp
 *
 *  Created on: Feb 4, 2016
 *      Author: luke
 */

#include<JK.hpp>



void cchem::rimp2_gradient::detail::buildJKMatrixFunctor::
operator()(double *pQmnSym, 
	   size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
	   size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols)
{
    const size_t naux = A1_BCStop - A1_BCStart;
    
    utility::timer timerTotal;
#pragma omp master
    timerTotal.reset();
    ////////////////
    //build J matrix
    ////////////////
    
#pragma omp single nowait
    for(int iF = 0; iF < pJ_.size(); iF++){
	
	//get factored density matrix (if applicable)
	double *D = pD_[iF];
	
	//create lower triangular density matrix
	for(size_t mu = 0, munu = 0; mu < N_; mu++){
	    for(size_t nu = 0; nu <= mu ; nu++, munu++){
		pDmn_[iF][munu] = ( mu==nu ? D[mu+nu*N_] : D[mu+nu*N_] + D[nu+mu*N_] );
	    }//nu
	}//mu
	
	//contract ERIs with pseudo-density matrix
	//     DAux(naux) = (L | mu nu ) * Dmn(mu nu)
	cblas_dgemv(CblasColMajor,CblasTrans,
		    mnSize_,naux,
		    1.0, pQmnSym, mnSize_,
		    pDmn_[iF], 1,
		    0.0, pDAux_[iF], 1);
	
	//contract auxiliary functions to get back to AO basis
	//     JTemp = ( mu nu | L) * DAux(naux)
	cblas_dgemv(CblasColMajor,CblasNoTrans,
		    mnSize_,naux,
		    1.0, pQmnSym, mnSize_,
		    pDAux_[iF], 1,
		    0.0, pJTemp_, 1);
	
	//triangular J to square J
	for(size_t mu = 0, munu = 0; mu < N_; mu++){
	    for(size_t nu = 0; nu <= mu ; nu++, munu++){
		pJ_[iF][mu + nu*N_ ] += pJTemp_[munu];
		if( mu != nu ) pJ_[iF][nu + mu*N_] += pJTemp_[munu];
	    }//nu
	}//mu
	
    }//iF
    
    ////////////////
    //build K matrix
    ////////////////
    
    //Each thread allocates storage for 3-center ERI
    double *aoBatch = new double[N_*naux];

    //loop over the number of factors
    for(int iF = 0; iF < pK_.size(); iF++){
	
	const int nocc = dims_[iF];
	
#pragma omp for schedule(dynamic)
	for (int mu = 0; mu < N_; mu++){
	    
	    //unpack (mu nu| L) integrals
	    for (int nu = 0; nu < N_; nu++){
		const size_t start = (mu >= nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu );
		cblas_dcopy(naux, &pQmnSym[start], mnSize_, &aoBatch[nu], N_);
	    }
	    
	    //expand left factor of K
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			nocc, naux, N_,
			1.0, pLeft_[iF], nocc,
			aoBatch, N_,
			0.0, &pExpLeft_[mu*nocc*naux], nocc);
	    //expand right factor of K
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			nocc, naux, N_,
			1.0, pRight_[iF], nocc,
			aoBatch, N_,
			0.0, &pExpRight_[mu*nocc*naux], nocc);
	    
	}//mu
	
	//construct K matrix
#pragma omp single
	cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		    N_, N_, nocc*naux,
		    1.0, pExpLeft_, nocc*naux,
		    pExpRight_, nocc*naux,
		    1.0, pK_[iF], N_);
	
    }//iF
    
    delete [] aoBatch;
    
    //#pragma omp master
    //std::cout << "CPU Chunk Time " << timerTotal << std::endl;
    
};//operator()

















