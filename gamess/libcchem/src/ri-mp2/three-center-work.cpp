/*
 * three-center-work.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: luke
 */

#include <three-center-work.hpp>

void ::cchem::rimp2_gradient::detail::ThreeCenterWork::
GetAOIntBatch(
	      std::vector<MapMatrixXd> &sph_c,
	      std::vector<MapMatrixXd> &threec_eri_mat_vec,
	      cchem::ri::AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri > &auxiliary_eri,
	      const Basis::Shell &L){

    BOOST_AUTO(const &basis, basis_.get());

    for(size_t s = 0; s < basis.shells().size(); s++){
	const Basis::Shell &S = basis.shells().at(s);
	
	for (size_t q = 0; q <= s; ++q) {
	    
	    const Basis::Shell &Q = basis.shells().at(q);
	    
	    if( S.size() >= Q.size()  ){
		double *aux_ptr = auxiliary_eri(Q,S,L);
		
		if(spherical_){
		    
		    MapMatrixXd threec_sph_batch(aux_ptr, Q.size()*S.size(), L.size());
		    Eigen::MatrixXd threec_sph_temp(L.sphsize(),Q.size()*S.size());
		    threec_sph_temp = threec_sph_batch*sph_c[L.Lmin()].block(0,0,L.size(),L.sphsize());
		    double *init_ptr = threec_sph_temp.data();
		    for(int l=0; l < L.sphsize(); l++){
			MapMatrixXd threec_batch(init_ptr, S.size(),Q.size());
			threec_eri_mat_vec.at(l).block(S.start(),Q.start(),S.size(),Q.size()) = threec_batch;
			threec_eri_mat_vec.at(l).block(Q.start(),S.start(),Q.size(),S.size()) = threec_batch.matrix().transpose();
			init_ptr += Q.size()*S.size();
		    }//l
		    
		}else{ //(spherical)
		    
		    for (int l=0; l<L.size(); l++){
			MapMatrixXd threec_batch(aux_ptr, S.size(), Q.size());
			threec_eri_mat_vec.at(l).block(S.start(),Q.start(),S.size(),Q.size()) = threec_batch;
			threec_eri_mat_vec.at(l).block(Q.start(),S.start(),Q.size(),S.size()) = threec_batch.matrix().transpose();
			aux_ptr += Q.size()*S.size();
		    }// l
		    
		} //(spherical)
		
	    }//( S.size() > Q.size() )
	    
	    if( Q.size() > S.size() ){
		double *aux_ptr = auxiliary_eri(S,Q,L);
		
		if(spherical_){
		    
		    MapMatrixXd threec_sph_batch(aux_ptr, Q.size()*S.size(), L.size());
		    Eigen::MatrixXd threec_sph_temp(L.sphsize(),Q.size()*S.size());
		    threec_sph_temp = threec_sph_batch*sph_c[L.Lmin()].block(0,0,L.size(),L.sphsize());
		    double *init_ptr = threec_sph_temp.data();
		    for(int l=0; l < L.sphsize(); l++){
			MapMatrixXd threec_batch(init_ptr, Q.size(),S.size());
			threec_eri_mat_vec.at(l).block(S.start(),Q.start(),S.size(),Q.size()) = threec_batch.matrix().transpose();
			threec_eri_mat_vec.at(l).block(Q.start(),S.start(),Q.size(),S.size()) = threec_batch;
			init_ptr += Q.size()*S.size();
		    }//l
		    
		}else{ //(spherical)
		    
		    for (int l=0; l<L.size(); l++){
			MapMatrixXd threec_batch(aux_ptr, Q.size(), S.size());
			threec_eri_mat_vec.at(l).block(S.start(),Q.start(),S.size(),Q.size()) = threec_batch.matrix().transpose();
			threec_eri_mat_vec.at(l).block(Q.start(),S.start(),Q.size(),S.size()) = threec_batch;
			aux_ptr += Q.size()*S.size();
		    }// l
		    
		}//(spherical)
		
	    }//( Q.size() > S.size() )
	    
	}// q
	
    }// s
    
    return;
    
}//GetAOIntBatch


void ::cchem::rimp2_gradient::detail::ThreeCenterWork::
BuildMOInts(std::vector<MapMatrixXd> &sph_c){

#pragma omp parallel
    if (pe_.node().rank() == 0) {
	
	utility::timer timerTest;
	
	detail::Thread::Task<Parallel::Task&> task(pe_.task());
#pragma omp barrier
	
	cchem::ri::AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri > 
	    auxiliary_eri(
			  boost::cref(auxbasis_),
			  boost::cref(basis_),
			  boost::cref(Ca_),
			  boost::cref(Cv_) );
	
	BOOST_AUTO(const &basis, basis_.get());
	BOOST_AUTO(const &auxbasis, auxbasis_.get());
	
	BOOST_AUTO(const &Ca, Ca_.get());
	BOOST_AUTO(const &Cv, Cv_.get());
	
	double *ptr_threec = auxiliary_eri.get_pointer(2); //N*N*max_aux_shell : for (mu nu|L) integrals
	std::vector<MapMatrixXd> threec_eri_mat_vec;
	for(int l=0;l<max_auxbasis_shell_size_;l++){
	    MapMatrixXd temp(ptr_threec +l*basis.size()*basis.size(), basis.size(), basis.size());
	    threec_eri_mat_vec.push_back(temp);
	}
	
	//ATTENTION, the get_pointer(1) does not reflect the fact that we are
	//     using spherical vs cartesian gaussians for auxiliary basis
	//     (assumes cartesian basis - larger space needed)
	//     LOW priority
	
	double *pSym = new double[auxbasis.max().size()*N_*(N_+1)/2];

	double *ptr_ab = auxiliary_eri.get_pointer(3); //nv*nv*max_aux_shell;
	MapMatrixXd t1_ab(  auxiliary_eri.get_pointer(1), nv_, basis.size() ); //N*nv : for 1/3 trans
	//reuse auxiliary_eri.get_pointer(1) for t1 (ia and ij) (its sure to be large enough)

	//pointer 1 is bigger than needed ( N*nv )
	MapMatrixXd t1(  auxiliary_eri.get_pointer(1), no_, basis.size() ); //N*nv;  : 1/3 trans

	       
	while (++task < auxbasis.shells().size() ) { //this is the L shell
	    const Basis::Shell &L = auxbasis.shells().at(task);
	    
	    //get a batch of integrals (over L shell)
	    this->GetAOIntBatch( sph_c, threec_eri_mat_vec, auxiliary_eri, L);

	    
	    //set up locations/size based on auxiliary basis set gaussians
	    int lstart = L.start();
	    int lsize = L.size();
	    int lstop = lstart+lsize;
	    if(spherical_){
		lstart = L.sphstart();
		lsize = L.sphsize();
		lstop = lstart+lsize;
	    }//(spherical_)
	    
			//pack up unique integerals
	    for(int l=0; l < lsize; l++){
		double *batch = threec_eri_mat_vec.at(l).data();
		double *pLSym = &pSym[ l * N_*(N_+1)/2 ];
		
		for(int mu = 0; mu < N_; mu++){
		    
		    size_t symStart = mu*(mu+1)/2 ;
		    cblas_dcopy( mu+1 , &batch[mu], N_, &pLSym[ symStart ], 1 );
		}//mu
		
	    }//l

	    size_t startSym[] = { 0, lstart };
	    size_t finishSym[] = { N_*(N_+1)/2, lstart+lsize };
	    BARE_MN_SYM_->put( pSym, startSym, finishSym );
	    
	    //VO vectors B_i(virt)^Q

	    double *ptr = ptr_ab;
	    for(int l=0; l<lsize; l++){
		
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			    no_, N_, N_,
			    1.0, Ca.data().begin(), no_,
			    threec_eri_mat_vec.at(l).data(), N_,
			    0.0, t1.data(), no_ );

		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
			    nv_, no_, N_,
			    1.0, Cv.data().begin(), nv_,
			    t1.data(), no_,
			    0.0, &ptr_ab[l*no_*(nv_+ns_)], nv_ );
		
	    } //l
	    size_t start[] = { 0, lstart };
	    size_t finish[] = { no_*(nv_+ns_), lstop };
	    BARE_IA_->put(ptr_ab,start,finish); // (ia|Q) - for gradient
	    
	}//++task : auxiliary shells

	delete [] pSym;

    }//(pe_.node().rank() == 0)

    return;
}//BuildMOInts






void ::cchem::rimp2_gradient::detail::ThreeCenterWork::
BuildMixedLagrangian(double *ptr_gamma_ia_P, double *ptr_gamma_inu_P,
		     std::vector<MapMatrixXd> &sph_c,
		     MapMatrixXd &lag_mu_i,
		     MapMatrixXd &lag_a_nu){

    double *pLagMuI = lag_mu_i.data();
    double *pLagANu= lag_a_nu.data();

#pragma omp parallel
    if (pe_.node().rank() == 0) { //this is the node process

	double *pThreadLagMuI = new double[N_*no_];
	memset(pThreadLagMuI,0,N_*no_*sizeof(double) );
	
	double *pThreadGammaINuP = new double[N_*no_];
	
	double *pThreadLagANu = new double[nv_*N_];
	memset(pThreadLagANu,0,nv_*N_*sizeof(double) );
	
	double *pThreadGammaIAP = new double[no_*nv_];
	
	detail::Thread::Task<Parallel::Task&> task(pe_.task());
#pragma omp barrier
	
	cchem::ri::AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri > 
	    auxiliary_eri(
			  boost::cref(auxbasis_),
			  boost::cref(basis_),
			  boost::cref(Ca_),
			  boost::cref(Cv_) );
	
	BOOST_AUTO(const &basis, basis_.get());
	BOOST_AUTO(const &auxbasis, auxbasis_.get());
	
	BOOST_AUTO(const &Ca, Ca_.get());
	BOOST_AUTO(const &Cv, Cv_.get());
	
	
	double *ptr_threec = auxiliary_eri.get_pointer(2); //N*N*max_aux_shell : for (mu nu|L) integrals
	std::vector<MapMatrixXd> threec_eri_mat_vec;
	for(int l=0;l<auxbasis.max().size();l++){
	    MapMatrixXd temp(ptr_threec +l*basis.size()*basis.size(), basis.size(), basis.size());
	    threec_eri_mat_vec.push_back(temp);
	}//l

	//pointer 1 is bigger than needed ( N*nv ) --- resused from t_ab transformation
	double *pT1 = auxiliary_eri.get_pointer(1);
	
	//get a batch of integrals
	while (++task < auxbasis.shells().size() ) { //this is the L shell
	    const Basis::Shell &L = auxbasis.shells().at(task);
	    
	    //get a batch of integrals (over L shell)
	    this->GetAOIntBatch( sph_c, threec_eri_mat_vec, auxiliary_eri, L);
	    
	    //set up locations/size based on auxiliary basis set gaussian type
	    int lstart = L.start();
	    int lsize = L.size();
	    if(spherical_){
		lstart = L.sphstart();
		lsize = L.sphsize();
	    }//(spherical_)

	    //the lag_mu_i and lag_a_nu are critical 
	    // (fix this for better parallel performance)
	    for(int l=0; l<lsize; l++){
		
		//accumlate L_mu_i(1) - no transformation needed
		size_t start2[] = {0, lstart+l};
		size_t finish2[] = {no_*N_, lstart +l +1 };
		GAMMA_INU_P_->get(pThreadGammaINuP, start2, finish2);
		
		// (N_,no_) = (N_,N_) * (N_, no_)
		//lag_mu_i += threec_eri_mat_vec.at(l)*gamma_inu_P;  
		double *pThreeCERI = threec_eri_mat_vec.at(l).data();
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			    N_,no_,N_,
			    1.0, pThreeCERI, N_,
			    pThreadGammaINuP, N_,
			    1.0, pThreadLagMuI, N_);
		
		//accumulate L_nu_a(2) - need 1/3 transformation
		//create 1/3 transform ints (mu nu|Q) -> (i nu|Q)
		
		//(no_,N_)=(no_,N_)*(N_,N_)
		//t1 = occ_coeff_mat * threec_eri_mat_vec.at(l); 
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			    no_,N_,N_,
			    1.0, Ca.data().begin(), no_,
			    pThreeCERI,N_,
			    0.0, pT1, no_);
		
		size_t start[] = {0, lstart+l};
		size_t finish[] = {no_*nv_, lstart +l +1 };
		GAMMA_IA_Q_->get(pThreadGammaIAP, start, finish);
		
		//lag_a_nu += gamma_ia_P * t1;
		// (nv_,N_) = (nv_,no_) * (no_,N_)
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
			    nv_,N_,no_,
			    -1.0, pThreadGammaIAP, nv_,
			    pT1,no_,
			    1.0, pThreadLagANu, nv_);
		
	    } //l
	    
	}//++task : auxiliary shells
	
#pragma omp critical(one)
	cblas_daxpy( N_*no_, 1.0, pThreadLagMuI, 1, pLagMuI, 1);
	
#pragma omp critical(two)
	cblas_daxpy( nv_*N_, 1.0, pThreadLagANu, 1, pLagANu, 1);
	
	delete [] pThreadLagMuI;
	delete [] pThreadLagANu;
	delete [] pThreadGammaIAP;
	delete [] pThreadGammaINuP;
	
	
    }//(pe_.node().rank() == 0)
    
    pe_.reduce("+",pLagMuI, (size_t)N_*no_); 
    pe_.reduce("+",pLagANu, (size_t)nv_*N_); 
    
    return;

}//BuildMixedLagrangian








void ::cchem::threeCenterInt::threeCenterEriWork::
getThreeCenterEriBatch(const size_t &domain,
		       const std::vector< size_t > &blockRanges,
		       cchem::ri::AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri > &auxiliary_eri,
		       std::vector<MapMatrixXd> &sph_c,
		       double *ptr_buff1,const size_t &offsetAO){
    
    //	size_t shellIndex = blockRanges[domain];
    //	const Basis::Shell &blockL = auxbasis_.shells().at( shellIndex );
    //	size_t offsetAO = (blockL.start())*N_*N_;
    //	if(spherical_) offsetAO = (blockL.sphstart())*N_*N_;
    
    
    //compute three center ERIs. Perform cartesian to spherical transformation
#pragma omp for schedule(dynamic)
    for(size_t task = blockRanges[domain]; task < blockRanges[domain+1]; task++ ){
	const Basis::Shell &L = auxbasis_.shells().at(task);

	for(size_t s = 0; s < basis_.shells().size(); s++){
	    const Basis::Shell &S = basis_.shells().at(s);
	    
	    for (size_t q = 0; q <= s; ++q) {
		const Basis::Shell &Q = basis_.shells().at(q);

		int braSize =  Q.size()*S.size();
		
		if(S.size() < Q.size()) {std::cout << "not good, S.size() is less than Q.size() "
						   << std::endl 
						   << "  did you turn off basis shell sorting?"
						   << " If not reactivate the code for the case when S.size() < Q.size()" << std::endl;}
		
		//get a triplet of ERIs ( AO AO | AUX )

		// if(S.L() >= Q.L()){
		double *aux_ptr = auxiliary_eri(Q,S,L);


		int lSize = L.size();
		if(spherical_)lSize = L.sphsize();
		
		if(spherical_){
		    double *ptr_sph = (sph_c[L.Lmin()].data());
		    //perform SALC transformation on auxiliary shell
		    //( AO AO | AUX(cartesian) ) -> ( AO AO | AUX(spherical) )
		    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				(int)braSize, (int)L.sphsize(), (int)L.size(),
				1.0, aux_ptr, (int)braSize,
				ptr_sph, (int)L.size(),
				0.0, threeCartToSph_, (int)braSize );
		    aux_ptr = threeCartToSph_;
		}//if(spherical)


		int lStart = L.start();
		if(spherical_) lStart = L.sphstart();
		//for now we do memcpy to incrementally fill ( MO MO | AUX ) matrix
		for(int l=0; l < lSize; l++){
		    
		    double *ptr_temp = ptr_buff1   + l*N_*N_ + lStart*N_*N_ - offsetAO;

		    for(int qput = 0; qput < Q.size(); qput++){
			
			std::memcpy( &ptr_temp[ (Q.start()+qput)*N_ + S.start() ],
				     &aux_ptr[l*braSize +qput*S.size()],
				     S.size()*sizeof(double) );
			
		    }//qput
		}//l

		// }else{ //(S.L() >= Q.L())

		//     double *aux_ptr = auxiliary_eri(S,Q,L);


		// int lSize = L.size();
		// if(spherical_)lSize = L.sphsize();
		
		// if(spherical_){
		//     double *ptr_sph = (sph_c[L.Lmin()].data());
		//     //perform SALC transformation on auxiliary shell
		//     //( AO AO | AUX(cartesian) ) -> ( AO AO | AUX(spherical) )
		//     cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
		// 		(int)braSize, (int)L.sphsize(), (int)L.size(),
		// 		1.0, aux_ptr, (int)braSize,
		// 		ptr_sph, (int)L.size(),
		// 		0.0, threeCartToSph_, (int)braSize );
		//     aux_ptr = threeCartToSph_;
		// }//if(spherical)


		// int lStart = L.start();
		// if(spherical_) lStart = L.sphstart();
		// //for now we do memcpy to incrementally fill ( MO MO | AUX ) matrix
		// for(int l=0; l < lSize; l++){
		//     double *ptr_temp = ptr_buff1   + l*N_*N_ + lStart*N_*N_ - offsetAO;
		//     for(int qput = 0; qput < Q.size(); qput++){
		// 	// std::memcpy( &ptr_temp[ (Q.start()+qput)*N_ + S.start() ],
		// 	// 	     &aux_ptr[l*braSize +qput*S.size()],
		// 	// 	     S.size()*sizeof(double) );
		// 	cblas_dcopy( Q.size() , &aux_ptr[l*braSize +qput], 1,
		// 		     &ptr_temp[ (Q.start()+qput)*N_ + S.start() ], N_ );
		// 	//&ptr_temp[ S.start()*N_ + Q.start()+qput ], N_ );
		//     }//qput
		// }//l

		// } //(S.L() >= Q.L())

	    }// q
	    
	}// s

    }//task : auxiliary shells
    
    
}//::cchem::threeCenterInt::threeCenterEriWork::getThreeCenterEriBatch






#if HAVE_CUBLAS

void ::cchem::threeCenterInt::threeCenterTransform::
GPUAllocate(){
    
    
    for(int idevice = 0; idevice < ndevice_; idevice++){
	
	gpu_counter.push_back(0);
	
	GPUerrchk(cudaSetDevice( idevice ));
	
	//device storage for ( mu nu |Q) integrals
	double *ptrMuNuQ = 0;
	int matrixSizeMuNuQ = N_*N_*nstreams_;
	GPUerrchk( cudaMalloc((void **)&ptrMuNuQ,matrixSizeMuNuQ*sizeof(double)) );
	dvMuNuQ.push_back(ptrMuNuQ);
	
	//device storage for (i mu|Q) 1/3 transformed integrals
	double *ptrINuQ = 0;
	int matrixSizeINuQ = no_*N_*nstreams_;
	GPUerrchk(cudaMalloc((void **)&ptrINuQ,matrixSizeINuQ*sizeof(double)));
	dvINuQ.push_back(ptrINuQ);
	
	//device storage for ( i a |Q)  2/3 transformed integrals
	double *ptrIAQ = 0;
	int matrixSizeIAQ = no_*nvs_*nstreams_;
	GPUerrchk( cudaMalloc((void **)&ptrIAQ,matrixSizeIAQ*sizeof(double)) );
	dvIAQ.push_back(ptrIAQ);
	
	//device storage for occupied LCAO coefficients
	double *ptrCoeffOcc = 0;
	int matrixSizeCoeffOcc = no_*N_;
	GPUerrchk( cudaMalloc((void **)&ptrCoeffOcc,matrixSizeCoeffOcc*sizeof(double)) );
	dvCoeffOcc.push_back(ptrCoeffOcc);
	GPUerrchk( cudaMemcpy(dvCoeffOcc[idevice],&(Ca_.data()[0]), no_*N_*sizeof(double),
			      cudaMemcpyHostToDevice) );    
	
	//device storage for singley occupied and virtual LCAO coefficients
	double *ptrCoeffSV = 0;
	int matrixSizeCoeffSV = nvs_*N_;
	GPUerrchk( cudaMalloc((void **)&ptrCoeffSV,matrixSizeCoeffSV*sizeof(double)) );
	dvCoeffSV.push_back(ptrCoeffSV);
	GPUerrchk( cudaMemcpy(dvCoeffSV[idevice],&(Cv_.data()[0]), nvs_*N_*sizeof(double),
			      cudaMemcpyHostToDevice) );    
	
    }//idevice
    
    
    
}//::cchem::threeCenterInt::threeCenterTransform::GPUAllocate





void ::cchem::threeCenterInt::threeCenterTransform::
GPUDeallocate(){
    
    for(int idevice = 0; idevice < ndevice_; idevice++){
	
	GPUerrchk( cudaSetDevice( idevice ) )
	    
	    GPUerrchk( cudaFree(dvMuNuQ[idevice]) );
	GPUerrchk( cudaFree(dvINuQ[idevice]) );
	GPUerrchk( cudaFree(dvIAQ[idevice]) );
	GPUerrchk( cudaFree(dvCoeffOcc[idevice]) );
	GPUerrchk( cudaFree(dvCoeffSV[idevice]) );
	
    }//idevice
    
}//::cchem::threeCenterInt::threeCenterTransform::GPUDeallocate


#endif //HAVE_CUBLAS


