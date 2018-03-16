
#include <config.h>

#if HAVE_CUBLAS

#include <deviceGradientEngine.hpp>


void ::cchem::rimp2_gradient::detail::GPUPmatGammaFunctor::operator()
    (const double *ptr_a1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
     size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
     const double *ptr_a2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
     size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols){

    const double one = 1.0;
    const double two = 2.0;
    const double zero = 0.0;
    const double minusOne = -1.0;
    
    int hardDevice = 0;
    const int tid = omp_get_thread_num();
    
    for(size_t iocc = A1_BRStart; iocc < A1_BRStop; iocc++){
	
	GPUerrchk(cudaSetDevice( hardDevice ));
	
	//this is only done once per iocc
	if(A2_BRStart == nf_)
	    {

		if(tid == 0)
		    {
			this->GPUInitialize(iocc, hardDevice, 
					    NA1BRows, ptr_a1,
					    one, zero,
					    A1_BRStop, A1_BRStart);
		    }else{

		}
		
		// const size_t hostChunkStart = a1off * NA1BRows;		
		// const size_t nBlocksA1 = A1_BRStop-A1_BRStart; 
	
		// //HOST --> GPU
		// GPUerrchk( cudaMemcpy2DAsync(dVecPiaQ_[hardDevice], 
		// 			     NA1BRows*sizeof(double),
		// 			     &ptr_a1[ hostChunkStart ],
		// 			     nBlocksA1*NA1BRows*sizeof(double),
		// 			     NA1BRows*sizeof(double), //width
		// 			     nl_,      //height
		// 			     cudaMemcpyHostToDevice) );
		
		// //GPU: transform ptr_a1 --> CiaQ
		// cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
		// 	    nv_, nl_, nl_, 
		// 	    &one, dVecPiaQ_[hardDevice], nv_, 
		// 	    dVecPPqm1_[hardDevice],nl_,
		// 	    &zero, dVecTrans_[hardDevice], nv_);
		
		// //GPU: initialize gamma
		// for(int istream = 0; istream < nstreams_; istream++){
		//     GPUerrchk( cudaMemsetAsync(dVecPGamma_[hardDevice][istream],
		// 			       0, 
		// 			       nv_*nl_*sizeof(double),
		// 			       dv1_streams_[hardDevice][istream]) );
		// }//istream

	    }//(A2_BRStart == nf_)


	//for this to work correctly (for now), you need the entire set of jocc
#pragma omp for schedule(dynamic)
	for(size_t jocc = A2_BRStart; jocc < A2_BRStop; jocc++){

	    this->GPUWork(jocc, A2_BRStart, NA2BRows, 
			  A2_BRStop, hardDevice,
			  NA1BRows, one, zero,
			  iocc, ptr_a2);

	    // const size_t a2off = jocc-A2_BRStart;
	    // const size_t hostCoeffStart = a2off * NA2BRows;
	    // const size_t nBlocksA2 = A2_BRStop-A2_BRStart;
	    // int istream = jocc%nstreams_;
	    
	    // //GPU:
	    // cublasSetStream( dv1_handle_[hardDevice],
	    // 		     dv1_streams_[hardDevice][istream] );
	    
	    // //HOST --> GPU
	    // GPUerrchk( cudaMemcpy2DAsync(dVecPCiaQ_[hardDevice][istream], 
	    // 				 NA1BRows*sizeof(double),
	    // 				 &ptr_a2[ hostCoeffStart ],
	    // 				 nBlocksA2*NA1BRows*sizeof(double),
	    // 				 NA1BRows*sizeof(double), //width
	    // 				 nl_,      //height
	    // 				 cudaMemcpyHostToDevice) );
	    
	    // //GPU: build ERIs
	    // cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_T, 
	    // 		nv_, nv_, nl_, 
	    // 		&one, dVecPiaQ_[hardDevice], nv_, 
	    // 		dVecPCiaQ_[hardDevice][istream],nv_,
	    // 		&zero, dVecPEri_[hardDevice][istream], nv_);
	    
	    
	    // //GPU: build intermediates
	    // deviceGradientPvv(nv_, 
	    // 		      iocc, 
	    // 		      jocc, 
	    // 		      dVecEA_[hardDevice],
	    // 		      dVecEV_[hardDevice],
	    // 		      dVecPAC_[hardDevice][istream],
	    // 		      dVecPBC_[hardDevice][istream],
	    // 		      dVecPEri_[hardDevice][istream],
			      
	    // 		      dVecPEnergy_[hardDevice][istream],
	    // 		      dVecPTijLong_[hardDevice][istream],
	    // 		      dVecPEriLong_[hardDevice][istream],
			      
	    // 		      cudaBlocks_,
	    // 		      cudaThreadsPB_,
	    // 		      dv1_streams_[hardDevice][istream] );
	    
	    // // //unified memory implementation
	    // // //build intermediates
	    // // deviceGradientPvv(nv_, 
	    // // 		      iocc, 
	    // // 		      jocc, 
	    // // 		      dVecEA_[hardDevice],
	    // // 		      dVecEV_[hardDevice],
	    // // 		      dVecPAC_[hardDevice][istream],
	    // // 		      dVecPBC_[hardDevice][istream],
	    // // 		      dVecPEri_[hardDevice][istream],
	    
	    // // 		      dVecPEnergy_[hardDevice][istream],
	    // // 		      &pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ + (no_-nf_)*nv_*nv_ ],
	    // // 		      &pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ ],
	    
	    // // 		      cudaBlocks_,
	    // // 		      cudaThreadsPB_,
	    // // 		      dv1_streams_[hardDevice][istream] );

	    // //GPU: accumlate P(vv) here
	    // cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_T, 
	    // 		nv_, nv_, nv_, 
	    // 		&one, dVecPAC_[hardDevice][istream], nv_, 
	    // 		dVecPBC_[hardDevice][istream], nv_,
	    // 		&one, dVecPPmatVirt_[hardDevice][istream], nv_);

	    // //GPU: build part of gamma
	    // cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
	    // 		nv_, nl_, nv_, 
	    // 		&one, dVecPTijLong_[hardDevice][istream], nv_, 
	    // 		dVecPCiaQ_[hardDevice][istream], nv_,
	    // 		&one, dVecPGamma_[hardDevice][istream], nv_);

	    // // //unified memory implementation
	    // // //build part of gamma
	    // // cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
	    // // 		nv_, nl_, nv_, 
	    // // 		&one, &pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ + (no_-nf_)*nv_*nv_ ], nv_, 
	    // // 		dVecPCiaQ_[hardDevice][istream], nv_,
	    // // 		&one, dVecPGamma_[hardDevice][istream], nv_);


	    // //GPU --> HOST: copy EriLong
	    // GPUerrchk( cudaMemcpyAsync(hVecPEriLong_[jocc-nf_], 
	    // 			       dVecPEriLong_[hardDevice][istream], 
	    // 			       nv_*nv_*sizeof(double),
	    // 			       cudaMemcpyDeviceToHost,
	    // 			       dv1_streams_[hardDevice][istream]) );

	    // //GPU --> HOST: copy TijLong
	    // GPUerrchk( cudaMemcpyAsync(hVecPTijLong_[jocc-nf_], 
	    // 			       dVecPTijLong_[hardDevice][istream], 
	    // 			       nv_*nv_*sizeof(double),
	    // 			       cudaMemcpyDeviceToHost,
	    // 			       dv1_streams_[hardDevice][istream]) );




	    // // //experiment
	    // // //copy EriLong form device to unified
 	    // // GPUerrchk( cudaMemcpyAsync(&pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ ], //pUnifiedEriLong_[jocc-nf_], 
	    // // 			  dVecPEriLong_[hardDevice][istream], 
	    // // 			  nv_*nv_*sizeof(double),
	    // // 			  cudaMemcpyDefault,
	    // // 			  dv1_streams_[hardDevice][istream]) );

	    // // //experiment
	    // // //copy TijLong form device to unified
	    // // GPUerrchk( cudaMemcpyAsync(&pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ + (no_-nf_)*nv_*nv_ ], 
	    // // 			       dVecPTijLong_[hardDevice][istream], 
	    // // 			       nv_*nv_*sizeof(double),
	    // // 			       cudaMemcpyDefault,
	    // // 			       dv1_streams_[hardDevice][istream]) );
	    
	}//jocc



	if(A2_BRStop == no_){


	    this->GPUFinalize(one, two,
			      hardDevice, iocc);
		
	    // //make sure all the streams are done
	    // pGPU_->synchronize();

	    // //GPU: we will just the master stream to avoid race conditions
	    // cublasSetStream( dv1_handle_[hardDevice],
	    // 		     dv1_streams_[hardDevice][0] );
	
	    // //GPU: accumulate Gamma on device
	    // for(int istream = 1; istream < nstreams_; istream++){
	    // 	cublasDgeam(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N,
	    // 		    nv_, nl_,
	    // 		    &one, dVecPGamma_[hardDevice][0], nv_,
	    // 		    &one, dVecPGamma_[hardDevice][istream], nv_,
	    // 		    dVecPGamma_[hardDevice][0], nv_);
	    // }//istream

	    // //GPU --> HOST
	    // //copy gamma[a1off] from device to host (might as well copy now)
	    // //the synchronization below will make sure the copy finishes
	    // GPUerrchk( cudaMemcpyAsync(pHostGammaLong_[iocc-nf_], 
	    // 			       dVecPGamma_[hardDevice][0], 
	    // 			       nv_*nl_*sizeof(double),
	    // 			       cudaMemcpyDeviceToHost,
	    // 			       dv1_streams_[hardDevice][0]) );

	    // // //unified memory implementation
	    // // //compute active-active block of P in unified memory
	    // // cublasDgemm(dv1_handle_[hardDevice],CUBLAS_OP_T,CUBLAS_OP_N,
	    // // 		no_-nf_, no_-nf_, nv_*nv_,
	    // // 		&minusOne, &pUnifiedStorage_[(no_-nf_)*nv_*nv_], nv_*nv_,
	    // // 		pUnifiedStorage_, nv_*nv_,
	    // // 		&one, pUnifiedPmatOcc_, no_-nf_);
		
	    
	    // //GPU: build part of Gammes^RS
	    // cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_T, CUBLAS_OP_N, 
	    // 		nl_, nl_, nv_,
	    // 		&two, dVecPGamma_[hardDevice][0], nv_,
	    // 		dVecTrans_[hardDevice], nv_,
	    // 		&one, dVecGammaRS_[hardDevice], nl_);

	    // //CPU: compute Paa block
	    // cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
	    // 	    no_-nf_, no_-nf_, nv_*nv_,
	    // 	    -1.0, hVecPTijLong_[0], nv_*nv_,
	    // 	    hVecPEriLong_[0], nv_*nv_,
	    // 	    1.0, pHostPmatOcc_, no_-nf_);
	    
	}//(A2_BRStop == no_)

	pGPU_->synchronize();
	
    }//iocc
    
  


} //GPUGammaPmatFunctor::operator()
#endif
