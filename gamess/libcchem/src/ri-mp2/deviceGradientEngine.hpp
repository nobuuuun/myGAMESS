


#ifndef SRC_RI_MP2_DEVICEGRADIENTENGINE_HPP
#define SRC_RI_MP2_DEVICEGRADIENTENGINE_HPP







#include "thread.hpp"

#include <cuda_runtime.h>
#include <cuda.h>
#include <device.hpp>




#include <Eigen/Dense>


#include <math.hpp>
extern "C" void deviceGradientPvv(const int nv,
				  const int iocc,
				  const int jocc,
				  const double *ea,
				  const double *ev,
				  double *pPvv,
				  double *pEnergy,
				  double *pTijLong,
				  double *pEriLong,
				  const int cudaBlocks, 
				  const int cudaThreadsPB, 
				  const cudaStream_t stream);


extern "C" void deviceJKCopy(int mu,
			     size_t N,
			     size_t naux,
			     size_t mnSize,
			     double *pQmnSym,
			     double *pAOBatch,
			     int cudaBlocks, 
			     int cudaThreadsPB, 
			     cudaStream_t stream);
		


/*
  @brief GPU: build relaxed density, build Gamma, compute the correlation energy
  @author LBR
  @detail 
  -build active-active/virtual/virtual block of relaxed density matrix
  -build gamma
  -compute the correlation energy
  @param nv number of virtual orbitals
  @param no number of occupied orbitals
  @param nl number of auxiliary functions
  @param nf number fo frozen orbitals
  @param ea active orbital energies
  @param ev virtual orbital energies
  @param pGPU pointer to gpu device object
  @param pPqm1 pointer to metric
  @param pmat relaxed density matrix
  @param pHostGammaRS host pointer to Gamma
  @param corrEnergy correlation energy
  @param pHostGammaLong host pointer to active orbital i component of Gamma  
*/

namespace cchem{
    namespace rimp2_gradient{
	namespace detail{

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
	    
class GPUPmatGammaFunctor : boost::noncopyable {

    size_t nv_,no_,nl_,nf_,na_;
    const double *ea_,*ev_;

    ::cchem::rimp2::device *pGPU_;

    int nstreams_;
    int ndevice_;
    int cudaBlocks_;
    int cudaThreadsPB_;
    
    std::vector< cudaStream_t * > dv1_streams_;
    std::vector<cublasHandle_t> dv1_handle_;


    std::vector< std::vector<double*> > dVecPEri_;
    std::vector< std::vector<double*> > dVecPGamma_;

    std::vector< std::vector<double*> > dVecPTijLong_;
    std::vector< std::vector<double*> > dVecPEriLong_;

    //std::vector< std::vector<double*> > dVecPAC_;
    //std::vector< std::vector<double*> > dVecPBC_;
    std::vector< std::vector<double*> > dVecPPmatVirt_;
    
    std::vector< std::vector<double*> > dVecPCiaQ_;
    
    std::vector< std::vector<double*> > dVecPEnergy_;
    
    std::vector< double* > dVecPiaQ_;
    
    std::vector< double* > dVecPPqm1_;
    
    std::vector< double* > dVecEA_;
    
    std::vector< double* > dVecEV_;
    
    std::vector< double* > dVecTrans_;
    
    std::vector< double* > dVecGammaRS_;

    double *pPqm1_;
    MapMatrixXd &pmat_;
    double *pHostGammaRS_;
    double &corrEnergy_;
    std::vector<double *> &pHostGammaLong_;

    //local (big!)
    double *pHostStorage_;
    std::vector<double *> hVecPEriLong_;
    std::vector<double *> hVecPTijLong_;
    
    
    double *pHostPmatOcc_;
    double *pHostPmatVirt_;
    
    // //unified memory implementation
    // double *pUnifiedStorage_;
    // double *pUnifiedPmatOcc_;
    // // std::vector<double*> pUnifiedTijLong_;
    // // std::vector<double*> pUnifiedEriLong_;
    

    // std::vector<cudaEvent_t> dv1_events_;


    void Allocate()
    {
	    
	pHostStorage_ = new double[ 2*(no_-nf_)*nv_*nv_ ];
	for(int i = 0; i < (no_-nf_); i++){
	    hVecPEriLong_.push_back(&pHostStorage_[i*nv_*nv_]);
	    hVecPTijLong_.push_back(&pHostStorage_[i*nv_*nv_ + (no_-nf_)*nv_*nv_ ]);
	}//i
	

	pHostPmatOcc_ = new double[ na_*na_ ] ;
	pHostPmatVirt_ = new double[ nv_*nv_ ] ;
	memset(pHostPmatOcc_, 0, na_*na_*sizeof(double) );
	memset(pHostPmatVirt_, 0, nv_*nv_*sizeof(double) );
	
	double *pTemp = 0;
	
	for (int idevice = 0; idevice < ndevice_; idevice++) {
	    
	    GPUerrchk ( cudaSetDevice( idevice ) );
	    
	    //inverse coulomb metric
	    GPUerrchk( cudaMalloc((void **)&pTemp, nl_*nl_ * sizeof(double)) );
	    dVecPPqm1_[idevice] = pTemp;
	    
	    //3-center ERIs (iocc)
	    GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nl_ * sizeof(double)) );
	    dVecPiaQ_[idevice] = pTemp;
	    
	    //orbitals energies
	    GPUerrchk( cudaMalloc((void **)&pTemp, nv_ * sizeof(double)) );
	    dVecEV_[idevice] = pTemp;
	    GPUerrchk( cudaMemcpy(dVecEV_[idevice], &ev_[0] , (nv_)*sizeof(double),
				  cudaMemcpyHostToDevice) );    
	    
	    //orbital energies
	    GPUerrchk( cudaMalloc((void **)&pTemp, no_ * sizeof(double)) );
	    dVecEA_[idevice] = pTemp;
	    GPUerrchk( cudaMemcpy(dVecEA_[idevice], &ea_[0] , (no_)*sizeof(double),
				  cudaMemcpyHostToDevice) );    
	    
	    //coefficients (iocc)
	    GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nl_ * sizeof(double)) );
	    dVecTrans_[idevice] = pTemp;
	    
	    //device gammaRS
	    GPUerrchk( cudaMalloc((void **)&pTemp, nl_*nl_ * sizeof(double)) );
	    dVecGammaRS_[idevice] = pTemp;
	    cudaMemset(dVecGammaRS_[idevice],0,nl_*nl_ * sizeof(double));
	    
	    
	    //allocate storage for cudaStreams
	    for(int istream=0; istream<nstreams_; istream++){
		
		//device storage for ri-4center ERIs (each stream creates unique integrals)
		GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nv_ * sizeof(double)) );
		dVecPEri_[idevice][istream] = pTemp;
		
		//gamma
		GPUerrchk( cudaMalloc((void **)&pTemp, nl_*nl_ * sizeof(double)) );
		dVecPGamma_[idevice][istream] = pTemp;

		//amplitudes
		GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nv_ * sizeof(double)) );
		dVecPTijLong_[idevice][istream] = pTemp;

		//modifided ERIs
		GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nv_ * sizeof(double)) );
		dVecPEriLong_[idevice][istream] = pTemp;
		
		//intermediates
		// GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nv_ * sizeof(double)) );
		// dVecPAC_[idevice][istream] = pTemp;
		
		//intermediates
		//GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nv_ * sizeof(double)) );
		//dVecPBC_[idevice][istream] = pTemp;
		
		//P(v-v)
		GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nv_ * sizeof(double)) );
		dVecPPmatVirt_[idevice][istream] = pTemp;
		cudaMemset(dVecPPmatVirt_[idevice][istream],0,nv_*nv_ * sizeof(double));
		
		//coefficients (jocc)
		GPUerrchk( cudaMalloc((void **)&pTemp, nv_*nl_ * sizeof(double)) );
		dVecPCiaQ_[idevice][istream] = pTemp;
		
		//for energy evaluation
		GPUerrchk( cudaMalloc((void **)&pTemp, cudaBlocks_*cudaThreadsPB_*sizeof(double)) );
		dVecPEnergy_[idevice][istream] = pTemp;
		cudaMemset(dVecPEnergy_[idevice][istream],0, cudaBlocks_*cudaThreadsPB_*sizeof(double));
		

		}//istream

	    }//idevice



	    // size_t memory;
	    // size_t total;
	    // GPUerrchk( cudaMemGetInfo(&memory,  &total ));
	    // std::cout << "    Free Memory : " << memory << std::flush << std::endl;

	    // //unified memory implementation
	    // //this is place here to make sure there is enough device memory for allocations above
	    // //  (i.e. the unified memory uses what is left over on the device)
	    // //highly experimental
	    // //	    cudaMallocManaged(&ptr, len);
	    // //GPUerrchk( cudaMalloc((void **)&pUnifiedStorage_, 2*na*nv_*nv_ * sizeof(double)) );
	    // GPUerrchk( cudaMallocManaged((void **)&pUnifiedStorage_, 2*na*nv_*nv_ * sizeof(double)) );

	    // // for(int i = 0; i < na; i++){
	    // // 	pUnifiedEriLong_.push_back(&pUnifiedStorage_[i*nv_*nv_]);
	    // // 	pUnifiedTijLong_.push_back(&pUnifiedStorage_[i*nv_*nv_ + na*nv_*nv_ ]);
	    // // }//i
	    // //GPUerrchk( cudaMalloc((void **)&pUnifiedPmatOcc_, na*na * sizeof(double)) );
	    // GPUerrchk( cudaMallocManaged((void **)&pUnifiedPmatOcc_, na*na * sizeof(double)) );
	    // cudaMemset(pUnifiedPmatOcc_, 0, na*na*sizeof(double)); 	    

	}//Allocate


	void Deallocate()
	{

	    delete [] pHostStorage_;


	    
	    // //unified memory implementation
	    // GPUerrchk( cudaFree(pUnifiedStorage_ ) );
	    // GPUerrchk( cudaFree(pUnifiedPmatOcc_ ) );

	    for (int idevice = 0; idevice < ndevice_; idevice++) {

		GPUerrchk( cudaSetDevice( idevice ) );

		GPUerrchk( cudaFree(dVecPPqm1_[idevice]) );
		GPUerrchk( cudaFree(dVecPiaQ_[idevice]) );
		GPUerrchk( cudaFree(dVecEV_[idevice]) );
		GPUerrchk( cudaFree(dVecEA_[idevice]) );
		GPUerrchk( cudaFree(dVecTrans_[idevice]) );
		GPUerrchk( cudaFree(dVecGammaRS_[idevice]) );
		
		//allocate storage for cudaStreams
		for(int istream=0; istream<nstreams_; istream++){
		    
		    GPUerrchk( cudaFree(dVecPEri_[idevice][istream]) );
		    GPUerrchk( cudaFree(dVecPGamma_[idevice][istream]) );
		    GPUerrchk( cudaFree(dVecPTijLong_[idevice][istream]) );
		    GPUerrchk( cudaFree(dVecPEriLong_[idevice][istream]) );
		    //GPUerrchk( cudaFree(dVecPAC_[idevice][istream]) );
		    //GPUerrchk( cudaFree(dVecPBC_[idevice][istream]) );
		    GPUerrchk( cudaFree(dVecPPmatVirt_[idevice][istream]) );
		    GPUerrchk( cudaFree(dVecPCiaQ_[idevice][istream]) );
		    GPUerrchk( cudaFree(dVecPEnergy_[idevice][istream]) );
		    
		}//istream

	    }//idevice


	}//Deallocate

    public:

    static size_t GPUMem(const size_t nl,
			 const size_t nv,
			 const size_t no,
			 ::cchem::rimp2::device *pGPU){
	
	//a mem limit function should not be needed since unified memory call is made
	//   after the the device allocations in "Allocate()"
	//this is just to indicate how much memory is used by the device
	// **does not include unified memory**
	
	
	size_t reqMem;
	//	    GPUerrchk( cudaMemGetInfo(&maxMemSize,  &total ));
	reqMem=0; 
	
	
	const int nstreams = pGPU->getNumStreams();
	const size_t cudaBlocks = pGPU->getCudaBlocks();
	const size_t cudaThreadsPB = pGPU->getCudaThreadsPB();

	/////////////////////////////////////////////
	// minimum memory requirements per device
	/////////////////////////////////////////////

	//////////////////////////////////
	// per device (1 copy)
	//////////////////////////////////
	//occupied orbital energies
	reqMem += no*sizeof(double);
	//singly occupied and virtual orbital energies
	reqMem += nv*sizeof(double);
	//for pPqm1_
	reqMem += nl*nl*sizeof(double);
	//for pVecPiaQ_
	reqMem += nv*nl*sizeof(double);
	//for pVecTrans_
	reqMem += nv*nl*sizeof(double);
	//for pVecGammaRS_
	reqMem += nl*nl*sizeof(double);

	//////////////////////////////////
	//per stream (nstreams copies)
	//////////////////////////////////
	//for dVecPEri_
	reqMem += nv*nv*nstreams*sizeof(double);
	//for dVecPGamma_ 
	reqMem += nl*nl*nstreams*sizeof(double);
	//dVecPTijLong_
	reqMem += nv*nv*nstreams*sizeof(double);
	//VecPEriLong_
	reqMem += nv*nv*nstreams*sizeof(double);
	// //for dVecPAC_
	// reqMem += nv*nv*nstreams*sizeof(double);
	////for dVecPBC_
	//reqMem += nv*nv*nstreams*sizeof(double);
	//for dVecPPmatVirt_
	reqMem += nv*nv*nstreams*sizeof(double);
	//for dVecPCiaQ_
	reqMem += nv*nl*nstreams*sizeof(double);
	//for dVecPEnergy_ (512 is hardwired value 64*8)
	reqMem += cudaBlocks*cudaThreadsPB*nstreams*sizeof(double);
	
	return reqMem;
	
    }
    

    GPUPmatGammaFunctor(size_t nv, size_t no, size_t nl, size_t nf,
			const double *ea, const double *ev,
			::cchem::rimp2::device *pGPU,
			double *pPqm1,
			MapMatrixXd &pmat,
			double *pHostGammaRS,
			double &corrEnergy,
			std::vector<double *> &pHostGammaLong)
	: nv_(nv),no_(no),nl_(nl),nf_(nf),na_(no_-nf_),
	  ea_(ea),ev_(ev),
	  pGPU_(pGPU),
	  nstreams_(pGPU_->getNumStreams()), ndevice_(pGPU_->getNumDevices()),
	  cudaBlocks_(pGPU_->getCudaBlocks()),cudaThreadsPB_(pGPU_->getCudaThreadsPB()),
	  dv1_streams_(pGPU_->cudaStreams() ),dv1_handle_(pGPU_->cublasHandles()),
	  dVecPEri_(std::vector< std::vector<double*> >(ndevice_, std::vector<double*>(nstreams_))),
	  dVecPGamma_(std::vector< std::vector<double*> >(ndevice_, std::vector<double*>(nstreams_))),
	  dVecPTijLong_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  dVecPEriLong_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  //dVecPAC_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  //dVecPBC_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  dVecPPmatVirt_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  dVecPCiaQ_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  dVecPEnergy_(std::vector< std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	  dVecPiaQ_(std::vector< double* > (ndevice_) ),
	  dVecPPqm1_(std::vector< double* > (ndevice_) ),
	  dVecEA_(std::vector< double* > (ndevice_) ),
	  dVecEV_(std::vector< double* > (ndevice_) ),
	  dVecTrans_(std::vector< double* > (ndevice_) ),
	  dVecGammaRS_(std::vector< double* > (ndevice_) ),
	  // dv1_events_(std::vector<cudaEvent_t> (nstreams_) ),
	  pPqm1_(pPqm1),
	  pmat_(pmat),
	  corrEnergy_(corrEnergy),
	  pHostGammaRS_(pHostGammaRS),
	  pHostGammaLong_(pHostGammaLong)
	  
	  //pUnifiedTijLong_(std::vector< double* > (ndevice_) ),
	  //pUnifiedEriLong_(std::vector< double* > (ndevice_) )

    {
	this->Allocate();
	
	//copy the coublomb metric to each device
	for(int idevice = 0; idevice < ndevice_; idevice++){
	    cudaMemcpy(dVecPPqm1_[idevice],pPqm1_, nl_*nl_*sizeof(double), cudaMemcpyHostToDevice);
	}//idevice
	corrEnergy_ = 0.0;
	
	// //create handles and events
	// for(int istream = 0; istream < nstreams_; istream++){
	// 	cudaEventCreate( &dv1_events_[istream] );
	// }//istream
	
    }
    
    
    ~GPUPmatGammaFunctor()
    {
	
	double one = 1.0;
	double zero = 0.0;
	int hardDevice = 0;
	for(int istream = 1; istream < nstreams_; istream++){
	    
	    //GPU: accumulate Pvv block
	    cublasDgeam(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N,
			nv_, nv_,
			&one, dVecPPmatVirt_[hardDevice][0], nv_,
			&one, dVecPPmatVirt_[hardDevice][istream], nv_, 
			dVecPPmatVirt_[hardDevice][0], nv_);
	    
	    //GPU: accumulate energy
	    cublasDgeam(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N,
			cudaBlocks_*cudaThreadsPB_, 1,
			&one, dVecPEnergy_[hardDevice][0], cudaBlocks_*cudaThreadsPB_,
			&one, dVecPEnergy_[hardDevice][istream], cudaBlocks_*cudaThreadsPB_, 
			dVecPEnergy_[hardDevice][0], cudaBlocks_*cudaThreadsPB_);
	    
	}//istream
	
	GPUerrchk( cudaMemcpy(pHostGammaRS_, 
			      dVecGammaRS_[0],
			      nl_*nl_*sizeof(double),
			      cudaMemcpyDeviceToHost) );
	
	
	
	//GPU --> HOST: copy Pvv block
	GPUerrchk( cudaMemcpy(pHostPmatVirt_, 
			      dVecPPmatVirt_[hardDevice][0], 
			      nv_*nv_*sizeof(double),
			      cudaMemcpyDeviceToHost) );
	
	//accumlate Pvv block
	MapMatrixXd pmatVirt(pHostPmatVirt_,nv_,nv_);
	pmat_.block(no_,no_,nv_,nv_) = pmatVirt;
	
	
	// //unified memory implementation
	// GPUerrchk( cudaMemcpy(pHostPmatOcc_, 
	// 		      pUnifiedPmatOcc_, 
	// 		      (no_-nf_)*(no_-nf_)*sizeof(double),
	// 		      cudaMemcpyDefault) );
	
	
	MapMatrixXd pmatOcc(pHostPmatOcc_,no_-nf_,no_-nf_);
	pmat_.block(nf_,nf_,no_-nf_,no_-nf_) = pmatOcc;
	
	
	
	double *pEnergy = new double[cudaBlocks_*cudaThreadsPB_];
	//GPU --> HOST: copy energy
	GPUerrchk( cudaMemcpy(pEnergy, 
			      dVecPEnergy_[hardDevice][0], 
			      cudaBlocks_*cudaThreadsPB_*sizeof(double),
			      cudaMemcpyDeviceToHost) );
	
	
	
	
	//accumlate energy
	for(int i = 0; i<cudaBlocks_*cudaThreadsPB_;i++){
	    corrEnergy_ += pEnergy[i];
	}//i
	
	//std::cout <<"correlation energy: " << corrEnergy_ << std::endl;
	
	delete [] pEnergy;
	    
	
	// //destroy events and handles
	// for(int istream = 0; istream < nstreams_; istream++){
	// 	cudaEventDestroy( dv1_events_[istream] );
	// }//istream
	
	this->Deallocate();

    }//~GPUPmatGammaFunctor()


    
    void operator()(const double *ptr_a1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
		    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
		    const double *ptr_a2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
		    size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols);


    void GPUInitialize(const size_t iocc, int hardDevice, 
		       const size_t NA1BRows, const double *ptr_a1,
		       const double one, const double zero,
		       const size_t A1_BRStop, const size_t A1_BRStart)
    {

	const size_t a1off = iocc-A1_BRStart;
	const size_t hostChunkStart = a1off * NA1BRows;		
	const size_t nBlocksA1 = A1_BRStop-A1_BRStart; 

	//HOST --> GPU
	GPUerrchk( cudaMemcpy2DAsync(dVecPiaQ_[hardDevice], 
				     NA1BRows*sizeof(double),
				     &ptr_a1[ hostChunkStart ],
				     nBlocksA1*NA1BRows*sizeof(double),
				     NA1BRows*sizeof(double), //width
				     nl_,      //height
				     cudaMemcpyHostToDevice) );
		
	//GPU: transform ptr_a1 --> CiaQ
	cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
		    nv_, nl_, nl_, 
		    &one, dVecPiaQ_[hardDevice], nv_, 
		    dVecPPqm1_[hardDevice],nl_,
		    &zero, dVecTrans_[hardDevice], nv_);
	
	//GPU: initialize gamma
	for(int istream = 0; istream < nstreams_; istream++){
	    GPUerrchk( cudaMemsetAsync(dVecPGamma_[hardDevice][istream],
				       0, 
				       nv_*nl_*sizeof(double),
				       dv1_streams_[hardDevice][istream]) );
	}//istream
    }//GPUInitialize


    void GPUWork(const size_t jocc, const size_t A2_BRStart, const size_t NA2BRows,
		 const size_t A2_BRStop, const int hardDevice,
		 const size_t NA1BRows, const double one, const double zero,
		 const size_t iocc, const double* ptr_a2 )
    {


	    const size_t a2off = jocc-A2_BRStart;
	    const size_t hostCoeffStart = a2off * NA2BRows;
	    const size_t nBlocksA2 = A2_BRStop-A2_BRStart;
	    const int istream = jocc%nstreams_;
	    
	    //GPU:
	    cublasSetStream( dv1_handle_[hardDevice],
			     dv1_streams_[hardDevice][istream] );
	    
	    //HOST --> GPU
	    GPUerrchk( cudaMemcpy2DAsync(dVecPCiaQ_[hardDevice][istream], 
					 NA1BRows*sizeof(double),
					 &ptr_a2[ hostCoeffStart ],
					 nBlocksA2*NA1BRows*sizeof(double),
					 NA1BRows*sizeof(double), //width
					 nl_,      //height
					 cudaMemcpyHostToDevice,
					 dv1_streams_[hardDevice][istream]) );
	    
	    //GPU: build ERIs
	    cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_T, 
			nv_, nv_, nl_, 
			&one, dVecPiaQ_[hardDevice], nv_, 
			dVecPCiaQ_[hardDevice][istream],nv_,
			&zero, dVecPEri_[hardDevice][istream], nv_);
	    
	    
	    //GPU: build intermediates
	    deviceGradientPvv(nv_, 
	    		      iocc, 
	    		      jocc, 
	    		      dVecEA_[hardDevice],
	    		      dVecEV_[hardDevice],
	    		      //dVecPAC_[hardDevice][istream],
	    		      //dVecPBC_[hardDevice][istream],
	    		      dVecPEri_[hardDevice][istream],
			      
	    		      dVecPEnergy_[hardDevice][istream],
	    		      dVecPTijLong_[hardDevice][istream],
	    		      dVecPEriLong_[hardDevice][istream],
			      
	    		      cudaBlocks_,
	    		      cudaThreadsPB_,
	    		      dv1_streams_[hardDevice][istream] );
	    
	    // //unified memory implementation
	    // //build intermediates
	    // deviceGradientPvv(nv_, 
	    // 		      iocc, 
	    // 		      jocc, 
	    // 		      dVecEA_[hardDevice],
	    // 		      dVecEV_[hardDevice],
	    // 		      dVecPAC_[hardDevice][istream],
	    // 		      dVecPBC_[hardDevice][istream],
	    // 		      dVecPEri_[hardDevice][istream],
	    
	    // 		      dVecPEnergy_[hardDevice][istream],
	    // 		      &pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ + (no_-nf_)*nv_*nv_ ],
	    // 		      &pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ ],
	    
	    // 		      cudaBlocks_,
	    // 		      cudaThreadsPB_,
	    // 		      dv1_streams_[hardDevice][istream] );

	    //GPU: accumlate P(vv) here
	    cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_T, 
			nv_, nv_, nv_, 
			//&one, dVecPAC_[hardDevice][istream], nv_, 
			&one, dVecPTijLong_[hardDevice][istream], nv_, 
			//dVecPBC_[hardDevice][istream], nv_,
			dVecPEriLong_[hardDevice][istream], nv_,
			&one, dVecPPmatVirt_[hardDevice][istream], nv_);

	    //GPU: build part of gamma
	    cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
	    		nv_, nl_, nv_, 
	    		&one, dVecPTijLong_[hardDevice][istream], nv_, 
	    		dVecPCiaQ_[hardDevice][istream], nv_,
	    		&one, dVecPGamma_[hardDevice][istream], nv_);

	    // //unified memory implementation
	    // //build part of gamma
	    // cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
	    // 		nv_, nl_, nv_, 
	    // 		&one, &pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ + (no_-nf_)*nv_*nv_ ], nv_, 
	    // 		dVecPCiaQ_[hardDevice][istream], nv_,
	    // 		&one, dVecPGamma_[hardDevice][istream], nv_);


	    //GPU --> HOST: copy EriLong
	    GPUerrchk( cudaMemcpyAsync(hVecPEriLong_[jocc-nf_], 
	    			       dVecPEriLong_[hardDevice][istream], 
	    			       nv_*nv_*sizeof(double),
	    			       cudaMemcpyDeviceToHost,
	    			       dv1_streams_[hardDevice][istream]) );

	    //GPU --> HOST: copy TijLong
	    GPUerrchk( cudaMemcpyAsync(hVecPTijLong_[jocc-nf_], 
	    			       dVecPTijLong_[hardDevice][istream], 
	    			       nv_*nv_*sizeof(double),
	    			       cudaMemcpyDeviceToHost,
	    			       dv1_streams_[hardDevice][istream]) );




	    // //experiment
	    // //copy EriLong form device to unified
 	    // GPUerrchk( cudaMemcpyAsync(&pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ ], //pUnifiedEriLong_[jocc-nf_], 
	    // 			  dVecPEriLong_[hardDevice][istream], 
	    // 			  nv_*nv_*sizeof(double),
	    // 			  cudaMemcpyDefault,
	    // 			  dv1_streams_[hardDevice][istream]) );

	    // //experiment
	    // //copy TijLong form device to unified
	    // GPUerrchk( cudaMemcpyAsync(&pUnifiedStorage_[ (jocc-nf_)*nv_*nv_ + (no_-nf_)*nv_*nv_ ], 
	    // 			       dVecPTijLong_[hardDevice][istream], 
	    // 			       nv_*nv_*sizeof(double),
	    // 			       cudaMemcpyDefault,
	    // 			       dv1_streams_[hardDevice][istream]) );
    }//GPUWork
    

    void GPUFinalize(const double one, const double two, 
		const int hardDevice, const size_t iocc)
    {

	    //make sure all the streams are done writing to their dVecPGamma_
	    pGPU_->synchronize();

	    //GPU: we will just the master stream to avoid race conditions
	    cublasSetStream( dv1_handle_[hardDevice],
			     dv1_streams_[hardDevice][0] );
	
	    //GPU: accumulate Gamma on device
	    for(int istream = 1; istream < nstreams_; istream++){
		cublasDgeam(dv1_handle_[hardDevice], CUBLAS_OP_N, CUBLAS_OP_N,
			    nv_, nl_,
			    &one, dVecPGamma_[hardDevice][0], nv_,
			    &one, dVecPGamma_[hardDevice][istream], nv_,
			    dVecPGamma_[hardDevice][0], nv_);
	    }//istream

	    //GPU --> HOST
	    //copy gamma[a1off] from device to host (might as well copy now)
	    //the synchronization below will make sure the copy finishes
	    GPUerrchk( cudaMemcpyAsync(pHostGammaLong_[iocc-nf_], 
				       dVecPGamma_[hardDevice][0], 
				       nv_*nl_*sizeof(double),
				       cudaMemcpyDeviceToHost,
				       dv1_streams_[hardDevice][0]) );

	    // //unified memory implementation
	    // //compute active-active block of P in unified memory
	    // cublasDgemm(dv1_handle_[hardDevice],CUBLAS_OP_T,CUBLAS_OP_N,
	    // 		no_-nf_, no_-nf_, nv_*nv_,
	    // 		&minusOne, &pUnifiedStorage_[(no_-nf_)*nv_*nv_], nv_*nv_,
	    // 		pUnifiedStorage_, nv_*nv_,
	    // 		&one, pUnifiedPmatOcc_, no_-nf_);
		
	    
	    //GPU: build part of Gammes^RS
	    cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_T, CUBLAS_OP_N, 
			nl_, nl_, nv_,
			&two, dVecPGamma_[hardDevice][0], nv_,
			dVecTrans_[hardDevice], nv_,
			&one, dVecGammaRS_[hardDevice], nl_);

	    //CPU: compute Paa block
	    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
	    	    no_-nf_, no_-nf_, nv_*nv_,
	    	    -1.0, hVecPTijLong_[0], nv_*nv_,
	    	    hVecPEriLong_[0], nv_*nv_,
	    	    1.0, pHostPmatOcc_, no_-nf_);

    }







};//GPUPmatGammaFunctor




	}//detail
    }//rimp2_gradient
}//cchem

















/*
  @brief JK matrix build
  @author LBR
  @detail this is used to build both the J(CPU) and K(GPU) matricies
  @param mnSize number of unique OBS pairs
  @param N number of OBS functions
  @param no number of occupied orbitals
  @param pJtemp scratch (host) pointer for building J matrix
  @param pD host pointer to AO density
  @param pJ host pointer to J matrix
  @param pDmn host pointer to unique OBS pairs
  @param pDAux scratch host pointer
  @param pK pointer to K matrix
  @param pLeft left decomposition
  @param pRight right decomposition
  @param dims dimension of J,K,pLeft,pRight
  @param maxNAB man number of auxiliary function uses per data read
  @param pGPU host pointer to device object
  @param myRank rank of calling process in global parallel environment
*/

namespace cchem{
    namespace rimp2_gradient{
	namespace detail{

class GPUBuildJKMatrixFunctor : boost::noncopyable {

    const size_t mnSize_, N_, no_;
    double *pJTemp_;
    std::vector<double *> &pD_, &pJ_, &pDmn_, &pDAux_,&pK_,&pLeft_,&pRight_;
    std::vector<int> &dims_;
    size_t maxNAB_;
    size_t maxNocc_;

    ::cchem::rimp2::device *pGPU_;
    
    int ndevice_;
    int nstreams_;
    int cudaBlocks_;
    int cudaThreadsPB_;
    std::vector< cudaStream_t * > dv1_streams_;
    std::vector<cublasHandle_t> dv1_handle_;
    
    
    int nFactors_;
    //device pointers
    std::vector< std::vector<double*> > dVecJ_;
    std::vector< std::vector<double*> > dVecK_;
    std::vector<double*> dVecpExpLeft_;
    std::vector<double*> dVecpExpRight_;
    std::vector< std::vector<double*> > dVecLeft_;
    std::vector< std::vector<double*> > dVecRight_;
    std::vector<double*> dVecPQmnSym_;
    std::vector< std::vector<double*> > dVecAOBatch_;
    std::vector<cudaEvent_t> dv1_events_;
    std::vector<cublasHandle_t> cublasStreamHandles_;
    int myRank_;
    
    std::vector<double*>hVecAOBatch_;
    std::vector<double*>hVecpExpLeft_;
    std::vector<double*>hVecpExpRight_;

public:
    
    //this function estimates the IOBufferSize for this functor
    static size_t GPUMemLimit(const size_t nl,
			      std::vector<int> &dims,
			      const size_t N,
			      const size_t mnSize,
			      const int nstreams,
			      const int nfactors){
	
	size_t maxMemSize, total, currentMem, memSet;
	GPUerrchk( cudaMemGetInfo(&maxMemSize,  &total ));
	//take 100 MB for gpu? maybe not needed?
	currentMem=100*(1000*1000);
	
	size_t maxNocc = std::max(dims[0],dims[1]);
	currentMem += nfactors*N*N*sizeof(double); //J
	currentMem += nfactors*N*N*sizeof(double); //K
	currentMem += nfactors*dims[0]*N*sizeof(double); //vec Left/Right
	currentMem += nfactors*dims[1]*N*sizeof(double); //vec Left/Right
	
	size_t maxAuxIO;
	for(size_t i = 0; i < nl; i++){
	    
	    currentMem += maxNocc*N*sizeof(double); //pExpRight
	    currentMem += maxNocc*N*sizeof(double); //pExpLeft
	    
	    currentMem += mnSize*sizeof(double); //vecPQmnSym
	    
	    currentMem += N*nstreams*sizeof(double); //dvecAOBatch
	    
	    if(currentMem > maxMemSize)break;
	    
	    //memSet = currentMem;
	    maxAuxIO = i+1;
	}
	
	//std::cout << "max aux for IO is found as: " << maxAuxIO << std::endl;
	
	
	memSet = maxAuxIO*mnSize*2*sizeof(double);
	
	//std::cout << maxMemSize << " " << memSet << std::endl;
	
	return memSet/1000/1000 +1;
	
    }//IOMemLimit
    
    
    void allocate(){
	
	double *pTemp = 0;
	
	for (int ithread = 1; ithread < omp_get_max_threads(); ithread++){
	    hVecAOBatch_[ithread] = new double[maxNAB_*N_];
	    hVecpExpLeft_[ithread] = new double[maxNAB_*maxNocc_];
	    hVecpExpRight_[ithread] = new double[maxNAB_*maxNocc_];
	}//ithread
	
	for (int idevice = 0; idevice < ndevice_; idevice++) {
	    
	    GPUerrchk ( cudaSetDevice( idevice ) );
	    
	    GPUerrchk( cudaMalloc((void **)&pTemp, maxNAB_*maxNocc_*N_ * sizeof(double)) );
	    dVecpExpLeft_[idevice] = pTemp;

	    GPUerrchk( cudaMalloc((void **)&pTemp, maxNAB_*maxNocc_*N_ * sizeof(double)) );
	    dVecpExpRight_[idevice] = pTemp;

	    GPUerrchk( cudaMalloc((void **)&pTemp, mnSize_*maxNAB_* sizeof(double)) );
	    dVecPQmnSym_[idevice] = pTemp;
	    
	    for(int ifactor = 0; ifactor < nFactors_; ifactor++){
		
		GPUerrchk( cudaMalloc((void **)&pTemp, N_*N_ * sizeof(double)) );
		dVecJ_[idevice][ifactor]  = pTemp;
		cudaMemset(dVecJ_[idevice][ifactor] ,0, N_*N_ * sizeof(double));
		
		GPUerrchk( cudaMalloc((void **)&pTemp, N_*N_ * sizeof(double)) );
		dVecK_[idevice][ifactor]  = pTemp;
		cudaMemset(dVecK_[idevice][ifactor] ,0, N_*N_ * sizeof(double));
		
		GPUerrchk( cudaMalloc((void **)&pTemp, dims_[ifactor]*N_* sizeof(double)) );
		dVecLeft_[idevice][ifactor]  = pTemp;
		
		GPUerrchk( cudaMalloc((void **)&pTemp, dims_[ifactor]*N_* sizeof(double)) );
		dVecRight_[idevice][ifactor] = pTemp;

	    }//ifactor

	    for (int istream = 0; istream < nstreams_; istream++) {
		
		GPUerrchk( cudaMalloc((void **)&pTemp, N_*maxNAB_* sizeof(double)) );
		dVecAOBatch_[idevice][istream] = pTemp;
		
	    }//istream


	}//idevice



    }//allocate




    void deallocate(){

	
	for (int ithread = 1; ithread < omp_get_max_threads(); ithread++){
	    delete [] hVecAOBatch_[ithread];
	    delete [] hVecpExpLeft_[ithread];
	    delete [] hVecpExpRight_[ithread];
	}//ithread
	
	
	for (int idevice = 0; idevice < ndevice_; idevice++) {
	    
	    GPUerrchk ( cudaSetDevice( idevice ) );

	    GPUerrchk( cudaFree(dVecpExpLeft_[idevice]) );
	    GPUerrchk( cudaFree(dVecpExpRight_[idevice]) );
	    GPUerrchk( cudaFree(dVecPQmnSym_[idevice]) );
	    
	    for(int ifactor = 0; ifactor < nFactors_; ifactor++){
		GPUerrchk( cudaFree(dVecJ_[idevice][ifactor]) );
		GPUerrchk( cudaFree(dVecK_[idevice][ifactor]) );
		GPUerrchk( cudaFree(dVecLeft_[idevice][ifactor]) );
		GPUerrchk( cudaFree(dVecRight_[idevice][ifactor]) );
	    }//ifactor
	    
	    for (int istream = 0; istream < nstreams_; istream++) {
		
		GPUerrchk( cudaFree(dVecAOBatch_[idevice][istream]) );
		
	    }//istream
	    
	}//idevice

    }//deallocate

    
    
    
    GPUBuildJKMatrixFunctor(const size_t mnSize, 
			    const size_t N,
			    const size_t no,
			    double *pJTemp, 
			    std::vector<double*> &pD, 
			    std::vector<double*> &pJ, 
			    std::vector<double*> &pDmn, 
			    std::vector<double*> &pDAux,
			    std::vector<double*> &pK, 
			    std::vector<double*> &pLeft, 
			    std::vector<double*> &pRight, 
			    std::vector<int> &dims, 
			    size_t maxNAB,
			    ::cchem::rimp2::device *pGPU,
			    int myRank):
	mnSize_(mnSize), 
	N_(N), 
	no_(no),
	pJTemp_(pJTemp), 
	pD_(pD), 
	pJ_(pJ), 
	pDmn_(pDmn), 
	pDAux_(pDAux),
	pK_(pK),
	pLeft_(pLeft),
	pRight_(pRight),
	dims_(dims),
	maxNAB_(maxNAB),
	maxNocc_(std::max(dims_[0],dims_[1])),
	pGPU_(pGPU),
	ndevice_( pGPU->getNumDevices() ),
	nstreams_( pGPU->getNumStreams() ),
	cudaBlocks_( pGPU->getCudaBlocks() ),
	cudaThreadsPB_( pGPU->getCudaThreadsPB() ),
	dv1_streams_( pGPU->cudaStreams() ),
	dv1_handle_( pGPU->cublasHandles() ),
	nFactors_(pK_.size()),
	dVecJ_(std::vector<std::vector< double* > >(ndevice_, std::vector<double*>(nFactors_) ) ),
	dVecK_(std::vector<std::vector< double* > >(ndevice_, std::vector<double*>(nFactors_) ) ),
	dVecpExpLeft_(std::vector< double* > (ndevice_) ),
	dVecpExpRight_(std::vector< double* > (ndevice_) ),
	dVecLeft_(std::vector<std::vector<double*> >(ndevice_,std::vector<double*>(nFactors_) ) ),
	dVecRight_(std::vector<std::vector<double*> >(ndevice_,std::vector<double*>(nFactors_) )),
	dVecPQmnSym_(std::vector< double* > (ndevice_) ),
	dVecAOBatch_(std::vector<std::vector<double*> >(ndevice_,std::vector<double*>(nstreams_))),
	dv1_events_(std::vector<cudaEvent_t> (nstreams_) ),
	cublasStreamHandles_(std::vector< cublasHandle_t > (nstreams_) ),
	myRank_(myRank),
	
	hVecAOBatch_ (std::vector<double*> (omp_get_max_threads()) ),
	hVecpExpLeft_ (std::vector<double*> (omp_get_max_threads()) ),
	hVecpExpRight_ (std::vector<double*> (omp_get_max_threads()) )
	
    {
	
	//create handles and events
	for(int istream = 0; istream < nstreams_; istream++){
	    cublasCreate( &cublasStreamHandles_[istream]);
	    cudaEventCreate( &dv1_events_[istream] );
	}//istream
	
	
	this->allocate();
	
	this->setLeftData();
	this->setRightData();
	
    };
    
    void setLeftData(){
	
	for(int idevice = 0; idevice < ndevice_; idevice++){
	    for(int ifactor = 0; ifactor < nFactors_; ifactor++){
	    	GPUerrchk( cudaMemcpy(dVecLeft_[idevice][ifactor],pLeft_[ifactor], 
	    			      dims_[ifactor]*N_*sizeof(double),
	    			      cudaMemcpyHostToDevice) );   
	    } //ifactor
	}//idevice
	
    }//setLeftData

    void setRightData(){
	
	
	for(int idevice = 0; idevice < ndevice_; idevice++){
	    for(int ifactor = 0; ifactor < nFactors_; ifactor++){
	    	GPUerrchk( cudaMemcpy(dVecRight_[idevice][ifactor],pRight_[ifactor], 
	    			      dims_[ifactor]*N_*sizeof(double),
	    			      cudaMemcpyHostToDevice) );  
	    }//ifactor
	}//idevice
	
    }//setRightData

    ~GPUBuildJKMatrixFunctor(){
	
	for(int ifactor = 0; ifactor < nFactors_; ifactor++){
	    //need to do a reduction!!!! (with more than one device)
	    if(ndevice_ > 1)std::cout << "you have a problem "<< __FILE__ << ":" << __LINE__ << std::endl;
	    GPUerrchk( cudaMemcpy(pK_[ifactor], 
				  dVecK_[0][ifactor], 
				  N_*N_*sizeof(double),
				  cudaMemcpyDeviceToHost) );
	}//ifactor
	
	//destroy events and handles
	for(int istream = 0; istream < nstreams_; istream++){
	    cublasDestroy( cublasStreamHandles_[istream]);
	    cudaEventDestroy( dv1_events_[istream] );
	}//istream
	
	this->deallocate();
	
    };
    
    void operator()(double *pQmnSym, 
		    size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
		    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols)
    {
	
	int tid = omp_get_thread_num();
	
	const size_t naux = A1_BCStop - A1_BCStart;
	
	utility::timer timerTotal;
	utility::timer timerTest;
	if(tid == 0){
	    timerTest.reset();
	    timerTotal.reset();
	}
	
	double one = 1.0;
	double zero = 0.0;
	
	int hardDevice = 0;
	cudaSetDevice( hardDevice );
	
	cudaEvent_t magicalEvent;
	if(tid == 0)cudaEventCreate( &magicalEvent );
	
	if(tid == 0)
	    for(int idevice = 0; idevice < ndevice_; idevice++)
		GPUerrchk( cudaMemcpyAsync(dVecPQmnSym_[idevice], 
					   pQmnSym , mnSize_*naux*sizeof(double),
					   cudaMemcpyHostToDevice) );    
	
	/////////////////////////////
	// build J matrix
	/////////////////////////////
	
#pragma omp single nowait
	for(int iF = 0; iF < pJ_.size(); iF++){
	    
	    //get factored density matrix (if applicable)
	    double *D = pD_[iF];
	    
	    //create lower triangular density matrix
	    for(size_t m = 0, mn = 0; m < N_; m++){
		for(size_t n = 0; n <=m ; n++, mn++){
		    pDmn_[iF][mn] = ( m==n ? D[m+n*N_] : D[m+n*N_] + D[n+m*N_] );
		}
	    }
	    
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
	    for(size_t m = 0, mn = 0; m < N_; m++){
		for(size_t n = 0; n <=m ; n++, mn++){
		    pJ_[iF][m + N_*n ] += pJTemp_[mn];
		    if( m!=n ) pJ_[iF][n + m*N_] += pJTemp_[mn];
		}//n
	    }//m
	    
	}//iF
	
	/////////////////////////////
	// build K matrix
	/////////////////////////////
	
	//make sure transfer is done
	//pGPU_->synchronize();
	size_t streamCounter = 0;
	if(tid == 0)pGPU_->synchronize();
	
	// if(tid == 0)
	//     std::cout << "time to transfer ints and compute J-matrix " 
	// 	      << timerTest << std::endl;
	
	for(int iF = 0; iF < pK_.size(); iF++){
	    
	    const int nocc = dims_[iF];
	    
#pragma omp for schedule(dynamic)		
	    for (int mu = 0; mu < N_; mu++){
		
		
		if(tid == 0){
		    
		    streamCounter++;
		    
		    cudaEventSynchronize (magicalEvent);
		    
		    int istream  =mu%nstreams_;
		    deviceJKCopy(mu,
				 N_,
				 naux,
				 mnSize_,
				 dVecPQmnSym_[hardDevice],
				 dVecAOBatch_[hardDevice][istream],
				 cudaBlocks_, 
				 cudaThreadsPB_,
				 dv1_streams_[hardDevice][istream]);
		    
		    
		    
		    //cudaStreamSynchronize(dv1_streams_[hardDevice][istream]);
		    
		    
		    //this is kind of wonky, but the stream must wait for teh dviceJKCopy to finish
		    //    I accomplised this with cudaEvents
		    
		    //IF YOU HAVE ISSUE HERE, make sure you set device when creating
		    // events on mutliple devices (IT MATTERS!)
		    cudaEventRecord( dv1_events_[istream], dv1_streams_[hardDevice][istream]);
		    
		    cudaStreamWaitEvent( dv1_streams_[hardDevice][istream],
					 dv1_events_[istream],
					 0) ;
		    
		    cublasSetStream( cublasStreamHandles_[istream],
				     dv1_streams_[hardDevice][istream] );
		    
		    
		    // for(int nu = 0; nu < N_; nu++){
		    // 	const unsigned int start = 
		    // 	    (mu >= nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu );
		    // 	cublasDcopy(cublasStreamHandles_[istream],
		    // 		    (int)naux,
		    // 		    &dVecPQmnSym_[hardDevice][start], (int)mnSize_,
		    // 		    &dVecAOBatch_[hardDevice][istream][nu], (int)N_);
		    // }
		    
		    
		    
		    //build left side
		    cublasDgemm(cublasStreamHandles_[istream], CUBLAS_OP_N, CUBLAS_OP_N, 
				nocc, naux, N_, 
				&one, dVecLeft_[hardDevice][iF], nocc, 
				dVecAOBatch_[hardDevice][istream], N_,
				&zero, &dVecpExpLeft_[hardDevice][mu*nocc*naux], nocc);
		    
		    //build right side
		    cublasDgemm(cublasStreamHandles_[istream], CUBLAS_OP_N, CUBLAS_OP_N, 
				nocc, naux, N_, 
				&one, dVecRight_[hardDevice][iF], nocc, 
				dVecAOBatch_[hardDevice][istream], N_,
				&zero, &dVecpExpRight_[hardDevice][mu*nocc*naux], nocc);
		    
		    if(streamCounter > nstreams_*nstreams_ &&
		       streamCounter%nstreams_ == 0){
			streamCounter = 0;
			cudaEventRecord (magicalEvent);
		    }
		    
		}else{
		    
		    //unpack (mu nu| L) integrals
		    for (int nu = 0; nu < N_; nu++){
			const size_t start = (mu >= nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu );
			cblas_dcopy(naux, &pQmnSym[start], mnSize_, &hVecAOBatch_[tid][nu], N_);
		    }
		    //build left factor of K
		    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				nocc, naux, N_,
				1.0, pLeft_[iF], nocc,
				hVecAOBatch_[tid], N_,
				0.0, hVecpExpLeft_[tid], nocc);
		    
		    
		    GPUerrchk( cudaMemcpy2DAsync(&dVecpExpLeft_[hardDevice][mu*nocc*naux], 
						 naux*sizeof(double),
						 hVecpExpLeft_[tid],
						 naux*sizeof(double),
						 naux*sizeof(double), //width
						 nocc,      //height
						 cudaMemcpyHostToDevice) );
		    
		    //build right factor of K
		    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				nocc, naux, N_,
				1.0, pRight_[iF], nocc,
				hVecAOBatch_[tid], N_,
				0.0, hVecpExpRight_[tid], nocc);
		    
		    
		    GPUerrchk( cudaMemcpy2D(&dVecpExpRight_[hardDevice][mu*nocc*naux], 
					    naux*sizeof(double),
					    hVecpExpRight_[tid],
					    naux*sizeof(double),
					    naux*sizeof(double), //width
					    nocc,      //height
					    cudaMemcpyHostToDevice) );
		    
		}//(tid == 0)
		
	    }//mu
	    
#pragma omp single
	    {
		pGPU_->synchronize();
		//build K
		cublasDgemm(dv1_handle_[hardDevice], CUBLAS_OP_T, CUBLAS_OP_N, 
			    N_, N_, nocc*naux, 
			    &one, dVecpExpLeft_[hardDevice], nocc*naux, 
			    dVecpExpRight_[hardDevice], nocc*naux, 
			    &one, dVecK_[hardDevice][iF], N_);
		pGPU_->synchronize();
	    }
	    
	}//iF
	
	
	if(tid == 0)cudaEventDestroy( magicalEvent );	    
	// if(tid == 0)
	//     std::cout << myRank_ 
	// 	      << " GPU Chunk time " 
	// 	      << timerTotal << std::endl;
	
    };//operator()
    
};//GPUBuildJKMatrixFunctor
	    
	    
	    
	}//detail
    }//rimp2_gradient
}//cchem


#endif //header guard

