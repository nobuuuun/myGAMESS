/*
 * ERIContractionEnergyAccum.hpp
 *
 *  Created on: May 16, 2016
 *      Author: luke
 */

#ifndef SRC_RI_MP2_ERICONTRACTIONENERGYACCUM_HPP_
#define SRC_RI_MP2_ERICONTRACTIONENERGYACCUM_HPP_

#include "god.hpp"

#include <ri-energy-terms.hpp>
#include <math.hpp>
#include <vector>
#include "thread.hpp"

 //lbr
//#undef HAVE_CUBLAS
//#define HAVE_CUBLAS 0


#if !HAVE_CUBLAS


namespace cchem{
namespace rimp2{
namespace detail{

using cchem::Thread;

class RIMP2Engine : boost::noncopyable {


	Thread thread_;

	const size_t &nd_,&ns_,&nv_,&nl_;
	const size_t nvs_;
	const int stride2_;
	const int stride_row_;
	double *data_twoc_;

    const double *ea_,*ev_;
	std::vector<double> &e_m_;
	std::vector<double *> pfock_;

	cchem::rimp2::detail::Energy_Term &energy_term_;

    omp_lock_t *pLocks_;

public:

    void reduceThreadEnergies(double *energyTerm,
			      double * global_pfock, Parallel &pe){

	if(ns_){
#pragma omp critical
	{

		for (size_t iocc = 0; iocc < nd_; iocc++){
		    double *glptr = &global_pfock[iocc*nv_];
		    //		    double *glptr = global_pfock.at(iocc);
		    double *thptr = &(thread_.data[1][iocc*nv_]);
		    for(int a = 0; a < nv_; a++){
						glptr[a] += thptr[a];
		    }
 		}

 	}//omp critical

#pragma omp barrier
	    
#pragma omp single
	{
	    pe.reduce("+",global_pfock, (size_t)(nd_*nv_) );
		     		    // energyTerm[4] = energy_term_.evaluate5(ea_,ev_,global_pfock,nd_,nv_,ns_);

	    MapMatrixXd pfocker(global_pfock, nv_,nd_);
	    if(pe.rank()==0){
		for (size_t d = 0; d < nd_; d++){
		    for(size_t a = 0; a <nv_; a++){
			
			energyTerm[4] +=  pow (pfocker(a,d),2)/(ea_[d]-ev_[ns_+a]);
			
		    }//a
		}//d
		energyTerm[4] = energyTerm[4]/2.0;
	    }
	    
	}//omp single
	}//ns
	
	//reduce all thread contributions to the energy(s)
	energy_term_(energyTerm);


    }//reduceThreadEnergies

	static
	boost::array<size_t,2> memory(const int stride2, const int stride_row, 
				      const size_t nv,const size_t ns,const size_t nd){

	    const size_t nvs = nv+ns;
		boost::array<size_t,2> memory;
		memory[0] = stride2*nvs*nvs*stride_row;  //storage for approximate 4center ERIs
		memory[1] = nd*nv;  //storage for approximate 4center ERIs
		return memory;
	}

	RIMP2Engine(size_t &nd, size_t &ns, size_t &nv, size_t &nl,
			int &stride2, int &stride_row, double *data_twoc,
			const double *ea, const double *ev,
			std::vector<double> &e_m,
		    cchem::rimp2::detail::Energy_Term &energy_term,
		    omp_lock_t *pLocks):
				nd_(nd),ns_(ns),nv_(nv),nl_(nl), nvs_(nv+ns), stride2_(stride2), stride_row_(stride_row),
				data_twoc_(data_twoc), ea_(ea), ev_(ev), e_m_(e_m),
				energy_term_(energy_term),
				pLocks_(pLocks)
    {
	     BOOST_AUTO(memory, this->memory(stride2_,stride_row_,nv_,ns_,nd_));
		for (int i = 0; i < memory.size(); ++i) {
			thread_.data[i] = thread_.malloc(memory[i]);
		}//memory

		double *ptr = thread_.data[1];
		for (size_t i = 0; i < nd_; i++){
		    pfock_.push_back(&ptr[i*nv_]);
		    std::fill(&ptr[i*nv_],&ptr[i*nv_+nv_],0);
		}//i

	 };//RIMP2Engine

    void resetLocks(){

	for(int i = 0; i < omp_get_max_threads(); i++){
	    omp_init_lock(&pLocks_[i]);
	    omp_set_lock(&pLocks_[i]);
	}
 
	omp_unset_lock(&pLocks_[0]);
	
    }//initializeLocks


	void operator()(double *ptr_bj, const int jblock, const int jblock_end,
			const int j_remain){


	    //IF YOU LIKE TO TEST THIS, YOU NEED TO USE 1 OMP THREAD!!!
	    //IF YOU LIKE TO TEST THIS, YOU NEED TO USE 1 OMP THREAD!!!

// 	    double *ptr_bi;

// #pragma omp single
// 	    ptr_bi = new double[j_remain*nl_*nvs_];

//  	    Eigen::MatrixXd pqm1(nl_,nl_);
//  	    MapMatrixXd lm1(data_twoc_,nl_,nl_);

//  #pragma omp master
//  	    std::cout << lm1.block(0,0,5,5) << std::endl << std::endl;

//  	    lm1 = lm1.triangularView<Eigen::Upper>();
//  	    pqm1 = lm1*lm1.transpose();

// 	    // 	    pqm1 = lm1.triangularView<Eigen::Upper>()*lm1.triangularView<Eigen::Upper>().transpose();


// #pragma omp single
// 		memcpy( ptr_bi, ptr_bj, nvs_*nl_*j_remain*sizeof(double) );


// #pragma omp for schedule(dynamic,1)
// 		for(int iocc = jblock; iocc < jblock_end; iocc++){

//this should give you numerically stable result
// 			// //form B_ia^Q = (ia|Q)L-1  -- L-1 is acutally upper (sorry)
// 			// int s = iocc - jblock;
// 			// cblas_dtrmm(CblasColMajor,
// 			// 		CblasRight,
// 			// 		CblasUpper,
// 			// 		CblasNoTrans,
// 			// 		CblasNonUnit,
// 			// 		nvs_, nl_,
// 			// 		1.0,data_twoc_, nl_,
// 			// 		&ptr_bi[s*nvs_], nvs_*j_remain );


//I question the numerical stability of this as the created of pqm1 was done via Lm1*Lm1^T,
//	    and hence the numerical error is compounded.
// 			// //form C_ia^Q = B_ia_Q*(L-1)^T
// 			// //			int s = iocc - jblock;
// 			// cblas_dtrmm(CblasColMajor,
// 			// 		CblasRight,
// 			// 		CblasUpper,
// 			// 		CblasTrans,
// 			// 		CblasNonUnit,
// 			// 		nvs_, nl_,
// 			// 		1.0,data_twoc_, nl_,
// 			// 		&ptr_bi[s*nvs_], nvs_*j_remain );

//                         int s = iocc - jblock; 
// 			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
// 				    nvs_,nl_,nl_,
// 				    1.0, &ptr_bj[s*nvs_],nvs_*j_remain,
// 				    pqm1.data(),nl_,
// 				    0.0, &ptr_bi[s*nvs_], nvs_*j_remain);

// 		}//iocc

		

// #pragma omp for schedule(dynamic,1)
// 		for(int iocc = jblock; iocc < jblock_end; iocc++){

// 			MapMatrixXd eri(thread_.data[0],nvs_,(iocc-jblock+1)*nvs_);

// 			//form ERIs (with Vp ordered array)
// 			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
// 					nvs_,(iocc-jblock+1)*nvs_, nl_,
// 					1.0,&ptr_bi[ (iocc-jblock)*nvs_], nvs_*j_remain,
// 					ptr_bj,nvs_*j_remain,
// 					0.0,thread_.data[0],(nvs_));
// 			//#pragma omp critical

// 			energy_term_(eri, jblock, iocc, ea_, ev_, e_m_, pfock_);

// 		}//iocc

// #pragma omp single
// 		memcpy( ptr_bj, ptr_bi, nvs_*nl_*j_remain*sizeof(double) );


// #pragma omp single
// 		delete [] ptr_bi;



//#pragma omp master
//	    omp_lock_t *lock = new omp_lock_t[ omp_get_num_threads() ];

//#pragma omp barrier















// #pragma omp for schedule(dynamic,1)
// 		for(int iocc = jblock; iocc < jblock_end; iocc++){

// 			//form B_ia^Q = (ia|Q)L-1  -- L-1 is acutally upper (sorry)
// 			int s = iocc - jblock;
// 			cblas_dtrmm(CblasColMajor,
// 					CblasRight,
// 					CblasUpper,
// 					CblasNoTrans,
// 					CblasNonUnit,
// 					nvs_, nl_,
// 					1.0,data_twoc_, nl_,
// 					&ptr_bj[s*nvs_], nvs_*j_remain );
// 		}//iocc

// #pragma omp for schedule(dynamic,1)
// 		for(int iocc = jblock; iocc < jblock_end; iocc++){

// 			MapMatrixXd eri(thread_.data[0],nvs_,(iocc-jblock+1)*nvs_);

// 			//form ERIs (with Vp ordered array)
// 			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
// 					nvs_,(iocc-jblock+1)*nvs_, nl_,
// 					1.0,&ptr_bj[ (iocc-jblock)*nvs_], nvs_*j_remain,
// 					ptr_bj,nvs_*j_remain,
// 					0.0,thread_.data[0],(nvs_));
// 			//#pragma omp critical

// 			energy_term_(eri, jblock, iocc, ea_, ev_, e_m_, pfock_);

// 		}//iocc







//#pragma omp for schedule(static,1)

//	    double *pThread = new double[nvs_*nl_];
#pragma omp for schedule(dynamic,1)
		for(int iocc = jblock; iocc < jblock_end; iocc++){

			//form B_ia^Q = (ia|Q)L-1  -- L-1 is acutally upper (sorry)
			int s = iocc - jblock;
			cblas_dtrmm(CblasColMajor,
					CblasRight,
					CblasUpper,
					CblasNoTrans,
					CblasNonUnit,
					nvs_, nl_,
					1.0,data_twoc_, nl_,
					&ptr_bj[s*nvs_], nvs_*j_remain );


// cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
// 	    nvs_, nl_, nl_,
// 	    1.0,&ptr_bj[s*nvs_], nvs_*j_remain,
// 	    data_twoc_, nl_,
// 	    0.0, pThread,nvs_);
// 	for(int i = 0; i < nl_; i++){
// cblas_dcopy(nvs_, &pThread[nvs_*i] ,1, &ptr_bj[s*nvs_ + nvs_*j_remain*i] ,1);
// 	}


			//must use schedule static to get this to work properly
			// omp_set_lock( &pLocks_[omp_get_thread_num()] );
			// omp_unset_lock( &pLocks_[ (omp_get_thread_num()+1)%omp_get_num_threads() ] );
			//this works with schedule dynamic
			omp_set_lock( &pLocks_[ (iocc-jblock)%omp_get_num_threads() ] );
			omp_unset_lock( &pLocks_[ (iocc-jblock+1)%omp_get_num_threads() ] );


// 		}//iocc

// #pragma omp for schedule(dynamic,1)
// 		for(int iocc = jblock; iocc < jblock_end; iocc++){



			MapMatrixXd eri(thread_.data[0],nvs_,(iocc-jblock+1)*nvs_);

			//form ERIs (with Vp ordered array)
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
					nvs_,(iocc-jblock+1)*nvs_, nl_,
					1.0,&ptr_bj[ (iocc-jblock)*nvs_], nvs_*j_remain,
					ptr_bj,nvs_*j_remain,
					0.0,thread_.data[0],(nvs_));

 			//#pragma omp critical

			//			for (int i = 0; i < 3; i++)std::cout << thread_.data[0][i] <<std::endl;
			//			for (int i = 0; i < 3; i++)std::cout << data_twoc_[i] <<std::endl;
			//			std::cout << std::endl;
			energy_term_(eri, jblock, iocc, ea_, ev_, e_m_, pfock_);

		}//iocc


#pragma omp for schedule(dynamic,1) nowait
		for(int iocc = jblock; iocc < jblock_end; iocc++){

			//form C_ia^Q = B_ia_Q*(L-1)^T
			int s = iocc - jblock;
			cblas_dtrmm(CblasColMajor,
					CblasRight,
					CblasUpper,
					CblasTrans,
					CblasNonUnit,
					nvs_, nl_,
					1.0,data_twoc_, nl_,
					&ptr_bj[s*nvs_], nvs_*j_remain );



// cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
// 	    nvs_, nl_, nl_,
// 	    1.0,&ptr_bj[s*nvs_], nvs_*j_remain,
// 	    data_twoc_, nl_,
// 	    0.0, pThread,nvs_);
// for(int i = 0; i < nl_; i++){
//     cblas_dcopy(nvs_, &pThread[nvs_*i] ,1, &ptr_bj[s*nvs_ + nvs_*j_remain*i] ,1);
// }

		}//iocc

//		delete [] pThread;

		//#pragma omp single
		//		delete [] lock;

	};//operator()


	void operator()(double *ptr_bj, double *ptr_bi, const int num_chunks,
			const int i_remain, const int j_remain,
			const int iblock, const int jblock,
			const int jblock_end){

#pragma omp for schedule(dynamic,1)
		for (int ichunk = 0; ichunk < num_chunks; ichunk++){
			int chunk_len = stride_row_;
			int chunk_pos = ichunk*stride_row_;
			if( (ichunk + 1) == num_chunks) chunk_len = i_remain - (num_chunks-1)*stride_row_;

			MapMatrixXd eri(thread_.data[0],chunk_len*nvs_,j_remain*nvs_);

			//form ERIs
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
					chunk_len*nvs_, stride2_*nvs_, nl_,
					1.0,&ptr_bi[ chunk_pos*nvs_], nvs_*i_remain,
					ptr_bj,nvs_*j_remain,
					0.0,thread_.data[0],nvs_ );

			int ipos = 0;
	for(int iocc = iblock+chunk_pos; iocc < iblock+chunk_pos+chunk_len; iocc++, ipos++){
	    energy_term_(eri,ipos, jblock, jblock_end, iocc, ea_, ev_, e_m_, pfock_);
	} //iocc



		}//ichunk


	};//operator()

}; //RIMP2Engine


}//namespace detail
}//namespace rimp2
}//namespace cchem



#endif //!HAVE_CUBLAS

















#if HAVE_CUBLAS


#include <cuda_runtime.h>
#include <cuda.h>
#include <device.hpp>

extern "C" void deviceCudaEnergyIJ(int iblock, int jblock, int jblock_end,
				   const double *d_C,
				   const double *d_ea,
				   const double *d_ev,
				   const double* d_EM,
				   int nd, int ns, int nvs,
				   double *d_E,
				   double* d_FOCK,
				   int cudaBlocks,
				   int cudaThreadsPB,
				   cudaStream_t stream);

extern "C" void deviceCudaEnergyII(int iocc, int jblock,
				   const double *d_C,
				   const double *d_ea,
				   const double *d_ev,
				   const double* d_EM,
				   int nd, int ns, int nvs,
				   double *d_E,
				   double* d_FOCK,
				   int cudaBlocks,
				   int cudaThreadsPB,
				   cudaStream_t stream);

extern "C" void deviceAccum(double *d_E, 
			    int cudaBlocks,
			    int cudaThreadsPB,
			    cudaStream_t stream);


namespace cchem{

namespace rimp2{
namespace detail{


class RIMP2Engine : boost::noncopyable {



    size_t nv_,ns_,nl_,nd_,stride2_,no_;
    const double *ea_,*ev_;
    std::vector<double> &e_m_;
    double *data_pqm1_;

    ::cchem::rimp2::device *pGPU_;

    int ndevice_;
    int nstreams_;
    int cudaBlocks_;
    int cudaThreadsPB_;


    //static load counter for gpu task assignment

    std::vector<cublasHandle_t> dv1_handle_;
    std::vector< cudaStream_t * > dv1_streams_;


    std::vector< std::vector<double*> > dv1_C_;
    std::vector< std::vector<double*> > dv1_E_;
    std::vector< std::vector<double*> > dv1_FOCK_;
    std::vector< std::vector<double*> > dv1_SB_;
    
    std::vector<double*> dv1_ATrans_; //newnewnewnew
    std::vector<double*> dv1_B_;
    std::vector<double*> dv1_ea_;
    std::vector<double*> dv1_ev_;
    std::vector<double*> dv1_EM_;

//vector of device pointers to cholesky decomposed coulomb metric
std::vector<double*> dv1_Lm1_;

    std::vector<int> gpu_counter_;

    std::vector<int> lastBlock_;
    std::vector<int> blocksToMove_;

    std::vector< std::vector<cudaEvent_t> > dv1_events_; 
    int cudaRuntimeVersion_;

    void Allocate(){


 for (int idevice = 0; idevice < ndevice_; idevice++) {

     GPUerrchk ( cudaSetDevice( idevice ) );

  //allocate storage for cudaStreams
  for(int i=0; i<nstreams_; i++){

      //device storage for ri-4center ERIs (each stream creates unique integrals)
      double *ptr_C = 0;
      int matrixSizeC = (nv_+ns_)*(nv_+ns_)*stride2_;
      GPUerrchk( cudaMalloc((void **)&ptr_C, matrixSizeC * sizeof(double)) );
      dv1_C_[idevice][i] = ptr_C;

      //device storage for ri-4center ERIs (each stream creates unique integrals)
      double *ptr_SB = 0;
      int matrixSizeSB = nl_*(nv_+ns_);
      GPUerrchk( cudaMalloc((void **)&ptr_SB, matrixSizeSB * sizeof(double)) );
      dv1_SB_[idevice][i] = ptr_SB;

      //device storage for thread results 
      double *ptr_E = 0;
      int vectorSizeC = cudaBlocks_*cudaThreadsPB_; //8blocks * 64 threads/block (threads/stream)
      GPUerrchk( cudaMalloc((void **)&ptr_E, vectorSizeC * sizeof(double)) );

      //initialize to zero
      GPUerrchk( cudaMemset( ptr_E, 0, vectorSizeC * sizeof(double)) );
      dv1_E_[idevice][i] = ptr_E;

      //fock-like matrix contribution if applicable (for ZAPT)
      double *ptr_FOCK = NULL;
      if(ns_){
  	  int matrixSizeFOCK = nd_*nv_; //no*(nv+ns);
  	  GPUerrchk( cudaMalloc((void **)&ptr_FOCK, matrixSizeFOCK * sizeof(double)) );

  	  //initialize to zero
  	  GPUerrchk( cudaMemset( ptr_FOCK, 0, matrixSizeFOCK * sizeof(double)) );      
      }//ns
      //push null ptr so we do not have any freakouts later
      //      dv_FOCK.push_back(ptr_FOCK);
      dv1_FOCK_[idevice][i] = ptr_FOCK;
      
  }//i

  //shared device storage for B_ia^Q 3-center terms
  double *d_ATrans = 0;
  int matrixSizeATrans = nl_*(nv_+ns_)*stride2_;
  GPUerrchk( cudaMalloc((void **)&d_ATrans, matrixSizeATrans * sizeof(double)) );
  dv1_ATrans_.push_back(d_ATrans);  

  //device storage for cholesky decomposed coulomb metric 
  double *ptr_Lm1 = 0;
  int matrixSizeLm1 = nl_*nl_;
  GPUerrchk( cudaMalloc((void **)&ptr_Lm1, matrixSizeLm1 * sizeof(double)) );
  dv1_Lm1_.push_back(ptr_Lm1);  
  GPUerrchk( cudaMemcpy(dv1_Lm1_[idevice], data_pqm1_, nl_*nl_*sizeof(double),
 			cudaMemcpyHostToDevice) );    


  //storage for active orbital energies (double + single (if applicable) occupied)
  double *d_ea = 0;
  int vectorSizeA = no_;
  GPUerrchk( cudaMalloc((void **)&d_ea, vectorSizeA * sizeof(double)) );
  dv1_ea_.push_back(d_ea);
  GPUerrchk( cudaMemcpy(dv1_ea_[idevice],&ea_[0], no_*sizeof(double),
 			cudaMemcpyHostToDevice) );    
  
   
  //storage for single (if applicable) occupied and virtual orbitals
  double *d_ev = 0;
  int vectorSizeSV = ns_+nv_;
  GPUerrchk( cudaMalloc((void **)&d_ev, vectorSizeSV * sizeof(double)) );
  dv1_ev_.push_back(d_ev);
  GPUerrchk( cudaMemcpy(dv1_ev_[idevice], &ev_[nd_] , (nv_+ns_)*sizeof(double),
 			cudaMemcpyHostToDevice) );    

  
  //modified energy denominators (for ZAPT)
  double *d_EM = 0;
  dv1_EM_.push_back(d_EM);
  if(ns_) {
      int vectorSizeEM = e_m_.size();
      GPUerrchk( cudaMalloc((void **)&d_EM, vectorSizeEM * sizeof(double)) );
      dv1_EM_[idevice] = d_EM;
      GPUerrchk( cudaMemcpy(dv1_EM_[idevice], &(e_m_[0]), e_m_.size()*sizeof(double),
 			    cudaMemcpyHostToDevice) );    
  }//ns


  }//idevice




    }//Allocate


	// static
	// boost::array<size_t,1> memory(const size_t nv, const size_t nd){

	// 	boost::array<size_t,1> memory;
	// 	memory[0] = nd*nv;  //storage for approximate 4center ERIs
	// 	return memory;
	// }


public:

    void Deallocate(double *e_term, Parallel &pe){


	if(ns_){
	Eigen::MatrixXd pfocker(nv_,nd_);
	pfocker.setZero();
	Eigen::MatrixXd pfocker_temp(nv_,nd_);

       for(int idevice = 0; idevice < ndevice_; idevice++) {
       	   for(int istream = 0; istream< nstreams_; istream++) {
       	       GPUerrchk( cudaMemcpy(pfocker_temp.data(),dv1_FOCK_[idevice][istream],
       				     nd_*nv_*sizeof(double),cudaMemcpyDeviceToHost) );
       	       pfocker += pfocker_temp;
       	   }//istream
       }//idevice

       pe.reduce("+",pfocker.data(), (size_t)(nd_*nv_) );

       if(pe.rank() == 0){
       for (size_t d = 0; d < nd_; d++){
	   for(size_t a = 0; a <nv_; a++){

	       e_term[4] +=  pow (pfocker(a,d),2)/(ea_[d]-ev_[nd_+ns_+a]);

	   }//a
       }//d
       e_term[4] = e_term[4]/2.0;
       }//
	}

	for(int idevice = 0; idevice < ndevice_; idevice++){

	    GPUerrchk( cudaSetDevice( idevice ) );
    
	    for(int istream = 0; istream < nstreams_; istream++){
	
		deviceAccum(dv1_E_[idevice][istream],
			    cudaBlocks_,
			    cudaThreadsPB_,
			    dv1_streams_[idevice][istream]);
		double temp_e;
		GPUerrchk( cudaMemcpy(&temp_e, dv1_E_[idevice][istream], 1*sizeof(double), cudaMemcpyDeviceToHost) );
	
		e_term[0] += temp_e;
	
		GPUerrchk( cudaFree(dv1_C_[idevice][istream]) );
		GPUerrchk( cudaFree(dv1_E_[idevice][istream]) );
		GPUerrchk( cudaFree(dv1_SB_[idevice][istream]) );

		if(ns_){
		    GPUerrchk( cudaFree(dv1_FOCK_[idevice][istream]) );
		}//ns
	    }//istream

	    GPUerrchk( cudaFree(dv1_ATrans_[idevice]) );    
	    GPUerrchk( cudaFree(dv1_Lm1_[idevice]) );    
	    GPUerrchk( cudaFree(dv1_ea_[idevice]) );    
	    GPUerrchk( cudaFree(dv1_ev_[idevice]) );    

	    if(ns_){
		GPUerrchk( cudaFree(dv1_EM_[idevice]) );    
	    }//ns

 }//idevice


    }//Deallocate


public:




    RIMP2Engine(int nv, int ns, int nl, int nd, int stride2, int no,
		const double *ea, const double *ev,
		std::vector<double> &e_m,
		double* data_pqm1,
		::cchem::rimp2::device * pGPU):
	nv_(nv),ns_(ns),nl_(nl),nd_(nd),stride2_(stride2),no_(no),
	ea_(ea),ev_(ev),
	e_m_(e_m),
	pGPU_(pGPU),
	data_pqm1_(data_pqm1),
	ndevice_(pGPU_->getNumDevices()),nstreams_(pGPU_->getNumStreams()),
	cudaBlocks_(pGPU_->getCudaBlocks()),cudaThreadsPB_(pGPU_->getCudaThreadsPB()),
	dv1_handle_(pGPU_->cublasHandles()),dv1_streams_(pGPU_->cudaStreams()),
	dv1_C_(std::vector< std::vector<double*> >(ndevice_, std::vector<double*>(nstreams_))),
	dv1_E_(std::vector< std::vector<double*> >(ndevice_, std::vector<double*>(nstreams_))),
 	dv1_FOCK_(std::vector< std::vector<double*> >(ndevice_, std::vector<double*>(nstreams_))),
	dv1_SB_(std::vector< std::vector<double*> >(ndevice_, std::vector<double*>(nstreams_))),
	gpu_counter_(std::vector<int>(ndevice_)),lastBlock_(std::vector<int>(ndevice_)),
	blocksToMove_(std::vector<int>(2*ndevice_)),
	dv1_events_(std::vector< std::vector<cudaEvent_t> >(ndevice_, std::vector<cudaEvent_t>(stride2_))),
	cudaRuntimeVersion_(pGPU->getCudaRuntimeVersion())

	//dv1_events_(std::vector< std::vector<cudaEvent_t> >(ndevice_, std::vector<cudaEvent_t>(nstreams_))),
    {

	//cudaRuntimeGetVersion(&cudaRuntimeVersion_);

	// for(int idevice = 0; idevice < ndevice_; idevice++){

	//     GPUerrchk( cudaSetDevice( idevice ) );
	//     for(int istream = 0; istream < nstreams_; istream++){
	// 	cudaEventCreateWithFlags( &dv1_events_[idevice][istream], cudaEventDisableTiming );
	// 	//cudaEventCreate( &dv1_events_[idevice][istream] );
	//     }//istream
	// }//idevice


     // for(int idevice = 0; idevice < ndevice_; idevice++){
	 
     // 	 GPUerrchk( cudaSetDevice( idevice ) );
     // 	 for(int istream = 0; istream < stride2_; istream++){
     // 	     cudaEventCreateWithFlags( &dv1_events_[idevice][istream], cudaEventDisableTiming );
     // 	 }//istream
     // }//idevice

	this->Allocate();
    };

    ~RIMP2Engine(){

	// for(int idevice = 0; idevice < ndevice_; idevice++){
	//     GPUerrchk( cudaSetDevice( idevice ) );
	//     for(int istream = 0; istream < nstreams_; istream++){
	// 	cudaEventDestroy( dv1_events_[idevice][istream] );
	//     }//istream
	// }//idevice

     // for(int idevice = 0; idevice < ndevice_; idevice++){
     // 	 GPUerrchk( cudaSetDevice( idevice ) );
     // 	 GPUerrchk( cudaDeviceSynchronize() );

     // 		 // for(int i = 0; i < nstreams_; i++){
     // 		 // std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
     // 		 // 	   << " " ;}
     // 		 // std::cout << std::endl;

     // 	 //for(int istream = 0; istream < nstreams_; istream++){
     // 	 for(int istream = 0; istream < stride2_; istream++){
     // 	     cudaEventDestroy( dv1_events_[idevice][istream] );
     // 	 }//istream
     // }//idevice


    // //make sure GPU is done forming C_ia^Q (stored in dv1_ATrans_)
    //  for(int idevice = 0; idevice < ndevice_; idevice++){
    //  	 GPUerrchk( cudaSetDevice( idevice ) );

    // 		 for(int i = 0; i < nstreams_; i++){
    // 		     std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
    // 		 	       << " " ;}
    // 		 std::cout << std::endl;

    //  	 GPUerrchk( cudaDeviceSynchronize() );

    // 		 for(int i = 0; i < nstreams_; i++){
    // 		     std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
    // 		 	       << " " ;}
    // 		 std::cout << std::endl;

    // 	 for(int istream = 0; istream < stride2_; istream++){
    // 	     GPUerrchk( cudaEventDestroy( dv1_events_[idevice][istream] ) );
    // 	 }//istream

    // 	  for(int i = 0; i < nstreams_; i++){
    // 	           std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
    // 			     << " " ;}
    // 	   std::cout << std::endl;
    // 	   std::cout << "------------" << std::endl;

    //  }//idevice

}


    void operator()(int jblock_end,int jblock, //std::vector<int> &lastBlock,
		    double *ptr_bj, size_t nvs,int j_remain)

    {

	//     std::vector<int> blocksToMove(2*ndevice_);

#pragma omp single
     {
     for (int idevice = 0; idevice < ndevice_; idevice++){
	 gpu_counter_[idevice] = 0;
	 lastBlock_[idevice] = 0;
     }


     // if( ndevice_ > (jblock_end-jblock) ){
     // 	 for(int idevice = 0; idevice < ndevice_; idevice ++){
     // 	     // GPUerrchk( cudaMemcpy(dv1_A[idevice],
     // 	     // 			   ptr_bj,
     // 	     // 			   j_remain*nl*nvs*sizeof(double),
     // 	     // 			   cudaMemcpyHostToDevice) );
     // 	     //	     	     std::cout << "you better worry! " << __FILE__ << __LINE__ << std::endl;
     // 	     //this is probably not an issue and if this happens, you are wasting devices.
     // 	     //if this error comes up you will need to set this up so that the device 
     // 	     //  which to not do work in the next section, that you properly populate
     // 	     //    its dv1_ATrans vector (here)
     // 	 }
     // }


     }//omp single

     // //create cuda events
     // for(int idevice = 0; idevice < ndevice_; idevice++){
     // 	 GPUerrchk( cudaSetDevice( idevice ) );
     // 	 for(int istream = 0; istream < stride2_; istream++){
     // 	     GPUerrchk( 
     // 		       cudaEventCreateWithFlags( 
     // 						&dv1_events_[idevice][istream], 
     // 						cudaEventDisableTiming ) 
     // 			);
     // 	 }//istream
     // }//idevice
     
     //std::cout << jblock << " " << jblock_end << std::endl;
     //#pragma omp for schedule(dynamic,1)  
     for(int iocc = jblock; iocc < jblock_end; iocc++){  

	 double alpha = 1.0;
	 double beta = 0.0;

	 int hard_device = (iocc-jblock)%ndevice_;
	 GPUerrchk( cudaSetDevice( hard_device ) );

	 int new_stream = gpu_counter_[hard_device]%nstreams_;

	 gpu_counter_[ hard_device ] += 1;

	 cublasSetStream(dv1_handle_[hard_device], dv1_streams_[hard_device][new_stream]);

	 int blocks = ndevice_;
	 if( (iocc + ndevice_) >= jblock_end) blocks = jblock_end - lastBlock_[hard_device]-jblock;




	 //this set up an array of I(nvs*nl) blocks to retrieve
	 // 1) these are stored per stream
	 // 2) transformed as (ia|Q) -> C_ia_Q to the dv1_ATrans vector.
	 // 3) the (ia|Q) block that is used here is the last one to be retieved 
	 for(int i = lastBlock_[hard_device], ipos=0; 
	     i < lastBlock_[hard_device] + blocks; 
	     i++,ipos++){
	     blocksToMove_[ipos]=i;
	 }//i
	 std::swap(blocksToMove_[hard_device],blocksToMove_[blocks-1]);

	 
	 for(int ipos = 0;
	     ipos < blocks;
	     ipos++){
	     int iblock = blocksToMove_[ipos];

	     //time to retrieve
	     GPUerrchk( cudaMemcpy2DAsync(dv1_SB_[hard_device][new_stream],
	 				  nvs*sizeof(double),
	 				  &ptr_bj[ (iblock)*nvs],
	 				  j_remain*nvs*sizeof(double),
	 				  nvs*sizeof(double),nl_,
	 				  cudaMemcpyHostToDevice,
	 				  dv1_streams_[hard_device][new_stream]) );

	     //transform now
	     cublasDsymm(dv1_handle_[hard_device], CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, 
	     		 (int)nvs, (int)nl_, 
	     		 &alpha, dv1_Lm1_[hard_device], nl_,
	     		 dv1_SB_[hard_device][new_stream], nvs, 
	     		 &beta, &dv1_ATrans_[hard_device][ (iblock)*nvs], j_remain*nvs);

	     // cublasDgemm(dv1_handle_[hard_device], CUBLAS_OP_N, CUBLAS_OP_N, 
	     // 		 (int)nvs, (int)nl_, nl_,
	     // 		 &alpha, dv1_SB_[hard_device][new_stream], nvs, 
	     // 		 dv1_Lm1_[hard_device], nl_,
	     // 		 &beta, &dv1_ATrans_[hard_device][ (iblock)*nvs], j_remain*nvs);
	     
	 }//iblock

	 //update how many blocks are on the device
	 lastBlock_[hard_device] += blocks;

	 //you have to synchronize streams if you have more than one	 
	 if(nstreams_ > 1)
	     GPUerrchk( cudaStreamSynchronize( dv1_streams_[hard_device][new_stream]) );



	 //if you ever have a problem with GPU ri-mp2 energies, suspect the next big of code.
	 //  it appear that cuda events did not work well with cuda 5.0, cuda 6.5 is OK
	 //  the device sychronizations below are your best bet to fix any issues
	 // if( (cudaRuntimeVersion_ < 6050) 
	 //     && (nstreams_ > 1) ){

	 //     GPUerrchk( cudaSetDevice( hard_device ) );
	 //     GPUerrchk( cudaDeviceSynchronize() );

	 // }else{
	 //     if(nstreams_ > 1){


	 // 	 // for(int i = 0; i < nstreams_; i++){
	 // 	 //     std::cout << cudaStreamQuery( dv1_streams_[hard_device][i]) 
	 // 	 // 	       << " " ;}
	 // 	 // std::cout << std::endl;
		 
	 // 	 // GPUerrchk( cudaStreamSynchronize( dv1_streams_[hard_device][new_stream]) );
		 
	 // 	 // for(int i = 0; i < nstreams_; i++){
	 // 	 //     std::cout << cudaStreamQuery( dv1_streams_[hard_device][i]) 
	 // 	 // 	       << " " ;}
	 // 	 // std::cout << std::endl;
		 
	 // 	 // std::cout << "------------------" << std::endl;


	 // 	 // for(int i = 0; i < nstreams_; i++){
	 // 	 //     std::cout << cudaStreamQuery( dv1_streams_[hard_device][i]) 
	 // 	 // 	       << " " ;}
	 // 	 // std::cout << std::endl;

	 // 	 // cudaEventRecord( dv1_events_[hard_device][ (new_stream+1)%nstreams_ ], 
	 //  	 //  		  dv1_streams_[hard_device][new_stream]);
	 // 	 // cudaStreamWaitEvent( dv1_streams_[hard_device][new_stream],
	 //  	 //  		      dv1_events_[hard_device][new_stream],
	 //  	 //  		      0) ;

	 // 	 //iocc-jblock
	 // 	 //gpu_counter_[ hard_device ]
	 // 	 GPUerrchk( 
	 // 		   cudaEventRecord( dv1_events_[hard_device][iocc-jblock ], 
	 // 				    dv1_streams_[hard_device][new_stream])
	 // 		    );
			   
	 // 	 if(iocc-jblock -ndevice_ +1 > 0)
	 // 	     GPUerrchk( 
	 // 		       cudaStreamWaitEvent( dv1_streams_[hard_device][new_stream],
	 // 					    dv1_events_[hard_device][iocc-jblock -ndevice_ ],
	 // 					    0)
	 // 			);


	 // 	 // for(int i = 0; i < nstreams_; i++){
	 // 	 //     std::cout << cudaStreamQuery( dv1_streams_[hard_device][i]) 
	 // 	 // 	       << " " ;}
	 // 	 // std::cout << std::endl;
		 
	 // 	 // std::cout << "------------------" << std::endl;


	 //  	 // cudaStreamWaitEvent( dv1_streams_[hard_device][new_stream],
	 //  	 // 		      dv1_events_[hard_device][new_stream],
	 //  	 // 		      0) ;
	 //  	 // cudaEventRecord( dv1_events_[hard_device][ (new_stream+1)%nstreams_ ], 
	 //  	 // 		  dv1_streams_[hard_device][new_stream]);
	 //     }//nstreams_
	 // }//(cudaRuntimeVersion_ < 6050)

	 // std::cout << iocc << " " << jblock_end<<" "<< iocc - jblock << std::endl;
	 // std::cout<<hard_device<<" "<< new_stream <<" "<<(new_stream+1)%nstreams_<<std::endl;

 	 //form ERIs 
	 cublasDgemm(dv1_handle_[hard_device], CUBLAS_OP_N, CUBLAS_OP_T, 
	 	     nvs, (iocc-jblock+1)*nvs, nl_, 
	 	     &alpha, dv1_SB_[hard_device][new_stream], nvs, 
	 	     dv1_ATrans_[hard_device],j_remain*nvs, 
	 	     &beta, dv1_C_[hard_device][new_stream],nvs);

	 //accumulate energy contributions on device
 	 deviceCudaEnergyII(iocc, jblock, (const double*)dv1_C_[hard_device][new_stream],
	 		    (const double*)dv1_ea_[hard_device],
	 		    (const double*)dv1_ev_[hard_device],
	 		    (const double*)dv1_EM_[hard_device],
	 		    (int)nd_, (int)ns_, (int)nvs,
	 		    dv1_E_[hard_device][new_stream],
	 		    dv1_FOCK_[hard_device][new_stream],
	 		    cudaBlocks_,
	 		    cudaThreadsPB_,
	 		    dv1_streams_[hard_device][new_stream] );
	 
     }//iocc


}//operator


    void operator()(int num_chunks, int i_remain, int j_remain,int stride_row,
		    double *ptr_bi, double *ptr_bj,int nvs, int iblock,int jblock)
{

    // //make sure GPU is done forming C_ia^Q (stored in dv1_ATrans_)
    //  for(int idevice = 0; idevice < ndevice_; idevice++){


    //  	 GPUerrchk( cudaSetDevice( idevice ) );

    // 		 for(int i = 0; i < nstreams_; i++){
    // 		     std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
    // 		 	       << " " ;}
    // 		 std::cout << std::endl;


    //  	 GPUerrchk( cudaDeviceSynchronize() );

    // 		 for(int i = 0; i < nstreams_; i++){
    // 		     std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
    // 		 	       << " " ;}
    // 		 std::cout << std::endl;

    // 	 for(int istream = 0; istream < stride2_; istream++){
    // 	     GPUerrchk( cudaEventDestroy( dv1_events_[idevice][istream] ) );
    // 	 }//istream

    // 	  for(int i = 0; i < nstreams_; i++){
    // 	           std::cout << cudaStreamQuery( dv1_streams_[idevice][i]) 
    // 			     << " " ;}
    // 	   std::cout << std::endl;
    // 	   std::cout << "------------" << std::endl;

    //  }//idevice


#pragma omp single
       {
	   for (int idevice = 0; idevice < ndevice_; idevice++){
	       gpu_counter_[idevice] = 0;
	   }//idevice

       }//omp single


#pragma omp for schedule(dynamic,1)
       for (int ichunk = 0; ichunk < num_chunks; ichunk++){
	   int chunk_len = stride_row;
	   int chunk_pos = ichunk*stride_row;
	   if( (ichunk + 1) == num_chunks) chunk_len = i_remain - (num_chunks-1)*stride_row;

	   double alpha = 1.0;
	   double beta = 0.0;

	   GPUerrchk( cudaSetDevice( ichunk%ndevice_ ) );

	   int new_stream = gpu_counter_[ichunk%ndevice_]%nstreams_;
	   gpu_counter_[ichunk%ndevice_] += 1;

	   GPUerrchk( cudaMemcpy2DAsync(dv1_SB_[ichunk%ndevice_][new_stream] ,
					nvs*sizeof(double),
					&ptr_bi[ichunk*nvs],
					i_remain*nvs*sizeof(double),
					nvs*sizeof(double),nl_, //width, height (rows)
					cudaMemcpyHostToDevice,
					dv1_streams_[ichunk%ndevice_][new_stream]) );

	   cublasSetStream(dv1_handle_[ichunk%ndevice_], dv1_streams_[ichunk%ndevice_][new_stream]);

	   cublasDgemm(dv1_handle_[ichunk%ndevice_], CUBLAS_OP_N, CUBLAS_OP_T, 
	   	       chunk_len*nvs, stride2_*nvs, nl_, 
	   	       &alpha, dv1_SB_[ichunk%ndevice_][ new_stream ], nvs, 
	   	       dv1_ATrans_[ichunk%ndevice_], j_remain*nvs, 
	   	       &beta, dv1_C_[ichunk%ndevice_][new_stream], nvs );

	   deviceCudaEnergyIJ( iblock+chunk_pos, jblock, jblock+stride2_, 
	   		       (const double*)dv1_C_[ichunk%ndevice_][new_stream],
	   		       (const double*)dv1_ea_[ichunk%ndevice_],
	   		       (const double*)dv1_ev_[ichunk%ndevice_],
	   		       (const double*)dv1_EM_[ichunk%ndevice_],
	   		       (int)nd_, (int)ns_, (int)nvs,
	   		       dv1_E_[ichunk%ndevice_][new_stream],
	   		       dv1_FOCK_[ichunk%ndevice_][new_stream],
	   		       cudaBlocks_,
	   		       cudaThreadsPB_,
	   		       dv1_streams_[ichunk%ndevice_][new_stream] );

       }//ichunk

       //       fineTimer.reset();
       for (int idevice = 0; idevice < ndevice_; idevice++){
	   GPUerrchk( cudaSetDevice( idevice ) );
	   GPUerrchk( cudaDeviceSynchronize() );
       }//idevice
       //       profile.gpuSyncEri += fineTimer;


    }//operator()




};//class RIMP2Engine



}//namespace detail
}//namespace rimp2
}//namespace cchem



#endif//HAVE_CUBLAS
















#endif /* SRC_RI_MP2_ERICONTRACTIONENERGYACCUM_HPP_ */
