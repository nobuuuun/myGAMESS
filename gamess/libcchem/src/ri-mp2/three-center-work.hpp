/*
 * thee-center-work.hpp
 *
 *  Created on: Nov 5, 2015
 *      Author: luke
 */

#ifndef SRC_RI_MP2_THEE_CENTER_WORK_HPP_
#define SRC_RI_MP2_THEE_CENTER_WORK_HPP_


#include "god.hpp"
#include "thread.hpp"
#include <Eigen/Dense>
#include <basis/basis.hpp>
#include "ri-mp2/ri-integrals.hpp"

#include "ri-async.hpp"
#include <math.hpp>

namespace cchem{
namespace rimp2_gradient{

namespace detail{

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;
typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef const Eigen::Map<const Eigen::MatrixXd,Eigen::AutoAlign> ConstMapMatrixXd;
typedef ::basis::Basis Basis;
using cchem::Thread;




/*
  @brief three-center ERIs
  @author LBR
  @detail interface to generated 3-center ERIs, perform AO->MO transformation,
  and build the mixed lagrangian (AO,virt) and (occupied,AO)
  -CPU only
  @param basis OBS
  @param auxbasis ABS
  @param Ca occupied LCAO coefficients
  @param Cv virtual LCAO coefficients
  @param pe parallel environment
  @param spherical flag to indicate whether ABS uses spherical guassians
  @param max_auxbasis_shell_size number of functions in highest moment shell
  @param no number of occupied MOs
  @param ns number of singly occupied MOs
  @param nv number of virtual MOs
  @param N number of OBS functions
  @param BARE_IA array pointer to (I A|L) integrals
  @param GAMMA_INU_P array pointer to mixed gamma
  @param GAMMA_IA_Q array pointer to MO gamma
  @param BARE_MN_SYM array pointer to (mu nu|L) integrals
*/
class ThreeCenterWork : boost::noncopyable {

	boost::reference_wrapper<const Basis> &basis_;
	boost::reference_wrapper<const Basis> &auxbasis_;
	boost::reference_wrapper<const Matrix> &Ca_;
	boost::reference_wrapper<const Matrix> &Cv_;
	Parallel &pe_;
	int spherical_;
	int max_auxbasis_shell_size_;
	size_t no_;
	size_t ns_;
	size_t nv_;
	size_t N_;
	Array<double> *BARE_IA_;
	Array<double> *GAMMA_INU_P_;
	Array<double> *GAMMA_IA_Q_;
//	Array<double> *BARE_MN_;
	Array<double> *BARE_MN_SYM_;

public:
	ThreeCenterWork(boost::reference_wrapper<const Basis> basis, 
			boost::reference_wrapper<const Basis> auxbasis,
			boost::reference_wrapper<const Matrix> Ca, 
			boost::reference_wrapper<const Matrix> Cv,
			Parallel &pe, int spherical,int max_auxbasis_shell_size, 
			size_t no, size_t ns, size_t nv, size_t N,
			Array<double> *BARE_IA,
			Array<double> *GAMMA_INU_P,  Array<double> *GAMMA_IA_Q, 
			Array<double> *BARE_MN_SYM):
	    basis_(basis),auxbasis_(auxbasis),Ca_(Ca),Cv_(Cv),
	    pe_(pe), spherical_(spherical),max_auxbasis_shell_size_(max_auxbasis_shell_size),
	    no_(no),ns_(ns),nv_(nv), N_(N), BARE_IA_(BARE_IA),
	    GAMMA_INU_P_(GAMMA_INU_P),GAMMA_IA_Q_(GAMMA_IA_Q), BARE_MN_SYM_(BARE_MN_SYM){}

    void GetAOIntBatch(std::vector<MapMatrixXd> &sph_c,
		       std::vector<MapMatrixXd> &threec_eri_mat_vec,
		       cchem::ri::AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri > &auxiliary_eri,
		       const Basis::Shell &L);
    
    void BuildMOInts(std::vector<MapMatrixXd> &sph_c);
    
    void BuildMixedLagrangian(double *ptr_gamma_ia_P, double *ptr_gamma_inu_P,
			      std::vector<MapMatrixXd> &sph_c,
			      MapMatrixXd &lag_mu_i,
			      MapMatrixXd &lag_a_nu);

}; //class ThreeCenterWork




}//namespace detail
}//namespace rimp2_gradient







namespace threeCenterInt {

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;
typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;

typedef ::basis::Basis Basis;
using cchem::Thread;



/*
  @brief CPU three-center ERI
  @author LBR
  @detail CPU computes ERI performs the transformation
  @param basis OBS
  @param auxbasis ABS
  @param spherical flag to indicate whether ABS uses spherical functions
  @param threeCartToSphere pointer to cartesian to spherical transformation matrix
  @param N number of OBS functions
*/

class threeCenterEriWork : boost::noncopyable {


	const Basis &basis_;
	const Basis &auxbasis_;
	const int spherical_;
	double *threeCartToSph_;
	const size_t N_;

public:
	threeCenterEriWork(boost::reference_wrapper<const Basis> basis,
			boost::reference_wrapper<const Basis> auxbasis,
			const int &spherical, double *threeCartToSph, const size_t &N):
				basis_(basis.get()),auxbasis_(auxbasis.get()),spherical_(spherical),
				threeCartToSph_(threeCartToSph),N_(N){
	}

	void getThreeCenterEriBatch(const size_t &domain,
			const std::vector< size_t > &blockRanges,
			cchem::ri::AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri > &auxiliary_eri,
			std::vector<MapMatrixXd> &sph_c,
			double *ptr_buff1, const size_t &offsetAO);

};//class threeCenterEriWork





typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;

#if !HAVE_CUBLAS

/*
  @brief CPU three-center AO->MO transformation
  @author LBR
  @detail CPU three-center AO->MO transformation
  @param Ca occupied LCAO coefficients
  @param Cv virtual LCAO coefficients
  @param auxbasis ABS
  @param spherical flag that indicates whether the ABS uses spherical guassians
  @param blockRankes vector of blocks of ABS shells
  @param nvs number of singly and virtual MOs
  @param no number of occupied orbitals
  @param N number of OBS functions
  @param async reference to object used for asynchronous operations
  @param Vp array pointer to transformed ERIs
*/

class threeCenterTransform : boost::noncopyable {

    Thread thread_;
    const Matrix &Ca_,&Cv_;
    const Basis &auxbasis_;
    const int &spherical_;
    const std::vector< size_t > &blockRanges_;
    const size_t &nvs_, &no_, &N_;
    cchem::ri::async &async_;
    Array<double> *Vp_;

public:
	static
	boost::array<size_t,1> memory(const size_t &nv,const size_t &N){

		boost::array<size_t,1> memory;
		memory[0] = N*nv;  //scratch for 1/3 transformation
		return memory;
	}

	threeCenterTransform(boost::reference_wrapper<const Matrix> Ca,
			boost::reference_wrapper<const Matrix> Cv,
			boost::reference_wrapper<const Basis> auxbasis,
			const int &spherical, const std::vector<size_t> &blockRanges,
			const size_t &nvs, const size_t &no, const size_t &N,
		    cchem::ri::async &async, Array<double> *Vp):
				Ca_(Ca.get()),Cv_(Cv.get()),auxbasis_(auxbasis.get()),
				spherical_(spherical),blockRanges_(blockRanges),
				nvs_(nvs),no_(no),N_(N),async_(async), Vp_(Vp){

	    // for(int i = 0; i < 10; i++)
	    // std::cout << "using cpu transform" << std::endl;
		BOOST_AUTO(memory, this->memory(nvs_,N_));
		for (int i = 0; i < memory.size(); ++i) {
			thread_.data[i] = thread_.malloc(memory[i]);
		}//memory

	}//threeCenterTransform

	~threeCenterTransform(){}


	void transformBlock(const size_t &domain,
			const size_t &offsetAO, const size_t &offsetMO,
			double *ptr_buff2 ,double *ptr_buff4)
	{
#pragma omp for
		for(size_t task = blockRanges_[domain]; task < blockRanges_[domain+1]; task++ ){
			const Basis::Shell &L = auxbasis_.shells().at(task);

			int lSize = L.size();
			int lStart = L.start();
			int lStop = L.stop();
			if(spherical_){
				lSize = L.sphsize();
				lStart = L.sphstart();
				lStop = L.sphstop();
			}

			double *ptr_scr = thread_.data[0];
			for(int l=0; l<lSize; l++){

				cblas_dsymm(CblasColMajor,CblasRight,CblasLower,
						(int)no_,(int)N_,
						1.0, &ptr_buff2[ l*N_*N_ + lStart*N_*N_ -offsetAO], (int)N_,
						Ca_.data().begin(), (int)no_,
						0.0,ptr_scr,(int)no_ );

				cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
						(int)(nvs_),(int)no_,(int)N_,
						1.0,Cv_.data().begin(), (int)(nvs_),
						ptr_scr,(int)no_,
						0.0,&ptr_buff4[ l*no_*nvs_ + lStart*no_*nvs_ - offsetMO ],(int)(nvs_) );
			}//l

		}//task :auxiliary shells


	};//transformBlock


	void writeBlock(const size_t &domain, double *ptr_buff3){

	    const Basis::Shell &lFirst = auxbasis_.shells().at(blockRanges_[domain]);
	    const Basis::Shell &lLast = auxbasis_.shells().at(blockRanges_[domain+1]-1);
	    size_t lPutStart = lFirst.start();
	    size_t lPutStop = lLast.stop();
	    if(spherical_){
	   	lPutStart = lFirst.sphstart();
	   	lPutStop = lLast.sphstop();
	    }//(spherical)



	    std::vector<size_t> startPut(2), finishPut(2);
	    startPut[0] = 0;
	    startPut[1] = lPutStart;
	    finishPut[0] = no_*nvs_;
	    finishPut[1] = lPutStop;

	    // std::cout << "about to send block off for writing " 
	    // 	      << startPut[0] << " "
	    // 	      << startPut[1] << " "
	    // 	      << finishPut[0] << " "
	    // 	      << finishPut[1] << " "
	    // 	      << std::endl << std::cout;

	    async_.put( startPut, finishPut, *Vp_, ptr_buff3);




	}


};//class threeCenterTransform

#endif //!HAVE_CUBLAS





#if HAVE_CUBLAS

    inline void GPUassert(cudaError_t code, char const * file, int line, bool Abort=true)
    {
	if (code != 0) {
	    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),file,line);
	    if (Abort) exit(code);
	}       
    }
    
#define GPUerrchk(ans) { GPUassert((ans), __FILE__, __LINE__); }



/*
  @brief GPU three-center AO->MO transformation
  @author LBR
  @detail GPU three-center AO->MO transformation
  @param Ca occupied LCAO coefficients
  @param Cv virtual LCAO coefficients
  @param auxbasis ABS
  @param spherical flag that indicates whether the ABS uses spherical guassians
  @param blockRankes vector of blocks of ABS shells
  @param nvs number of singly and virtual MOs
  @param no number of occupied orbitals
  @param N number of OBS functions
  @param pGPU pointer to device constructs
  @param async reference to object used for asynchronous operations
  @param Vp array pointer to transformed ERIs
*/

class threeCenterTransform : boost::noncopyable {

    Thread thread_;
    const Matrix &Ca_,&Cv_;
    const Basis &auxbasis_;
    const int &spherical_;
    const std::vector< size_t > &blockRanges_;
    const size_t &nvs_, &no_, &N_;
    const int ndevice_,nstreams_;
    
    //vector of device pointers to (mu nu|Q) integrals (per stream)
    std::vector<double*> dvMuNuQ;
    
    //vector of device pointers to 1/3 transformed integrals (per stream)
    std::vector<double*> dvINuQ;

    //vector of device pointers to (I A|Q) integrals (per stream)
    std::vector<double*> dvIAQ;

    //vector of device pointers to occupied LCAO coefficients(per device)
    std::vector<double*> dvCoeffOcc;

    //vector of device pointers to singly occupied & virtual LCAO coefficients(per device)
    std::vector<double*> dvCoeffSV;

    //static load counter for gpu task assignment
    std::vector<int> gpu_counter;

    std::vector<cublasHandle_t> dv1_handle_;

    std::vector< cudaStream_t * > dv1_streams_;

    cchem::ri::async &async_;

    Array<double> *Vp_;


    void GPUAllocate();

    void GPUDeallocate();

public:

	threeCenterTransform(boost::reference_wrapper<const Matrix> Ca,
			     boost::reference_wrapper<const Matrix> Cv,
			     boost::reference_wrapper<const Basis> auxbasis,
			     const int &spherical, const std::vector<size_t> &blockRanges,
			     const size_t &nvs, const size_t &no, const size_t &N,
			     const int ndevice, const int nstreams,
			     const std::vector<cublasHandle_t> &dv1_handle,
			     std::vector< cudaStream_t * > &dv1_streams,
			     cchem::ri::async &async, Array<double> *Vp):
				Ca_(Ca.get()),Cv_(Cv.get()),auxbasis_(auxbasis.get()),
				spherical_(spherical),blockRanges_(blockRanges),
				nvs_(nvs),no_(no),N_(N),ndevice_(ndevice), nstreams_(nstreams),
				dv1_handle_(dv1_handle), dv1_streams_(dv1_streams),
				async_(async), Vp_(Vp)
    {

	    // for(int i = 0; i < 10; i++)
	    // 	std::cout << "using GPU transform " << std::endl;

	GPUAllocate();

	}//threeCenterTransform

    ~threeCenterTransform(){GPUDeallocate();}


	void transformBlock(const size_t &domain,
			const size_t &offsetAO, const size_t &offsetMO,
			double *ptr_buff2 ,double *ptr_buff4)
	{

			    for(int idevice = 0; idevice < ndevice_; idevice++){
				gpu_counter[idevice] = 0;
			    }//idevice

			    for(size_t task = blockRanges_[domain]; task < blockRanges_[domain+1]; task++ ){
				const Basis::Shell &L = auxbasis_.shells().at(task);



				int lSize = L.size();
				int lStart = L.start();
				int lStop = L.stop();
				if(spherical_){
				    lSize = L.sphsize();
				    lStart = L.sphstart();
				    lStop = L.sphstop();
				}

				for(int l=0; l<lSize; l++){

				    int hard_device = (task+l)%ndevice_;
				    GPUerrchk( cudaSetDevice( hard_device ) );

				    int new_stream = gpu_counter[hard_device]%nstreams_;
				    gpu_counter[ hard_device ] += 1;
				    cublasSetStream(dv1_handle_[hard_device], dv1_streams_[hard_device][new_stream]);

				    double alpha = 1.0;
				    double beta = 0.0;

				    GPUerrchk( cudaMemcpyAsync(&dvMuNuQ[hard_device][new_stream*N_*N_],
							       &ptr_buff2[ l*N_*N_ + lStart*N_*N_ -offsetAO], 
							       N_*N_*sizeof(double),
							       cudaMemcpyHostToDevice,
							       dv1_streams_[hard_device][new_stream]) );;    
				    
   cublasDsymm(dv1_handle_[hard_device], CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, 
	       (int)no_, (int)N_, 
	       &alpha, &dvMuNuQ[hard_device][new_stream*N_*N_], (int)N_, 
	       dvCoeffOcc[hard_device], (int)no_, 
	       &beta, &dvINuQ[hard_device][new_stream*no_*N_], (int)no_);
   
   cublasDgemm(dv1_handle_[hard_device], CUBLAS_OP_N, CUBLAS_OP_T, 
	       (int)nvs_, (int)no_, (int)N_, 
	       &alpha, dvCoeffSV[hard_device], (int)nvs_,
	       &dvINuQ[hard_device][new_stream*no_*N_], (int)no_,
	       &beta, &dvIAQ[hard_device][new_stream*no_*nvs_], nvs_);

				    GPUerrchk( cudaMemcpyAsync(&ptr_buff4[ l*no_*nvs_ + lStart*no_*nvs_ - offsetMO ], 
							       &dvIAQ[hard_device][new_stream*no_*nvs_],       		  
							       no_*nvs_*sizeof(double),
							       cudaMemcpyDeviceToHost,
							       dv1_streams_[hard_device][new_stream]) );    

				}//l

			    }//task :auxiliary shells


	};//transformBlock


	void writeBlock(const size_t &domain, double *ptr_buff3){


	    const Basis::Shell &lFirst = auxbasis_.shells().at(blockRanges_[domain]);
	    const Basis::Shell &lLast = auxbasis_.shells().at(blockRanges_[domain+1]-1);
	    size_t lPutStart = lFirst.start();
	    size_t lPutStop = lLast.stop();
	    if(spherical_){
	   	lPutStart = lFirst.sphstart();
	   	lPutStop = lLast.sphstop();
	    }//(spherical)

	    std::vector<size_t> startPut(2), finishPut(2);
	    startPut[0] = 0;
	    startPut[1] = lPutStart;
	    finishPut[0] = no_*nvs_;
	    finishPut[1] = lPutStop;

	    // std::cout << "about to send block off for writing " 
	    // 	      << startPut[0] << " "
	    // 	      << startPut[1] << " "
	    // 	      << finishPut[0] << " "
	    // 	      << finishPut[1] << " "
	    // 	      << std::endl << std::cout;

	    async_.put( startPut, finishPut, *Vp_, ptr_buff3);




	}

};//class threeCenterTransform

#endif //HAVE_CUBLAS










}//namespace ThreeCenterInt


}//namespace cchem



#endif /* SRC_RI_MP2_THEE_CENTER_WORK_HPP_ */
