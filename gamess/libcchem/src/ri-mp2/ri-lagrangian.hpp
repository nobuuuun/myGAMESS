/*
 * ri-lagrangian.hpp
 *
 *  Created on: Nov 11, 2015
 *      Author: luke
*/

#ifndef SRC_RI_MP2_RI_LAGRANGIAN_HPP_
#define SRC_RI_MP2_RI_LAGRANGIAN_HPP_


#include "god.hpp"

#include <JK.hpp>
#include <boost/ref.hpp>

#include "thread.hpp"

#include <Eigen/Dense>


#include <ri-async.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/ref.hpp>

#include <device.hpp>
#if HAVE_CUBLAS
#include <deviceGradientEngine.hpp>
#endif

//#include <math.hpp>
namespace cchem {
namespace rimp2_gradient {
namespace detail {

/*
  @brief CPU: build relaxed density, build Gamma, compute the correlation energy
  @author LBR
  @detail 
  -build active-active/virtual/virtual block of relaxed density matrix
  -build gamma
  -compute the correlation energy
  @param nf number fo frozen orbitals
  @param no number of occupied orbitals
  @param nv number of virtual orbitals
  @param nl number of auxiliary functions
  @param ea active orbital energies
  @param ev virtual orbital energies
  @param energy correlation energy
  @param pmat relaxed density matrix
  @param pGammaLong host pointer to active orbital i component of Gamma  
  @param gamma_rs host pointer to Gamma
  @param data_pqm1 pointer to metric
*/
class CPUPmatGammaFunctor : boost::noncopyable {
private:

    //dimensions
    const size_t nf_,no_,nv_,nl_,na_;
    //energy
    double &energy_;
    //orbital energies
    const double *ea_,*ev_;
    //[na][nv*nv]
    //std::vector<double *> &pTijLong_;
    // [ (no+nv), (no+nv) ]
    MapMatrixXd &pmat_;
    //[na][nv*nl]
    std::vector<double *> &pGammaLong_;
    //gamma_rs array
    // [nl,nl]
    double *gamma_rs_;
    //coulomb metric
    double *data_pqm1_;


    //number of CPU threads allowed
    int max_threads_;

    //the following pointers are declared here to make deletion (deallocate) easier
    //temporary storage of 3/3 transformed integrals
    double *pTrans_;
    //pointer eri storage
    double *pEri_;
    //pointer to virtual block of pmat storage
    double *pPmatVirt_;
    //pointer to active block of pmat storage
    double *pPmatOcc_;
    //pointer to gamma storage
    double *pGamma_;
    //scratch space for pPmatVirt contribution
    //double *pAC_;
    //scratch space for pPmatVirt contribution
    //double *pBC_;
    //storage for modified ERIs
    double *pEriLong_;
    //storage for tijab amplitudes
    double *pTijLong_;         
    //shared storage for modified ERIs
    std::vector<double *> pVecEriLong_;
    //shared storage for modified ERIs
    std::vector<double *> pVecTijLong_;


    //these arrays are for thread work
    //thread pointer  storage for constructed ERIs
    std::vector<double *> pVecEri_;
    //thread pointer  storage for virtual block of pmat
    std::vector<double *> pVecPmatVirt_;
    //thread pointer storage for active block of pmat
    std::vector<double *> pVecPmatOcc_;
    //thread pointer storage for Gamma_ia^Q
    std::vector<double *> pVecGamma_;
    //thread storage for energy contribution
    std::vector<double> threadEnergy_;

public:
    
    static size_t CPUMem(const size_t nv, const size_t na, 
			 const size_t nl, const size_t no){

	int maxThreads = omp_get_max_threads(); 
	size_t reqMem = 0;

	// pEri_
	reqMem += maxThreads*nv*nv*sizeof(double);
	// pEriLong_
	reqMem += na*nv*nv*sizeof(double);
	// pTijLong_
	reqMem += na*nv*nv*sizeof(double);
	// pPmatVirt_
	reqMem += maxThreads*nv*nv*sizeof(double);
	// pPmatOcc_
	reqMem += maxThreads*nv*no*sizeof(double);
	// pGamma_
	reqMem += maxThreads*nv*nl*sizeof(double);
	// pTrans_
	reqMem += nv*nl*sizeof(double);

	// +1 for good measure
	return reqMem+1;

    }

    CPUPmatGammaFunctor(const size_t nf, const size_t no, const size_t nv, const size_t nl,
			const double *ea, const double *ev,
			double &energy,
			MapMatrixXd &pmat,
			std::vector<double *> &pGammaLong,
			double *gamma_rs,
			double *data_pqm1):
	nf_(nf),no_(no),nv_(nv),nl_(nl),na_(no_-nf_),
	ea_(ea),ev_(ev),
	energy_(energy),
	pmat_(pmat),
	pGammaLong_(pGammaLong),
	gamma_rs_(gamma_rs),
	data_pqm1_(data_pqm1),
	max_threads_(omp_get_max_threads() )
    {

	pEriLong_ = new double[na_*nv_*nv];
	pTijLong_ = new double[na_*nv_*nv];
	for(int i = 0; i <na_; i++){
	    pVecEriLong_.push_back(&pEriLong_[i*nv_*nv_]);
	    pVecTijLong_.push_back(&pTijLong_[i*nv_*nv_]);
	}//i

	pEri_ = new double[max_threads_*nv_*nv];
	pPmatVirt_ = new double[max_threads_*nv_*nv];

	pPmatOcc_ = new double[max_threads_*na_*na_];

	pGamma_ = new double[max_threads_*nv_*nl_];
	pTrans_ = new double[nv_*nl_];
	
	for(int i = 0; i < max_threads_; i++){
	    
	    pVecEri_.push_back(&pEri_[i*nv_*nv_]);
	    
	    pVecPmatVirt_.push_back(&pPmatVirt_[i*nv_*nv_]);
	    memset(pVecPmatVirt_[i], 0, nv_*nv_*sizeof(double) );
	    
	    pVecPmatOcc_.push_back(&pPmatOcc_[i*na_*na_]);
	    memset(pVecPmatOcc_[i], 0, na_*na_*sizeof(double) );
	    
	    pVecGamma_.push_back(&pGamma_[i*nv_*nl_]);
	    
	    threadEnergy_.push_back( (double)0.0 );
	    
	}//i
	
    }//constructor
    
    ~CPUPmatGammaFunctor(){
	
	for(int i = 0; i < max_threads_; i++){
	    
	    MapMatrixXd pmatOcc(pVecPmatOcc_[i],no_-nf_,no_-nf_);
	    pmat_.block(nf_,nf_,no_-nf_,no_-nf_) += pmatOcc;
	    
	    MapMatrixXd pmatVirt(pVecPmatVirt_[i],nv_,nv_);
	    pmat_.block(no_,no_,nv_,nv_) += pmatVirt;
	    
	    energy_ += threadEnergy_[i];
	    
	}//i
	
	delete [] pEri_;
	delete [] pEriLong_;
	delete [] pTijLong_;
	delete [] pPmatVirt_;
	delete [] pPmatOcc_;
	delete [] pGamma_;
	delete [] pTrans_;
	
    }//destructor
    
    void operator()(const double *ptr_a1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
		    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
		    const double *ptr_a2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
		    size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols);
    
}; // CPUPmatGammaFunctor






/*
  @brief W[III] and MO lagrangian
  @author LBR
  @detail this routine can be used to build W[III] (energy weighted density) 
  matrix and the MO lagrangian using J/K matricies with factorization of 
  the appropriate matrix
  @param N number of OBS functions
  @param no number of occupied MOs
  @param nv number of virtual MOs
  @param nl number of ABS functions
  @param pCa pointer to occupied LCAO coefficients
  @param pCv pointer to virtual LCAO coefficients
  @param LM1_MN_SYM array pointer to unique (mu nu|L) integrals
  @param pmat relaxed density matrix
  @param pe parallel environment
*/
    class matrix_factor : public ::cchem::rimp2_gradient::detail::JK {

	const double *pCa_;
	const double *pCv_;

	size_t N_;
	const size_t no_;
	const size_t nv_;
	const size_t nl_;

	Array<double> *LM1_MN_SYM_;
	MapMatrixXd &densityMatrix_;

	double *pScrNno_;

	int Npi_,Nni_;

	Parallel &pe_;
	::cchem::rimp2::device *pGPU_;


    public:
	matrix_factor(size_t N, const size_t no,
		      const size_t nv, const size_t nl,
		      const double *pCa, const double *pCv,

		      Array<double> *LM1_MN_SYM,
		      MapMatrixXd &densityMatrix,
		      Parallel &pe,
		      ::cchem::rimp2::device *pGPU = NULL):
	    ::cchem::rimp2_gradient::detail::JK(nl,N),
	    pCa_(pCa),pCv_(pCv),
	    N_(N), no_(no), nv_(nv),nl_(nl),	    
	    LM1_MN_SYM_(LM1_MN_SYM),densityMatrix_(densityMatrix),pe_(pe),pGPU_(pGPU)
	{

	    pScrNno_ = new double[N_*no_];

	    //JK inheritence allocates a J and K, but we need another of each
	    pJ_.push_back( new double[N_*N_] );
	    pK_.push_back( new double[N_*N_] );

	    for(int i = 0; i < 2; i++){
		pDAux_.push_back( new double[ nl_ ] );
		pDmn_.push_back( new double[ mnSize_ ] );
		pD_.push_back( new double[ N_*N_ ] );
		pLeft_.push_back( new double[ (no_+nv_)*N_ ] );
		pRight_.push_back( new double[ (no_+nv_)*N_ ] );
	    }//i

	};

	~matrix_factor(){

	    delete [] pScrNno_;

	    delete [] pJ_[1];
	    delete [] pK_[1];
	    for(int i = 0; i < 2; i++){
		delete [] pDAux_[i];
		delete [] pDmn_[i];
		delete [] pD_[i];
		delete [] pLeft_[i];
		delete [] pRight_[i];
	    }//i


	}//~matrix_factor

	template<typename T>
	void print_matrix(T & mat);

	void buildLagMO(MapMatrixXd &lag_mo, MapMatrixXd &ap);

	void buildWmat3(MapMatrixXd &wmat3, MapMatrixXd &ap);

	void buildJK();

	void pseudoSquareRoot(const double &cutoff);

};

}//namespace detail
}//namespace rimp2_gradient
}//namespace cchem










//#include <device.hpp>

// namespace cchem{
// namespace rimp2_gradient{
// namespace detail{


// typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;



// class rimp2_lagrangian : boost::noncopyable {


// 	Thread thread_;
// 	size_t nd_,ns_,nv_,nvs_,nl_,no_;
// 	Array<double> *CIAQ_;
// 	Array<double> *BARE_IJ_;
// 	Array<double> *BARE_AB_;


// public:

// 	//	static
// 	//	boost::array<size_t,4> memory(const size_t &nv, const size_t &no, const size_t &nl){
// 	//		boost::array<size_t,4> memory;
// 	//		memory[0] = nv*nl;
// 	//		memory[1] = no*nl;
// 	//		memory[2] = nv*nl;
// 	//		memory[3] = nv*nv; //worse case
// 	//		return memory;
// 	//	}

// 	rimp2_lagrangian(const size_t &nd, const size_t &ns, const size_t &nv, const size_t &nl,
// 			 Array<double> *CIAQ,Array<double> *BARE_IJ,Array<double> *BARE_AB)
			 
// :nd_(nd),ns_(ns),nv_(nv),nl_(nl),nvs_(nv+ns),no_(nd+ns),
//  CIAQ_(CIAQ), BARE_IJ_(BARE_IJ), BARE_AB_(BARE_AB)
// {
// 		//		BOOST_AUTO(memory, this->memory(nv_, no_, nl_));
// 		//		for (int i = 0; i < memory.size(); ++i) {
// 		//			thread_.data[i] = thread_.malloc(memory[i]);
// 		//		}
// };

// 	~rimp2_lagrangian(){thread_.free();};

// 	double *operator()(int i){return thread_.data[i];}

// 	void vooo(const MapMatrixXd &eri_vooo, const MapMatrixXd &pij, MapMatrixXd &lag, const int i, const int j);

// 	void vvvo(const MapMatrixXd &eri_vvvo, const MapMatrixXd &pab, MapMatrixXd &lag, const int i, const int b);

// 	//	void set_up_block_info();

// 	void build_pmat_contributions(const MapMatrixXd &pmat, MapMatrixXd &lag_mo);


// };//rimp2_lagrangian


// }//namespace detail
// }//namespace rimp2_gradient
// }//namespace cchem






// namespace cchem {
// namespace rimp2_gradient {
// namespace detail {


// struct lag_vovo_functor{
// private:
// 	size_t no_,nv_;
// 	MapMatrixXd &u_;
// 	MapMatrixXd &unxt_;
// public:
// 	lag_vovo_functor(size_t no, size_t nv, MapMatrixXd &u, MapMatrixXd &unxt):
// 		no_(no),nv_(nv), u_(u), unxt_(unxt) {}
// 	void operator() (MapMatrixXd &eri_vooo,
// 			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 			MapMatrixXd &bjk, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
// 			bool &loop_restriction);
// }; //lag_vovo_functor






// struct lag_vvoo_functor{
// private:
// 	size_t no_,nv_;
// 	MapMatrixXd &u_;
// 	MapMatrixXd &unxt_;
// public:
// 	lag_vvoo_functor(size_t no, size_t nv, MapMatrixXd &u, MapMatrixXd &unxt):
// 		no_(no),nv_(nv), u_(u), unxt_(unxt) {}
// 	void operator() (MapMatrixXd &eri_vooo,
// 			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 			MapMatrixXd &bjk, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
// 			bool &loop_restriction);
// }; //lag_vvoo_functor

// struct lag_vooo_functor{
// private:
// 	size_t no_,nv_;
// 	MapMatrixXd &pmat_;
// 	MapMatrixXd &lag_mo_;
// public:
// 	lag_vooo_functor(size_t no, size_t nv, MapMatrixXd &pmat, MapMatrixXd &lag_mo):
// 		no_(no),nv_(nv), pmat_(pmat), lag_mo_(lag_mo) {}
// 	void operator() (MapMatrixXd &eri_vooo,
// 			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 			MapMatrixXd &bjk, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
// 			bool &loop_restriction);
// }; //lag_vooo_functor

// struct lag_vvvo_functor{
// private:
// 	size_t no_,nv_;
// 	MapMatrixXd &pmat_;
// 	MapMatrixXd &lag_mo_;
// public:
// 	lag_vvvo_functor(size_t no, size_t nv, MapMatrixXd &pmat, MapMatrixXd &lag_mo):
// 		no_(no),nv_(nv), pmat_(pmat), lag_mo_(lag_mo) {}
// 	void operator() (MapMatrixXd &eri_vvvo,
// 			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
// 			MapMatrixXd &bba, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
// 			bool &loop_restriction);
// }; //lag_vvvo_functor




// }//namespace detail 
// }//namespace rimp2_gradient
// }//namespace cchem






#endif /* SRC_RI_MP2_RI_LAGRANGIAN_HPP_ */
