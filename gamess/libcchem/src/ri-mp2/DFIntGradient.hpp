/*
 * DFIntGradient.hpp
 *
 *  Created on: Feb 12, 2016
 *      Author: luke
 */

#ifndef SRC_RI_MP2_DFIntGRADIENT_HPP_
#define SRC_RI_MP2_DFIntGRADIENT_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/ref.hpp>
#include "array/hdf5.hpp"







#include "core/wavefunction.hpp"
#include "runtime.hpp"








#include <Eigen/Dense>
#include <ri-async.hpp>
#include "gradient-helpers.hpp"

namespace cchem{
namespace rimp2_gradient{
namespace DFIntGradient{

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;
typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef const Eigen::Map<const Eigen::MatrixXd,Eigen::AutoAlign> ConstMapMatrixXd;

class gradient {

	//occupied MOs
	boost::reference_wrapper<const Matrix> &Ca_;
	//virtual MOs
	boost::reference_wrapper<const Matrix> &Cv_;
	//SCF density (MO basis)
	MapMatrixXd &Pscf_;
	//total density : P(2) (MO basis)
	Eigen::MatrixXd &Pmo_;
	//number of basis functions, occupied MOs, virtuals MOs, number of auxiliary basis functions
	size_t N_,no_,nv_,nl_;


	//pointer to (ao ao|aux) integral storage
	Array<double> *BARE_MN_SYM_;

	//inverse coulomb metric (P|Q)^-1
	MapMatrixXd &pqm1_;
	//auxiliary basis set
	boost::reference_wrapper<const Basis> auxbasis_;
	//orbital basis set
	boost::reference_wrapper<const Basis> basis_;
	//spherical to cartesian transformation matrix
	std::vector<MapMatrixXd> sph_c_;
	//spherical to cartesian transformation matrix
	std::vector<double *> &pSphTrans_;
	//energy gradient
	Eigen::MatrixXd &eg_;

	//number of atoms
	int &natoms_;
	//are we using spherical auxiliary functions
	const int &spherical_;

	//parallel environment
	Parallel &pe_;

	Runtime &rt_;

	//positive factor
	std::vector<Eigen::MatrixXd> posFactor_;

	//negative factor
	std::vector<Eigen::MatrixXd> negFactor_;

    //(naux,naux)
	std::vector<Eigen::MatrixXd> auxForce_;


	//reference density in AO basis
	std::vector<Eigen::MatrixXd> D2AO_;

	//full density in AO basis
	std::vector<Eigen::MatrixXd> P2AO_;


    //3-center ERIs contacted with reference density
    std::vector<Eigen::VectorXd> DAux_;    

    //3-center ERIs contacted with correlated density
    std::vector<Eigen::VectorXd> PAux_;    

	//number of positive, negative eigenvalues for factorization
	int Npi_;
	int Nni_;


public:
	gradient(boost::reference_wrapper<const Matrix> Ca,
			boost::reference_wrapper<const Matrix> Cv,
			MapMatrixXd &Pscf,
			Eigen::MatrixXd &Pmo,
			size_t N,
			size_t no,
			size_t nv,
			size_t nl,
			Array<double> *BARE_MN_SYM,
			MapMatrixXd &pqm1,
			boost::reference_wrapper<const Basis> auxbasis,
			boost::reference_wrapper<const Basis> basis,
			std::vector<MapMatrixXd> sph_c,
			std::vector<double *> &pSphTrans,
			Eigen::MatrixXd &eg,
			int natoms,
			const int &spherical,
			Parallel &pe,
			Runtime &rt
	)
 :Ca_(Ca), Cv_(Cv), Pscf_(Pscf),Pmo_(Pmo),N_(N),no_(no),nv_(nv),nl_(nl),
  BARE_MN_SYM_(BARE_MN_SYM),pqm1_(pqm1),auxbasis_(auxbasis),basis_(basis),sph_c_(sph_c),
  pSphTrans_(pSphTrans),eg_(eg),
  natoms_(natoms),spherical_(spherical),pe_(pe),rt_(rt)
    {
	Npi_ = 0;
	Nni_ = 0;
    };

    ~gradient(){};

    void compute_DFIntGradient();
    
    void build_three_center_terms(Array<double> * BARE_PN_);

    //decompose P(2) into +/- factors and transform
    // e.g. (MO,MO) -> (MO,+/-) -> (AO,+/-)
    void setup();

    void setupAsyncMNAuxAccess(size_t &NMB, std::vector <std::pair<size_t,size_t> >&AsShWork,
			       std::vector <std::pair<size_t,size_t> >&AsBasWork,
			       size_t &MaxBuffer);
    
    void TwoC_DERI_contribution();

    void ThreeC_DERI_contribution(Array<double> * BARE_PN_, Array<double> * MNAUX_FORCE_);







    struct DAuxPAuxFunctor : boost::noncopyable{
	private:
	//reference density 
	const double *pDmn_;
	//correlated density
	const double *pPmn_;
	Eigen::VectorXd &DAux_;
	Eigen::VectorXd &PAux_;
	const size_t N_;
	
    public:
		DAuxPAuxFunctor(const double *pDmn, const double *pPmn,
				Eigen::VectorXd &DAux, Eigen::VectorXd &PAux, const size_t N):
					pDmn_(pDmn), pPmn_(pPmn), DAux_(DAux), PAux_(PAux), N_(N){}


	    void operator()(double *ao_ints, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
			    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols)
	    {
		
		const size_t mnSize = N_*(N_+1)/2;
		
#pragma omp for schedule(dynamic)
		for(size_t i = A1_BCStart; i < A1_BCStop; i++){
		    const double *ptr_temp = ao_ints + mnSize*(i-A1_BCStart);
		    
		    DAux_(i) = cblas_ddot(mnSize,pDmn_, 1,ptr_temp,1);
		    PAux_(i) = cblas_ddot(mnSize,pPmn_, 1,ptr_temp,1);
		    
		}//i
		
	    }; //operator()
	    
	    
	}; //DAuxPAuxFunctor





























struct PNFunctor : boost::noncopyable{
private:
    Eigen::MatrixXd &Pos_;
    Eigen::MatrixXd &Neg_;
    const size_t N_, no_;
    const int Npi_, Nni_;
    ConstMapMatrixXd &occ_coeff_mat_;
    
    double *ptr_scratch_Ain_;
    std::vector<double *> pScrAin_;
public:
    PNFunctor(Eigen::MatrixXd &Pos, Eigen::MatrixXd &Neg,
	      const size_t N, const size_t no, const int Npi, const int Nni,ConstMapMatrixXd &occ_coeff_mat):
	Pos_(Pos), Neg_(Neg),N_(N), no_(no), Npi_(Npi), Nni_(Nni),occ_coeff_mat_(occ_coeff_mat){
	
	const int max_threads = omp_get_max_threads();
		
	ptr_scratch_Ain_ = new double[max_threads*no_*N_];
	
	double * temp_ptr = NULL;
	for(int i = 0; i < max_threads; i++)
	    pScrAin_.push_back( &ptr_scratch_Ain_[i*no_*N_] );
    }//constructor
    
    ~PNFunctor(){
	delete [] ptr_scratch_Ain_;
    }//~
    
    
    void operator()(double *ao_ints, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
		    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
		    double *result, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
		    size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols)
    {
	const size_t mnSize = N_*(N_+1)/2;
	double *pScratchAin = pScrAin_[omp_get_thread_num()];
	const double *pPos_ = Pos_.data();
	const double *pNeg_ = Neg_.data();
	const double *pOccCoeffMat_ = occ_coeff_mat_.data();
	
	double *pIntScr = new double [N_*N_];
#pragma omp for schedule(dynamic)
	for(size_t i = A1_BCStart; i < A1_BCStop; i++){
	    
	    double *ptr_ao = &ao_ints[(i-A1_BCStart)*mnSize];
	    
	    for(size_t mu = 0; mu < N_; mu++){
		const size_t nInts = mu +1;
		const double *pGetPos = &ptr_ao[ mu*(mu+1)/2 ];
		double *pPutPos = &pIntScr[ N_*mu ];
		memcpy(pPutPos, pGetPos, nInts*sizeof(double) );
		
	    }//mu
	    
	    
	    
	    cblas_dsymm(CblasColMajor,CblasRight, CblasUpper, //Lower,
			no_, N_,
			1.0, pIntScr, N_,
			pOccCoeffMat_, no_,
			0.0, pScratchAin, no_);
	    
	    double *ptr_la = &result[(i-A1_BCStart)*no_*(Npi_+Nni_)];
	    double *ptr_ra = &result[(i-A1_BCStart)*no_*(Npi_+Nni_) + Npi_*no_];
	    
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
			Npi_,no_,N_,
			1.0, pPos_, Npi_,
			pScratchAin,no_,
			0.0, ptr_la, Npi_);
	    
	    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
			Nni_,no_,N_,
			1.0, pNeg_, Nni_,
			pScratchAin,no_,
			0.0, ptr_ra, Nni_);
	    
	}//i
	
	delete [] pIntScr;
    }; //operator()
}; //PNFunctor







































    struct auxForceFunctor : boost::noncopyable{
    private:
    	const size_t N_, no_, nl_;
    	const int Npi_, Nni_;
    	Eigen::MatrixXd &auxForce_;
	
    	//number of CPU threads allowed
    	int maxThreads_;

    	double *threadData_;
    	//	std::vector<double *> vecThreadData_;

    public:
    	auxForceFunctor(const size_t N, const size_t no, const size_t nl, 
    		 const int Npi, const int Nni, Eigen::MatrixXd &V):
    	    N_(N), no_(no), nl_(nl), Npi_(Npi), Nni_(Nni), auxForce_(V)
    	{

	    // maxThreads_ = omp_get_max_threads();    
	    // maxThreads_ = 1;
    	    // threadData_ = new double[nl_*nl_];
    	    // memset( threadData_, 0, nl_*nl_*sizeof(double) );
    	    // for(int i = 0; i < max_threads_; i++){
    	    // 	vecThreadData_.pushBack( &threadData[nl_*nl_*i] );
    	    // }//i
	    


	    
    	}
	

    	~auxForceFunctor()
    	{

    	    //V_.setZero();
    	    // if(nl_*nl > (2**31-1) )
    	    // 	std::cout << "there may be an integer problem at" 
    	    // 		  << __LINE__  << " of file " 
    	    // 		  << __FILE__ <<std::endl; 

    	    // for(int i = 1; i < maxThreads_; i++){
    	    // 	cblas_daxpy(nl_*nl_, 1.0,
    	    // 		    vecThreadData_[i],1,
    	    // 		    vecThreadData_[0],1);
    	    // }//i

    	    //std::memcpy(V_.data(), threadData_, nl_*nl_*sizeof(double) );

    	    //delete [] threadData_;
    	}
	
    	void operator()	(const double *pFactor1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
    			 size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
    			 const double *pFactor2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
    			 size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols)
    		{

		    //#pragma omp parallel
		    {

                    const size_t maxThreads = omp_get_num_threads();
    		    const int tid = omp_get_thread_num();
    		    //			const int cols1 = A1_BCStop - A1_BCStart;
    		    const int cols2 = A2_BCStop - A2_BCStart;

    		    //get chunk size (for cols1 only)
    		    //size_t chunkSize = (A1_BCStop - A1_BCStart)/maxThreads_;
    		    size_t chunkSize = (A1_BCStop - A1_BCStart)/maxThreads;

    		    const size_t posChunkStart = tid*chunkSize*no_*(Npi_+Nni_);
    		    const size_t negChunkStart = tid*chunkSize*no_*(Npi_+Nni_) + no_*Npi_;

    		    const size_t VStart = tid*chunkSize + A1_BCStart;

    		    //assign remainder (if any) to last thread
    		    if(tid == maxThreads-1 )
			//chunkSize += (A1_BCStop - A1_BCStart)%maxThreads_;
    		    	chunkSize += (A1_BCStop - A1_BCStart)%maxThreads;

    		    //postive part
    		    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
    		    		chunkSize, cols2, no_*Npi_,
    		    		2.0, &pFactor1[posChunkStart], no_*(Npi_+Nni_),
    		    		pFactor2, no_*(Npi_+Nni_),
    		    		0.0, &auxForce_.data()[VStart + A2_BCStart*nl_], nl_);

    		    //negative part
    		    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
    		    		chunkSize, cols2, no_*Nni_,
    		    		-2.0, &pFactor1[negChunkStart], no_*(Npi_+Nni_),
    		    		&pFactor2[no_*Npi_], no_*(Npi_+Nni_),
    		    		1.0, &auxForce_.data()[VStart + A2_BCStart*nl_], nl_);

		    }//omp parallel

    		}; //operator()

    	}; //auxForceFunctor






    struct MNAuxForceFuctor: boost::noncopyable {

private:
    const double *pPosFactor_;
    const double *pNegFactor_;
    const size_t N_, no_;
    const int Npi_, Nni_;
    const double *pOccCoeff_;
    
    double *pScrMNAux_;
    std::vector<double *> scrMNAuxForceVec_;
public:
    
	MNAuxForceFuctor(const double *pPosFactor, const double *pNegFactor,
			 const size_t N, const size_t no, const int Npi, const int Nni, const double *pOccCoeff):
	    pPosFactor_(pPosFactor), pNegFactor_(pNegFactor),N_(N), no_(no), Npi_(Npi), Nni_(Nni),pOccCoeff_(pOccCoeff){
	    
	int max_threads = omp_get_max_threads();
	
	pScrMNAux_ = new double[max_threads*no_*N_];
	for(int i = 0; i < max_threads; i++)
	    scrMNAuxForceVec_.push_back(pScrMNAux_ + i*no_*N_);
	
    }//constructor
    
    ~MNAuxForceFuctor(){
	delete [] pScrMNAux_;
    }//~MNAuxForceFuctor
    
    void operator()(double *ao_ints, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
		    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
		    double *result, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
		    size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols)
    {
	
	double *pScrMNAux = scrMNAuxForceVec_[omp_get_thread_num()];

#pragma omp for schedule(dynamic)
	for(size_t i = A1_BCStart; i < A1_BCStop; i++){
	    
	    //(mu "+" |Q)(Q|P)-1
	    const double *pMuPos = &ao_ints[(i-A1_BCStart)*no_*(Npi_+Nni_)];
	    //(mu "-" |Q)(Q|P)-1
	    const double *pMuNeg = &ao_ints[(i-A1_BCStart)*no_*(Npi_+Nni_)+no_*Npi_];
    
	    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
	    		no_,N_,Npi_,
	    		1.0, pMuPos, Npi_,
	    		pPosFactor_,Npi_,
			0.0, pScrMNAux, no_);
	    
	    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
	    		no_,N_,Nni_,
	    		-1.0, pMuNeg, Nni_,
	    		pNegFactor_,Nni_,
	    		1.0, pScrMNAux, no_);
	    
	    //AO->MO
	    double *pMNAuxForce = &result[(i-A1_BCStart)*N_*N_];
	    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
	    		N_,N_,no_,
	    		2.0, pOccCoeff_, no_,
	    		pScrMNAux, no_,
	    		0.0, pMNAuxForce, N_);
	    
	}//i
	
    }; //operator()
    
}; //MNAuxForceFuctor





	template<typename T>
	void print_matrix(T & mat){

		std::cout << mat.rows() << "x" << mat.cols() << std::endl;
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

};

}//namespace DFIntGradient
}//namespace rimp2_gradient
}//namespace cchem




#endif /* SRC_RI_MP2_DFIntGRADIENT_HPP_ */
