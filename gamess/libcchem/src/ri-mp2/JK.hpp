/*
 * JK.hpp
 *
 *  Created on: Feb 4, 2016
 *      Author: luke
 */

#ifndef SRC_RI_MP2_JK_HPP_
#define SRC_RI_MP2_JK_HPP_

#include "god.hpp"

#include <Eigen/Core>

#include "core/wavefunction.hpp"
#include "thread.hpp"

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef Eigen::Map<Eigen::VectorXd,Eigen::AutoAlign> MapVectorXd;

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;



#include <math.hpp>
namespace cchem{
namespace rimp2_gradient{
namespace detail{




/*
  @brief JK matrix base class
  @author LBR
  @detail build JK matrices with CPU only
  @param 
*/

class JK : boost::noncopyable{

    //auxiliary basis size
    const size_t nl_;
    //OBS size
    const size_t N_;

protected:

    //N*(N+1)/2
    const size_t mnSize_;
    // Pseudo-occupied C matrices, left side
    std::vector<double *> pLeft_;
    // Pseudo-occupied C matrices, right side
    std::vector<double *> pRight_;
    // Pseudo-density matrices (N_,N_)
    std::vector<double *> pD_;
    // J matrices (N_,N_)
    std::vector<double *> pJ_;
    // K matrices (N_,N_)
    std::vector<double *> pK_;


    std::vector<double *> pDmn_;
    std::vector<double *> pDAux_;
    
    //temporary storage for N*(N+1)/2 J matrix
    double *pJTemp_;
    
    
    
public:
    JK(const size_t nl, size_t N):
	nl_(nl),N_(N), mnSize_( (N_*(N_+1)/2) )
    {
	//allocate some memory
	pJTemp_ = new double[mnSize_];
	pJ_.push_back( new double[N_*N_] );
	pK_.push_back( new double[N_*N_] );
    }
    
    ~JK(){
	delete [] pJTemp_;
	delete [] pJ_[0];
	delete [] pK_[0];
    }

// the following two routine can speed up z-vector iterations since the left side is the same (i.e. only the right side changes)
//        -BUT this would require storing the full left side, which can be expensive. so abandon this...
//	//build left factor of K matrix
//	void buildKLeft(size_t N, Eigen::MatrixXd &Qmn, double *lEl);
//
//	//build right factor of K matrix then K matrix
//	void buildKRight(size_t N, Eigen::MatrixXd &Qmn, double *lEl);

	// //build J matrix
	// void buildJMatrix(size_t N, double *pQmnSym, const size_t naux); //, Eigen::VectorXd &J_temp_);

	// //build K matrix
	// void buildKMatrix(size_t N, int same, double *pQmnSym, const size_t naux,
	// 		std::vector<int> &dims);

}; //JK









class buildJKMatrixFunctor : boost::noncopyable {

    const size_t mnSize_,N_;
    double *pJTemp_;
    std::vector<double *> &pD_, &pJ_, &pDmn_, &pDAux_,&pK_,&pLeft_,&pRight_;
    std::vector<int> &dims_;
    double *pExpLeft_;
    double *pExpRight_;
    size_t maxNAB_;
    size_t maxNocc_;
    
public:
    
    buildJKMatrixFunctor(const size_t mnSize, const size_t N,
			 double *pJTemp, std::vector<double *> &pD, 
			 std::vector<double *> &pJ, std::vector<double *> &pDmn, 
			 std::vector<double *> &pDAux, std::vector<double *> &pK, 
			 std::vector<double *> &pLeft, std::vector<double *> &pRight, 
			 std::vector<int> &dims, size_t maxNAB):
	mnSize_(mnSize), N_(N), pJTemp_(pJTemp), pD_(pD), pJ_(pJ), 
	pDmn_(pDmn), pDAux_(pDAux),pK_(pK),pLeft_(pLeft),pRight_(pRight),dims_(dims),
	maxNAB_(maxNAB),maxNocc_(std::max(dims_[0],dims_[1])),
	pExpLeft_(NULL), pExpRight_(NULL)
    {
	pExpLeft_ = new double[maxNAB_*maxNocc_*N_];
	pExpRight_ = new double[maxNAB_*maxNocc_*N_];
    };
    
    ~buildJKMatrixFunctor(){
	delete [] pExpLeft_ ;
	delete [] pExpRight_ ;
    };
    
    void operator()(double *pQmnSym, 
		    size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
		    size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols);
    
};//buildJKMatrixFunctor
    
    
}//namespace detail
}//namespace rimp2_gradient
}//namespace cchem



#endif /* SRC_RI_MP2_JK_HPP_ */
