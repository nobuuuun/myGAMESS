/*
 * metric.hpp
 *
 *  Created on: Nov 4, 2015
 *      Author: luke
 */

#ifndef SRC_RI_MP2_METRIC_HPP_
#define SRC_RI_MP2_METRIC_HPP_


#include "thread.hpp"
#include <Eigen/Dense>
#include <basis/basis.hpp>
//#include <gradient-helpers.hpp>



namespace cchem{
namespace rimp2_gradient{

namespace detail{

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef ::basis::Basis Basis;
using cchem::Thread;


/*
  @brief build coulomb metric
  @author LBR
  @detail either use cholesky or pseudoinversion
  @param auxbasis ABS
  @param pe parallel environment
  @param spherical flag to indicate if ABS is spherical
  @param nl number of ABS functions
  @param data_twoc two-center ERI matrix
  @param data_pqm1 inverted coulomb matrix
*/
class Metric : boost::noncopyable {

	Thread thread_;
	boost::reference_wrapper<const Basis> &auxbasis_;
	Parallel &pe_;
	size_t nl_;
	int spherical_;
	double * data_twoc_;
	double * data_pqm1_;

public:

	Metric(boost::reference_wrapper<const Basis> auxbasis, Parallel &pe, int spherical,
			size_t nl, double *data_twoc, double *data_pqm1) :
				pe_(pe),nl_(nl), spherical_(spherical),
				auxbasis_(auxbasis),data_twoc_(data_twoc), data_pqm1_(data_pqm1){};

	void BuildTwoCenterEri(std::vector<MapMatrixXd> &sph_c_);

	void DecomposeLDLT();

	void DecomposeLLT();

	//this is for L -> L-1 then  (L-1).transpose().L to form pqm1.
	void InvertL();

	void Lm1_to_pqm1(){
	    MapMatrixXd twoc_pqm1(data_pqm1_, nl_, nl_);
	    twoc_pqm1 = twoc_pqm1.transpose()*twoc_pqm1;
	}

/*
   @brief create (P|Q)^-1/2 via pseudo-inversion with eigen decomposition
   @author LBR
   @param tol [in] tolerance to keep eigen vectors
   @param pPqm12 [out] sqrt inverse metric pointer
   @param pTwoC [in] two-center two-electron matrix pointer
   @param doPrint [in] indication to which process is master
*/
    void pseudoInvertedSqrt(const double &tol, double *pPqm12, double *pTwoC, bool doPrint);

/*
   @brief create L^-1 or (P|Q)^-1 via inversion of choleksy decomposed matrix
   @author LBR
   @param data_pqm1 [out] inverse coulomb metric
   @param data_twoc [in/out] two-center two-electron matrix -> transformed into L^-1
   @param doPrint [in] indication to which process is master
*/
    void choleskyInversion(double *pqm1, double *twoc_eri, bool doPrint);

};//Metric

}//Namespace detail

}//namespace rimp2_gradient
}//namespace cchem


#endif /* SRC_RI_MP2_METRIC_HPP_ */
