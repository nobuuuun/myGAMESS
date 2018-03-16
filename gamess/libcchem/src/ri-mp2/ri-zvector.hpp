/*
 * ri-zvector.hpp
 *
 *  Created on: Nov 5, 2015
 *      Author: luke
 */

#ifndef SRC_RI_MP2_RI_ZVECTOR_HPP_
#define SRC_RI_MP2_RI_ZVECTOR_HPP_



#include "thread.hpp"
#include <boost/noncopyable.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/ref.hpp>



#include <ri-async.hpp>

//#include <ri-mp2/ri-lagrangian.hpp>
#include <ri-lagrangian.hpp>
//#include <ri-integrals.hpp>
//#include "ri-mp2/ri-integrals.hpp"

namespace cchem{
namespace rimp2_gradient{
namespace detail{

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::vector_range<const vector_type> vector_range;

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef Eigen::Map<Eigen::VectorXd,Eigen::AutoAlign> MapVectorXd;
typedef Eigen::Map<Eigen::VectorXi,Eigen::AutoAlign> MapVectorXi;

using cchem::Thread;

    //#include <math.hpp>
#include <JK.hpp>

#include "god.hpp"


/*
  @brief z-vector base class
  @author LBR
  @detail plug in appropriate electronic hessian
  @param nd number of doubly occupied MOs
  @param ns number of singly occpied MOs
  @param nv number of virtual MOs
  @param maxitc maxium zvector iterations
  @param uconv convergence criteria
  @param small convergence criteria
  @param ea active orbital energies
  @param ev virtual orbtital energies
  @param pe parallel environment  
*/


class ZVector : boost::noncopyable {

protected:
    Thread thread_;
    size_t nd_,ns_,nv_,no_,maxitc_;
    
    double uconv_,small_;
    boost::reference_wrapper<const vector_range> ea_;
    boost::reference_wrapper<const vector_range> ev_;
    
    //parallel environment
    Parallel &pe_;
    
	static
	boost::array<size_t,5> memory(const size_t &nv, const size_t &no, // const size_t &nl,
			const size_t &maxitc){
		boost::array<size_t,5> memory;

		memory[0] = maxitc; //B
		memory[1] = maxitc*maxitc; //UAU
		memory[2] = maxitc; //CC
		memory[3] = maxitc; //UU
		memory[4] = maxitc*maxitc; //ALPHA

		return memory;
	}//memory

public:

    ZVector(const size_t &nd, const size_t &ns, const size_t &nv, //const size_t &nl,
	    const size_t &maxitc,double uconv, double small,
	    boost::reference_wrapper<const vector_range> ea,
	    boost::reference_wrapper<const vector_range> ev,
	    Parallel &pe) :
	nd_(nd),ns_(ns),nv_(nv),no_(nd+ns),
	maxitc_(maxitc),uconv_(uconv),small_(small),
	ea_(ea),ev_(ev),pe_(pe)
    {
	BOOST_AUTO(memory, this->memory(nv_, no_, maxitc_));
	for (int i = 0; i < memory.size(); ++i) {
	    thread_.data[i] = thread_.malloc(memory[i]);
	}//i
    }//ZVector

	virtual ~ZVector(){ thread_.free();}// std::cout << "deleting ZVector" << std::endl; }

	void sym_eig(MapMatrixXd &lag, const vector_range &ea, const vector_range &ev, std::vector<MapMatrixXd> &u);

	void sym_eig(const vector_range &ea,const vector_range &ev, MapMatrixXd &unxt, int iter);

	void build_uau(int iter, std::vector<MapMatrixXd> &u, MapMatrixXd &unxt);

	void form_alpha(int iter, std::vector<MapMatrixXd> &u);

	void lu(int iter, MapVectorXi &ipivot);

	void lus(const int iter, MapVectorXi &ipivot );

	void form_new_solution(int iter, std::vector<MapMatrixXd> &u, MapMatrixXd &rhs);

	void update_unxt(int iter, MapMatrixXd &prhs, MapMatrixXd &unxt, std::vector<MapMatrixXd> &u);

	void solve(MapMatrixXd &lag_mo, MapMatrixXd &pmat, MapMatrixXd &zai);

	virtual void build_orbital_hessian(MapMatrixXd &u, MapMatrixXd &unxt) = 0;

};//zvector





/*
  @brief RI z-vector derived class
  @author LBR
  @detail this builds the orbital hessian without JK matrices
  @param nd number of doubly occupied MOs
  @param ns number of singly occpied MOs
  @param nv number of virtual MOs
  @param nl number of ABS functions
  @param maxitc maxium zvector iterations
  @param uconv convergence criteria
  @param small convergence criteria
  @param ea active orbital energies
  @param ev virtual orbtital energies
  @param BARE_IA array pointer to (ia|L) integrals
  @param BARE_AB array pointer to (ab|L) integrals
  @param CIAQ array pointer to C_IA^Q coefficients
  @param CIJQ array pointer to C_IJ^Q coefficients
  @param pe parallel environment  
*/

// class RIZVector : public ::cchem::rimp2_gradient::detail::ZVector {

// 	const size_t nl_;
// 	Array<double> *BARE_IA_;
// 	Array<double> *BARE_AB_;
// 	Array<double> *CIAQ_;
// 	Array<double> *CIJQ_;

// 	std::vector < std::pair <std::vector <size_t>, std::vector <size_t> > >	CIAReadPos_,ABReadPos_,CReadPos_;
// 	std::vector < std::pair < size_t, size_t> > CIABlocks_,ABBlocks_,CBlocks_;
// 	size_t MaxNCIA_ ,MaxNAB_, MaxNC_, NCIAReads_, NABReads_, NCReads_;

// 	::cchem::rimp2_gradient::detail::RIAsync test_lag;
// 	::cchem::rimp2_gradient::detail::RIAsync test_lag2;

// public:

// 	RIZVector(const size_t &nd, const size_t &ns, const size_t &nv, const size_t &nl,
// 			const size_t &maxitc, double uconv, double small,
// 			boost::reference_wrapper<const vector_range> ea,boost::reference_wrapper<const vector_range> ev,
// 			Array<double> *BARE_IA, Array<double> *BARE_AB, Array<double> *CIAQ, Array<double> *CIJQ,
// 		  const size_t &NMB, Parallel &pe):
// 	    ::cchem::rimp2_gradient::detail::ZVector(nd, ns, nv, maxitc,uconv,small,ea,ev,pe),
// 				 nl_(nl),BARE_IA_(BARE_IA),BARE_AB_(BARE_AB),CIAQ_(CIAQ), CIJQ_(CIJQ),
// 				 MaxNCIA_(0), MaxNAB_(0), MaxNC_(0),NCIAReads_(0), NABReads_(0), NCReads_(0),
// 				 test_lag( ::cchem::rimp2_gradient::detail::RIAsync(NMB, CIAQ_, no_, BARE_IA_,no_,1,1) ),
// 				 test_lag2( ::cchem::rimp2_gradient::detail::RIAsync(NMB, BARE_AB_,0,nv_,CIJQ_,0, no_, 0,1,1) )
// 	//				 test_lag2( ::cchem::rimp2_gradient::detail::RIAsync(NMB, CIJQ_, no_, BARE_AB_,nv_,0,1) )
// 	{
// 		//		this->set_up_block_information();
// 		//		::cchem::rimp2_gradient::detail::RIAsync test_lag_1(NMB, CIAQ_, no_, BARE_IA_,no_, 0, 1);
// 	};

// 	~RIZVector(){}//std::cout << "deleting RIZVector" << std::endl;};

// 	void build_orbital_hessian(MapMatrixXd &u, MapMatrixXd &unxt)
// 	{

// 		utility::timer timer;

// 		timer.reset();
// 		test_lag.DoOperation2Async(  ::cchem::rimp2_gradient::detail::lag_vovo_functor(no_, nv_, u, unxt) );
// 		std::cout << "vovo_functor " << timer << std::endl;

// 		timer.reset();
// 		test_lag2.DoOperation2Async(	::cchem::rimp2_gradient::detail::lag_vvoo_functor(nv_, no_, u, unxt) );
// 		std::cout << "vvoo_functor " << timer << std::endl;
// 		//		test_lag2.DoOperation2Async(	::cchem::rimp2_gradient::detail::lag_vvoo_functor(no_, nv_, u, unxt) );

// 	}//build_orbtial_hessian_alt

// };//class RIZVector





 }//namespace detail
 }//namespace rimp2_gradient
 }//namespace cchem




#include <device.hpp>
#if HAVE_CUBLAS
#include <deviceGradientEngine.hpp>
#endif
namespace cchem{
namespace rimp2_gradient{
namespace detail{


/*
  @brief RI z-vector derived class
  @author LBR
  @detail this builds the orbital hessian with JK matrices
  @param N number of OBS functions
  @param nd number of doubly occupied MOs
  @param ns number of singly occpied MOs
  @param nv number of virtual MOs
  @param nl number of ABS functions
  @param maxitc maxium zvector iterations
  @param uconv convergence criteria
  @param small convergence criteria
  @param ea active orbital energies
  @param ev virtual orbtital energies
  @param NMB number of MB allocated from IO
  @param Ca occupied LCAO coefficients
  @param Cv virtual LCAO coefficients
  @param LM1_MN_SYM array pointer to (mu nu|L) integrals
  @param pOcc occupied LCAO coefficents (probably redundant)
  @param pe parallel environment  
  @param pGPU optional pointer to device class
*/

class JK_RIZVector : public ::cchem::rimp2_gradient::detail::ZVector,
		     public ::cchem::rimp2_gradient::detail::JK
{
    
    const size_t nl_, no_;
    const size_t N_;
    boost::reference_wrapper<const Matrix> &Ca_;
    boost::reference_wrapper<const Matrix> &Cv_;
    
    //(mu nu|P)Lm1
    Array<double> *LM1_MN_SYM_;

    //for half transformed JK matrices
    double *pHalf_;

    ::cchem::rimp2::device *pGPU_;

public:


    JK_RIZVector(size_t N,const size_t &nd, const size_t &ns, const size_t &nv, const size_t &nl,
		 const size_t &maxitc, double uconv, double small,
		 boost::reference_wrapper<const vector_range> ea,
		 boost::reference_wrapper<const vector_range> ev,
		 const size_t &NMB,
		 boost::reference_wrapper<const Matrix> Ca,
		 boost::reference_wrapper<const Matrix> Cv,
		 Array<double> *LM1_MN_SYM,double *pOcc,
		 Parallel &pe,
		 ::cchem::rimp2::device *pGPU = NULL
		 ):
	::cchem::rimp2_gradient::detail::ZVector(nd, ns, nv, maxitc,uconv,small,ea,ev,pe),
	::cchem::rimp2_gradient::detail::JK(nl,N),
	nl_(nl),no_(nd+ns),N_(N),
	Ca_(Ca),Cv_(Cv), LM1_MN_SYM_(LM1_MN_SYM),pGPU_(pGPU)
    {
	const int max_threads_ = omp_get_max_threads();

	//half transformed J/K matrix (scratch space)
	pHalf_ = new double[N_*no_];

	//always Occupied LCAO coefficients
	pLeft_.push_back(pOcc);

	//zvector update * virt_coeff_mat
	pRight_.push_back(new double[no_*N_]);
	
	//densty(ao,ao) update = pLeft_*pRight_
	pD_.push_back(new double[N_*N_]);
	
	//lower trangular density matrix
	pDmn_.push_back(new double[mnSize_]);

	//density (ao,ao) contracted with (ao ao|aux) integrals
	pDAux_.push_back(new double[nl_]);
	
    };
    
    ~JK_RIZVector(){
	delete [] pHalf_;
	for(int i = 0; i < pD_.size(); i++){
	    delete [] pD_[i];
	    delete [] pRight_[i];
	    delete [] pDmn_[i];
	    delete [] pDAux_[i];
	}//i
    }//~
    
    virtual void build_orbital_hessian(MapMatrixXd &u, MapMatrixXd &unxt);
    
};//class JK_RIZVector

    
}//namespace detail
}//namespace rimp2_gradient
}//namespace cchem


#endif /* SRC_RI_MP2_RI_ZVECTOR_HPP_ */
