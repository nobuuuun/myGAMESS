/*
 * ri-driver.hpp
 *
 *  Created on: Jun 16, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_MP2_RI_DRIVER_HPP_
#define LIBCCHEM_SRC_MP2_RI_DRIVER_HPP_

#include "thread.hpp"
#include "core/wavefunction.hpp"
#include "integrals/eri.hpp"
//#include "integrals/twoc-eri.hpp"
#include "integrals/ri-int.hpp"
//#include "ri-integral-interface.hpp"

#include "blas.hpp"  //this is need for boost::numeric::ublas::matrix

#include <boost/noncopyable.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/thread.hpp>

//#include </opt/local/include/eigen3/Eigen/Core>

#include <Eigen/Dense>
namespace cchem {
namespace ri {

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> AlignedMapMatrixXd;
typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;

using cchem::Thread;

}
}

namespace cchem{
namespace ri{

template<typename T>
struct AuxiliaryTwoCenterInt : boost::noncopyable {
private:
	Thread thread_;

	boost::reference_wrapper<const Basis> auxbasis_; //if you change the order (auxbasis_ <--> twoc_eri_) thngs will mess up
	integrals::TwoCenterInt<T> twoc_work_;           //  --the members are initialized in the order they are declared, not the way they are listed in the constructor initializition list

public:

	typename integrals::TwoCenterInt<T>::Doublets doublets;

	static
	boost::array<size_t,1> memory(const Basis & auxbasis);

	AuxiliaryTwoCenterInt(boost::reference_wrapper<const Basis> auxbasis)
	:	auxbasis_(auxbasis), twoc_work_(auxbasis_.get())
	{
		BOOST_AUTO(memory, this->memory(auxbasis_.get()));
		for (int i = 0; i < memory.size(); ++i) {
			thread_.data[i] = thread_.malloc(memory[i]);
		}
	}//constructor

	~AuxiliaryTwoCenterInt() {thread_.free();}

	double *operator()(const Basis::Shell &Q, const Basis::Shell &S) {

		BOOST_AUTO(const &auxbasis, auxbasis_.get());

		double *G = thread_.data[0];    //maxs*maxs

		boost::array<int,2> doublet = {{S.index() , Q.index()}};
		doublets.push_back(doublet);
		twoc_work_(S, Q, doublets, G);
		doublets.pop_back(); //this is really ghetto, but i do this so doublets does not get crazy big

		return G;

	} //*operator()

};//AuxiliaryTwoCenterInt

//explicit (full) template specialization
template<>
	inline boost::array<size_t,1> AuxiliaryTwoCenterInt< ::rysq::TwoCenterEri >::memory(const Basis & auxbasis){
	size_t maxs = auxbasis.max().size(); // #BF in largest shell (ie f-shell: maxs=10)
	boost::array<size_t,1> memory;
	memory[0] = maxs*maxs;
//	std::cout << "Auxiliary_TwoC_Work<integrals::TwoC_Eri>::memory" << std::endl;
	return memory;
}

//explicit (full) template specialization
template<>
inline boost::array<size_t,1> AuxiliaryTwoCenterInt< ::rysq::TwoCenterDerivativeEri >::memory(const Basis & auxbasis){
	size_t maxs = auxbasis.max().size(); // #BF in largest shell (ie f-shell: maxs=10)
	boost::array<size_t,1> memory;
	memory[0] = maxs*maxs*6;
//	std::cout << "Auxiliary_TwoC_Work<integrals::TwoC_DEri>::memory" << std::endl;
	return memory;
}

//explicit (full) template specialization
template<>
inline boost::array<size_t,1> AuxiliaryTwoCenterInt< ::rysq::TwoCenterDerivativeOverlap >::memory(const Basis & basis){
	size_t maxs = basis.max().size(); // #BF in largest shell (ie f-shell: maxs=10)
	boost::array<size_t,1> memory;
	memory[0] = maxs*maxs*6;
//	std::cout << "Auxiliary_TwoC_Work< ::rysq::TwoCenterDerivativeOverlap >::memory" << std::endl;
	return memory;
}

























template <typename T>
struct AuxiliaryThreeCenterInt : boost::noncopyable{

	typedef Eigen::Map<const Eigen::MatrixXd> ConstMapMatrixXd;


private:
	Thread thread_;
	boost::reference_wrapper<const Basis> basis_;
	boost::reference_wrapper<const Basis> auxbasis_;
	boost::reference_wrapper<const Matrix> C_;
	boost::reference_wrapper<const Matrix> Cv_;
	integrals::ThreeCenterInt<T> threec_eri_;                   //  --the members are initialized in the order they are declared, not the way they are listed in the constructor initialization list


public:

	typename integrals::ThreeCenterInt<T>::Triplets triplets;
	static
	boost::array<size_t,4> memory(const Basis &auxbasis,
			const Basis &basis,
			const Matrix &C,
			const Matrix &Cv);


	AuxiliaryThreeCenterInt(boost::reference_wrapper<const Basis> auxbasis,
			boost::reference_wrapper<const Basis> basis,
			boost::reference_wrapper<const Matrix> C,
			boost::reference_wrapper<const Matrix> Cv)
	:	basis_(basis),auxbasis_(auxbasis),C_(C), Cv_(Cv),
		threec_eri_(auxbasis_.get(),basis_.get())
	{
		BOOST_AUTO(memory, this->memory(auxbasis_.get(), basis_.get(),C_.get(),Cv_.get()));
		for (int i = 0; i < memory.size(); ++i) {
			thread_.data[i] = thread_.malloc(memory[i]);
		}

	}

	~AuxiliaryThreeCenterInt(){thread_.free();}

	double *get_pointer(int i){return thread_.data[i];}


	/** (sq|l) -> (ai|l) */
	void contract(MapMatrixXd threec_eri_mat, double *buffer,
			ConstMapMatrixXd occ_coeff_mat, ConstMapMatrixXd virt_coeff_mat,
			MapMatrixXd t1){

		BOOST_AUTO(const &C, C_.get());
		BOOST_AUTO(const &Cv, Cv_.get());

		t1 = occ_coeff_mat*threec_eri_mat;

		MapMatrixXd t2(buffer, Cv.size1(), C.size1() ); //t2(nvir,nocc)

		t2 = virt_coeff_mat * t1.transpose();  //t2(nvir,nocc)

	}//contract

	/** (sq|l) -> (ab|l) */
	void contract_ab(MapMatrixXd threec_eri_mat,
			double *buffer,
			ConstMapMatrixXd virt_coeff_mat,
			MapMatrixXd t1_ab){

		BOOST_AUTO(const &Cv, Cv_.get());

		t1_ab = virt_coeff_mat * threec_eri_mat;

		MapMatrixXd t2(buffer, Cv.size1(), Cv.size1() ); //t2(nv,nv)

		t2 = virt_coeff_mat * t1_ab.transpose();

	}//contract_ab

	/** (sq|l) -> (ij|l) */
	void contract_ij(MapMatrixXd threec_eri_mat,
			double *buffer,
			ConstMapMatrixXd occ_coeff_mat,
			MapMatrixXd t1_ij){

		BOOST_AUTO(const &C, C_.get());

		t1_ij = occ_coeff_mat * threec_eri_mat;

		MapMatrixXd t2(buffer, C.size1(), C.size1() ); //t2(no,no)

		t2 = occ_coeff_mat * t1_ij.transpose();

	}//contract_ij



	double* operator()(const Basis::Shell Q, const Basis::Shell S, const Basis::Shell L) {

		double *G = thread_.data[0]; //max_obs_shell*max_obs_shell*max_aux_shell

		boost::array<int,3> triplet = {{S.index(), Q.index(), L.index()}};
		triplets.push_back(triplet);
		threec_eri_(S,Q,L,triplets, G );
		triplets.pop_back();

		return G;
	};

};//Auxiliary_ThreeC_ERI

template<>
inline boost::array<size_t,4> AuxiliaryThreeCenterInt< ::rysq::ThreeCenterEri >::memory(const Basis &auxbasis,
		const Basis &basis,
		const Matrix &C,
		const Matrix &Cv) {

	size_t max_obs_shell = basis.max().size();    // #functions in largest obs shell
	size_t max_aux_shell = auxbasis.max().size(); // #functions in largest aux shell
	size_t no = C.size1();   //number of occupied orbitals
	size_t nv = Cv.size1();  //number of virtual orbitals
	size_t N = basis.size(); //size of obs basis

	boost::array<size_t,4> memory;
	memory[0] = max_obs_shell*max_obs_shell*max_aux_shell; //G (mu nu|P) triplet
	memory[1] = N*nv;                                      //t1 & t1_ab 1/3 transformed buffer
	memory[2] = N*N*max_aux_shell;                         //ptr_threec (mu nu |P) integrals
	memory[3] = nv*nv*max_aux_shell;                       //ptr_ab (for ij, ia, and ab) 2/3 transformed ints (rs|P)

	return memory;
}



template<>
inline boost::array<size_t,4> AuxiliaryThreeCenterInt< ::rysq::ThreeCenterDerivativeEri >::memory(const Basis &auxbasis,
		const Basis &basis,
		const Matrix &C,
		const Matrix &Cv) {

	size_t max_obs_shell = basis.max().size();    // #functions in largest obs shell
	size_t max_aux_shell = auxbasis.max().size(); // #functions in largest aux shell
	size_t no = C.size1();   //number of occupied orbitals
	size_t nv = Cv.size1();  //number of virtual orbitals
	size_t N = basis.size(); //size of obs basis

	boost::array<size_t,4> memory;
	memory[0] = 9*max_obs_shell*max_obs_shell*max_aux_shell; //G (mu nu|P) triplet
	memory[1] = 0;                                           //not used
	memory[2] = 0;                                           //not used
	memory[3] = 0;                                           //not used



//std::cout << "AuxiliaryThreeCenterInt< ::rysq::ThreeCenterDerivativeEri >::memory" << std::endl;
return memory;
}




} // namespace ri
} //namespace cchema



#endif /* LIBCCHEM_SRC_MP2_RI_DRIVER_HPP_ */
