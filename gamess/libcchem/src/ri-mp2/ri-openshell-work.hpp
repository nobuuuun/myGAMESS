/*
 * ri-openshell-work.hpp
 *
 *  Created on: Jun 17, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_MP2_RI_OPENSHELL_WORK_HPP_
#define LIBCCHEM_SRC_MP2_RI_OPENSHELL_WORK_HPP_

#include <Eigen/Dense>

namespace cchem{
namespace rimp2{
namespace detail {

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;

struct ri_openshell_work : boost::noncopyable {

private:
    	double * pfock_;
	size_t nd_,ns_,nv_,nl_;
	Runtime::Memory *memory_;

public:

	ri_openshell_work(const size_t &nd, const size_t &ns,
			const size_t &nv, const size_t &nl,
			Runtime::Memory &memory):
		nd_(nd), ns_(ns), nv_(nv), nl_(nl),memory_(&memory)
	{
	    pfock_ = memory.malloc<double>(nv_*nd_);
	    std::fill(&pfock_[0],&pfock_[nv_*nd_], 0); //initialize vector
   	};//ri_openshell_work


	~ri_openshell_work(){ //release memory
	    memory_->free(pfock_);
	}//~ri_openshell_work

    double * get_pfock() {return pfock_;};


	template<class T>
	void ri_exchange_eri(Array<double> *Vp, double *ptr_bi, double *ptr_eri, T &e_m, double *data_pqm1){

		MapMatrixXd bi(ptr_bi, (nv_+ns_), nl_);   //nv*auxbasis.size()
		MapMatrixXd lm1(data_pqm1, nl_, nl_ );   //auxbasis.size()*auxbasis.size()
		for(int s = 0; s<ns_; s++){

			MapMatrixXd eri(ptr_eri,(s+1),(s+1));

			size_t start[] = { (s+nd_)*(nv_+ns_), 0 };
			size_t finish[] = { (s+nd_+1)*(nv_+ns_), nl_ };
			Vp->get(ptr_bi, start, finish);

			bi = bi*lm1.triangularView<Eigen::Upper>();
			eri = (bi.block(0,0,s+1,nl_) ) * ( (bi.block(0,0,s+1,nl_)).transpose()) ;// form eris

			for(int t = 0; t<=s ; t++){
				e_m[t] += eri(t,t)/2;
				if( t != s )e_m[s] += eri(t,t)/2;
			}//t

		}//s

	} //ri_exhange_eri




	template<class T>
	void ri_exchange_eri2(Array<double> *VpT, double *ptr_bi, double *ptr_eri, T &e_m){

		MapMatrixXd bi(ptr_bi, nl_, (nv_+ns_) );   //nv*auxbasis.size()

		for(int s = 0; s<ns_; s++){

			MapMatrixXd eri(ptr_eri,(s+1),(s+1));

			size_t start[] = { 0, (s+nd_)*(nv_+ns_) };
			size_t finish[] = { nl_, (s+nd_+1)*(nv_+ns_) };

			VpT->get(ptr_bi, start, finish);

			eri = (bi.block(0,0,nl_,s+1) ).transpose() * ( (bi.block(0,0,nl_,s+1)) ) ;// form eris

			for(int t = 0; t<=s ; t++){
				e_m[t] += eri(t,t)/2;
				if( t != s )e_m[s] += eri(t,t)/2;
			}//t

		}//s

	} //ri_exhange_eri2



};//ri_openshell_work

}//detail
}//rimp2
}//cchem


#endif /* LIBCCHEM_SRC_MP2_RI_OPENSHELL_WORK_HPP_ */
