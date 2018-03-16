/*
 * overlap-integrals.hpp
 *
 *  Created on: Sep 1, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_INTEGRALS_OVERLAP_INTEGRALS_HPP_
#define LIBCCHEM_SRC_INTEGRALS_OVERLAP_INTEGRALS_HPP_

#include <boost/noncopyable.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/thread.hpp>

#include "hermite.hpp"

namespace cchem{
namespace overlap{


struct TwoCenterDerivativeOverlap : boost::noncopyable {

	typedef std::vector< boost::array<int,2> > Doublets;



	TwoCenterDerivativeOverlap(){}//constructor

	~TwoCenterDerivativeOverlap(){};

	template<typename Q, typename S>
	double *operator()(const Q &q, const S &s, double *data){

//		typedef const basis::shell_data data_wrapper;

		const basis::shell_data a = (*q.data());
		const basis::shell_data b = (*s.data());
		std::cout << q.data()->K() << " " << q.data()->L() <<" ||||| " << s.data()->K() << " " << s.data()->L()<< std::endl;


		int ak = a.K();
		int bk = b.K();

		double xi = q.center().at(0);
		double yi = q.center().at(0);
		double zi = q.center().at(0);

		double xj = s.center().at(0);
		double yj = s.center().at(0);
		double zj = s.center().at(0);

		double rr = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);

		for (int Ki = 0; Ki < ak; ++Ki){
			double ai = exp(-a(Ki));
			double arri = rr*ai;
			double axi = ai*xi;
			double ayi = ai*yi;
			double azi = ai*zi;


			for (int Kj = 0; Kj < bk; ++Kj){
				double aj = exp(-b(Kj));

				double AB1 = 1/(ai+aj);
				double jrriab = aj*arri*AB1;
				double fac = exp(-jrriab);
				double ax = (axi+aj*xj)*AB1;
				double ay = (ayi+aj*yj)*AB1;
				double az = (azi+aj*zj)*AB1;


				double T = std::sqrt(AB1);
				double x0 = ax;
				double y0 = ay;
				double z0 = az;
				int Li1d = a.L() +2; // +1 for derivative
				int Lj1 = b.L()+1;

				for (int j = 0; j < Lj1; j++){
					int nj = j;

					for (int i = 0; i < Li1d; i++){
						int ni = i;

						double xint = 0.0;
						double yint = 0.0;
						double zint = 0.0;
						int npts = (ni+nj-2)/2+1;

					}//i
				}//j

			}//Kj

		} //Ki




	}//double *operator()

}; //struct TwoCenterDerivativeOverlap









template<typename T>
struct TwoCenterInt : boost::noncopyable {
private:
	Thread thread_;

	boost::reference_wrapper<const Basis> basis_;

	T twoc_work_;
public:

	typedef std::vector< boost::array<int,2> > Doublets;

	static
	boost::array<size_t,1> memory(const Basis & auxbasis);

	TwoCenterInt(boost::reference_wrapper<const Basis> basis)
	:	basis_(basis)
	{
		BOOST_AUTO(memory, this->memory(basis_.get()));
		for (int i = 0; i < memory.size(); ++i) {
			thread_.data[i] = thread_.malloc(memory[i]); //maxb*maxb*maxs*maxs;
		}
	}//constructor

	~TwoCenterInt() {thread_.free();}


	double *operator()(const Basis::Shell &Q, const Basis::Shell &S) {

		BOOST_AUTO(const &basis, basis_.get());

		double *G = thread_.data[0];    //maxs*maxs

//		boost::array<int,2> doublet = {{S.index() , Q.index()}};
		twoc_work_(S, Q, G);
		return G;

	} //*operator()




};//struct OverlapInt


////explicit (full) template specialization
//template<>
//	inline boost::array<size_t,1> OverlapInt< TwoCenterOverlap >::memory(const Basis & auxbasis){
//	size_t maxs = auxbasis.max().size(); // #BF in largest shell (ie f-shell: maxs=10)
//	boost::array<size_t,1> memory;
//	memory[0] = maxs*maxs;
//	std::cout << "Auxiliary_TwoC_Work<integrals::TwoC_Eri>::memory" << std::endl;
//	return memory;
//}

//explicit (full) template specialization
template<>
inline boost::array<size_t,1> TwoCenterInt< TwoCenterDerivativeOverlap >::memory(const Basis & basis){
	size_t maxs = basis.max().size(); // #BF in largest shell (ie f-shell: maxs=10)
	boost::array<size_t,1> memory;
	memory[0] = maxs*maxs*6;
	std::cout << "TwoCenterInt< TwoCenterDerivativeOverlap >::memory" << std::endl;
	return memory;
}




}//overlap
}//cchem


#endif /* LIBCCHEM_SRC_INTEGRALS_OVERLAP_INTEGRALS_HPP_ */
