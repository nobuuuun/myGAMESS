/*
 * ri-int.hpp
 *
 *  Created on: Aug 14, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_INTEGRALS_RI_INT_HPP_
#define LIBCCHEM_SRC_INTEGRALS_RI_INT_HPP_


#include "basis/basis.hpp"
#include "integrals/quartet.hpp"
#include "integrals/screening.hpp"
#include "integrals/ri-rysq.hpp"
#include "integrals/eri.hpp"
#include "boost/utility/profiler.hpp"

#include <vector>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "foreach.hpp"

namespace integrals {





template<typename T>
struct TwoCenterInt {

	typedef Basis::Shell Shell;

	typedef std::vector< boost::array<int,2> > Doublets;

	TwoCenterInt(const Basis &auxbasis)
	: centers_(auxbasis.centers())
	{
		max_ = auxbasis.data().size();
	}

//	template<typename Q, typename S>
//	void operator()(const Q &q, const S &s, const Doublets &doublets, double *data)
//	{
//		Doublet<const Shell::Data&> doublet =
//		{ detail::shell(q), detail::shell(s) }; // Q >= S !!! when sorted !!!
//		std::fill_n(data, doublet.size(), 0);
//		evaluate(doublet, centers_, doublets, data);
//	}

	template<typename Q, typename S>
	void operator()(const Q &q, const S &s, const Doublets &doublets, double *data);

private:

	typedef boost::array<Shell::Data::key_type,2> key_type;
	boost::ptr_map<key_type, integrals::rysq::TwoCenterInt<T> > cache_;
	std::vector<Basis::Center> centers_;
	size_t max_; //number of unique shells --> can have difference centers
	void evaluate(Doublet<const Shell::Data&> doublet,
			const std::vector<Basis::Center> &centers,
			const Doublets &doublets,
			double *G) {

		key_type key = {{ doublet.q.key(), doublet.s.key()}};

		if (cache_.size() > max_) {
			cache_.clear();
		}
		if (!cache_.count(key)) {
			cache_.insert(key, new integrals::rysq::TwoCenterInt<T>(doublet));
		}

		double cutoff = 0.0;
		integrals::rysq::TwoCenterInt<T> &rysq = cache_.at(key);
		rysq(centers, doublets, G, cutoff);
	}

}; //namespace TwoCenterInt



template<>
template<typename Q, typename S>
void TwoCenterInt< ::rysq::TwoCenterEri >::operator()(const Q &q, const S &s, const Doublets &doublets, double *data){

	Doublet<const Shell::Data&> doublet =
	{ detail::shell(q), detail::shell(s) }; // Q >= S !!! when sorted !!!
	std::fill_n(data, doublet.size(), 0);
	evaluate(doublet, centers_, doublets, data);

}

//yeah this sucks but the factor of "3" is the only difference
template<>
template<typename Q, typename S>
void TwoCenterInt< ::rysq::TwoCenterDerivativeEri >::operator()(const Q &q, const S &s, const Doublets &doublets, double *data){

	Doublet<const Shell::Data&> doublet =
	{ detail::shell(q), detail::shell(s) }; // Q >= S !!! when sorted !!!
	std::fill_n(data, doublet.size()*6, 0);
	evaluate(doublet, centers_, doublets, data);

}


//yeah this sucks but the factor of "3" is the only difference
template<>
template<typename Q, typename S>
void TwoCenterInt< ::rysq::TwoCenterDerivativeOverlap >::operator()(const Q &q, const S &s, const Doublets &doublets, double *data){

	Doublet<const Shell::Data&> doublet =
	{ detail::shell(q), detail::shell(s) }; // Q >= S !!! when sorted !!!
	std::fill_n(data, doublet.size()*6, 0);
	evaluate(doublet, centers_, doublets, data);

}






template<typename T>
	struct ThreeCenterInt {

	typedef Basis::Shell Shell;

	typedef std::vector< boost::array<int,3> > Triplets;

	ThreeCenterInt(const Basis &obsbasis, const Basis &auxbasis)
	: obscenters_(obsbasis.centers()),auxcenters_(auxbasis.centers())
	{
		obsmax_ = obsbasis.data().size(); //number of unique obs shells
		auxmax_ = auxbasis.data().size(); //number of unique aux shells
	}

//	template<typename L, typename Q, typename S>
//	void operator()(const Q &q, const S &s, const L &l, const Triplets &triplets, double *data)
//	{
//		Triplet<const Shell::Data&> triplet = { detail::shell(q), detail::shell(s), detail::shell(l) };
//		std::fill_n(data, triplet.size(), 0);
//		evaluate(triplet, obscenters_, auxcenters_, triplets, data);
//	}

		template<typename L, typename Q, typename S>
		void operator()(const Q &q, const S &s, const L &l, const Triplets &triplets, double *data);

private:

	typedef boost::array<Shell::Data::key_type,3> key_type;
	boost::ptr_map<key_type, integrals::rysq::ThreeCenterInt<T> > cache_;
	std::vector<Basis::Center> obscenters_;
	std::vector<Basis::Center> auxcenters_;
	size_t auxmax_; //number of unique shells --> can have difference centers
	size_t obsmax_;
	void evaluate(Triplet<const Shell::Data&> triplet,
			const std::vector<Basis::Center> &obscenters,
			const std::vector<Basis::Center> &auxcenters,
			const Triplets &triplets,
			double *G) {

		key_type key = {{ triplet.q.key(), triplet.s.key(), triplet.l.key()}};
		if (cache_.size() > auxmax_) { //maybe obs max is better here?
			cache_.clear();
		}
		if (!cache_.count(key)) {
			cache_.insert(key, new integrals::rysq::ThreeCenterInt<T>(triplet));
		}
		double cutoff = 0.0;
		integrals::rysq::ThreeCenterInt<T> &rysq = cache_.at(key);
		rysq(auxcenters, obscenters, triplets, G, cutoff);
	}

}; //namespace ThreeCenterInt

template<>
template<typename L, typename Q, typename S>
void ThreeCenterInt< ::rysq::ThreeCenterEri >::operator()(const Q &q, const S &s, const L &l, const Triplets &triplets, double *data)
{
	Triplet<const Shell::Data&> triplet = { detail::shell(q), detail::shell(s), detail::shell(l) };
	std::fill_n(data, triplet.size(), 0);
	evaluate(triplet, obscenters_, auxcenters_, triplets, data);
}

template<>
template<typename L, typename Q, typename S>
void ThreeCenterInt< ::rysq::ThreeCenterDerivativeEri >::operator()(const Q &q, const S &s, const L &l, const Triplets &triplets, double *data)
{
	Triplet<const Shell::Data&> triplet = { detail::shell(q), detail::shell(s), detail::shell(l) };
	std::fill_n(data, triplet.size()*9, 0);
	evaluate(triplet, obscenters_, auxcenters_, triplets, data);
}




} //namespace integrals


#endif /* LIBCCHEM_SRC_INTEGRALS_RI_INT_HPP_ */
