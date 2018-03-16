/*
 * twoc-rysq-eri.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_RI_RYSQ_INT_HPP_
#define LIBCCHEM_RYSQ_SRC_RI_RYSQ_INT_HPP_

#include "rysq-core.hpp"  //need this for shell struct
#include "ri-rysq-core.hpp"  //has doublet and triplet interfaces

#include <vector>
#include <memory>


namespace rysq {



    //general two-center integral type
    class TwoCenterInt {
    public:

	struct Parameters{
	    double cutoff, scale;
	    Parameters(double cutoff = 1.0e-10, double scale = 1)
		: cutoff(cutoff), scale(scale) {}
	};//

	TwoCenterInt(){};

	~TwoCenterInt(); //destroy kernel

	void operator()(const Doublet<Center> &centers, double *I,   //business end -> call to kernel
			const Parameters &parameters = Parameters());

	void operator()(const std::vector<Center> &centers,  //intermediate
			const Int2 &doublet, double *I,
			const Parameters &parameters = Parameters()) {

	    Doublet<Center> r(centers.at(doublet[0]), centers.at(doublet[1]));
	    (*this)(r, I, parameters);  //calls business end
	}

	//start here
	void operator()(const std::vector<Center> &centers,   //declare looper : calls intermediate
			const std::vector<Int2> &doublet,
			double *Eri,
			const Parameters &parameters = Parameters());

    protected:
	void *kernel_;
	size_t size_;

    }; //class TwoCenterInt

    //two-center eri type
    class TwoCenterEri : public TwoCenterInt {
    public:
	TwoCenterEri(const Doublet<Shell> &shells);
    };


    //two-center derivative eri type
    class TwoCenterDerivativeEri : public TwoCenterInt {
    public:
	TwoCenterDerivativeEri(const Doublet<Shell> &shells);
    };

    //two-center derivative overlap type
    class TwoCenterDerivativeOverlap : public TwoCenterInt {
    public:
	TwoCenterDerivativeOverlap(const Doublet<Shell> &shells);
    };




























    class ThreeCenterInt {
    public:
	class Impl;
	size_t size;

	struct Parameters {
	    double cutoff, scale;
	    Parameters(double cutoff = 1.0e-10, double scale = 1)
		: cutoff(cutoff), scale(scale) {}
	};

	ThreeCenterInt(){}; //basically sets up kernel

	~ThreeCenterInt(); //destroy kernel

	void operator()(const Triplet<Center> &centers, //business end -> call to kernel
			double *I,
			const Parameters &parameters = Parameters());

	void operator()(const std::vector<Center> &obscenters,  //intermediate
			const std::vector<Center> &auxcenters,
			const Int3 &triplet, double *I,
			const Parameters &parameters = Parameters()) {

	    Triplet<Center> r(obscenters.at(triplet[0]), obscenters.at(triplet[1]), auxcenters.at(triplet[2]));
	    (*this)(r, I, parameters);  //calls business end
	}

	void operator()(const std::vector<Center> &obscenters,   //declare looper : calls intermediate
			const std::vector<Center> &auxcenters,
			const std::vector<Int3> &triplet,
			double *Eri,
			const Parameters &parameters = Parameters());


    protected:
	void *kernel_;
	size_t size_;

    }; //class ThreeCenterInt


    //three-center eri type
    class ThreeCenterEri : public ThreeCenterInt {
    public:
	ThreeCenterEri(const Triplet<Shell> &shells);
    };//ThreeCenterEri

    //three-center derivative eri type
    class ThreeCenterDerivativeEri : public ThreeCenterInt {
    public:
	ThreeCenterDerivativeEri(const Triplet<Shell> &shells);
    };//ThreeCenterDerivativeEri










}//namespace rysq



#endif /* LIBCCHEM_RYSQ_SRC_RI_RYSQ_INT_HPP_ */
