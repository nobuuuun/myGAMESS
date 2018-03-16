/*
 * twoc-eri-transform.hpp
 *
 *  Created on: Apr 2, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_RI_INT_TRANSFORM_HPP_
#define LIBCCHEM_RYSQ_SRC_RI_INT_TRANSFORM_HPP_

namespace rysq {

    namespace twoc_transform {


	template<class bra = void> //primary struct template
	struct TwoCenterTransform; //forward declaration, empty base class


	template<>
	struct TwoCenterTransform<> //explicit specialization (no template parameters)
	{
	    struct Data : kernel::TwoCenterTransform<>::Data {
		Data() {}
		Data(double *Q) : data_(Q) {}
		double* begin() { return data_; }

	    private:
		double *data_;

	    };
	};


	template<class bra> //primary struct template
	struct TwoCenterEriTransform : kernel::TwoCenterTransform<bra> {
	    typedef kernel::TwoCenterTransform<bra> Base;
	    typedef TwoCenterTransform<>::Data Data;
	    Base& operator()(typename Base::Data &data) {
		this->data_ = static_cast<Data&>(data);
		return *this;
	    }//Base& operator()
	    void operator()(int kl, const double *Q, double scale) {
		const size_t ni = bra::A::size;
		const size_t size = ni; //*nj;
		double *data = data_.begin() + kl*size;
		for (size_t i = 0; i < size; ++i) {
		    data[i] += scale*Q[i];
		}
	    }// void operator()
	private:
	    Data data_;
	};









	template<class bra,class ket> //primary struct template (derived type)
	struct TwoCenterDerivativeEriTransform : kernel::TwoCenterTransform<bra> {
	    typedef kernel::TwoCenterTransform<bra> Base;
	    typedef TwoCenterTransform<>::Data Data;
	    Base& operator()(typename Base::Data &data) {
		this->data_ = static_cast<Data&>(data);
		return *this;
	    }//Base& operator()

	    void operator()(int kl, const double *Q, double scale) {
		const size_t ni = bra::A::size;
		const size_t nj = ket::A::size;
		const size_t size = ni;
		double *data = data_.begin();

		size_t DStart = kl*size;
		size_t QStart = 0;

		for (size_t nderiv = 0; nderiv<3; nderiv++){
		    for (size_t i = 0; i < ni; i++){
			data[DStart+i] += scale*Q[QStart+i];
		    } //i
		    DStart +=nj*size;
		    QStart +=size;
		} //nderiv

		size_t AStart = kl*size;
		size_t BStart = 3*ni*nj + kl*size;
		for(size_t i=0; i<3; i++){
		    for(size_t b = 0; b< ni; b++){
			data[BStart +b] -= data[AStart+b];
		    }//b
		    AStart += ni*nj;
		    BStart += ni*nj;
		}//i




	    }// void operator()
	private:
	    Data data_;
	};



    } //namespace twoc_transform
} //namespace rysq











namespace rysq {
    namespace threec_transform {





	template<class bra = void>
	struct ThreeCenterTransform;

	template<>
	struct ThreeCenterTransform<> {
	    struct Data : kernel::ThreeCenterTransform<>::Data {
		Data() {}
		Data(double *Q) : data_(Q) {}
		double* begin() { return data_; }
	    private:
		double *data_;
	    };
	};


	template<class bra>
	struct ThreeCenterEriTransform : kernel::ThreeCenterTransform<bra> {
	    typedef kernel::ThreeCenterTransform<bra> Base;
	    typedef ThreeCenterTransform<>::Data Data;

	    Base& operator()(typename Base::Data &data) {
		this->data_ = static_cast<Data&>(data);
		return *this;
	    }

	    void operator()() {};
	    void operator()(int k, const double *Q, double scale) {
		const size_t ni = bra::A::size;
		const size_t nj = bra::B::size;
		const size_t size = ni*nj;
		double *data = data_.begin() + k*size;
		for (size_t i = 0; i < size; ++i) {
		    data[i] += scale*Q[i];
		}
	    }

	private:
	    Data data_;

	};


	template<class bra,class ket>
	struct ThreeCenterDerivativeEriTransform : kernel::ThreeCenterTransform<bra> {
	    typedef kernel::ThreeCenterTransform<bra> Base;
	    typedef ThreeCenterTransform<>::Data Data;

	    Base& operator()(typename Base::Data &data) {
		this->data_ = static_cast<Data&>(data);
		return *this;
	    }

	    void operator()(int k, const double *Q, double scale) {
		const size_t ni = bra::A::size;
		const size_t nj = bra::B::size;
		const size_t nk = ket::A::size;
		const size_t size = ni*nj;

		double *data = data_.begin();

		size_t DStart = k*size;
		size_t QStart = 0;
		for (size_t nderiv = 0; nderiv<6; nderiv++){

		    for (size_t i = 0; i < ni*nj; i++){
			data[DStart+i] += scale*Q[QStart+i];
		    }
		    DStart +=nk*size;
		    QStart +=size;
		}

		// size_t AStart = k*size;
		// size_t BStart = 3*ni*nj*nk + k*size;
		// size_t CStart = 6*ni*nj*nk + k*size;
		// for (size_t i = 0; i < 3; ++i) {
		//     for( size_t c = 0; c < ni*nj; c++){
		// 	data[CStart +c] -= data[AStart+c] + data[BStart+c];
		//     }//c
		//     AStart += ni*nj*nk;
		//     BStart += ni*nj*nk;
		//     CStart += ni*nj*nk;
		// }//i

	    }


	    void operator()(){
		const size_t ni = bra::A::size;
		const size_t nj = bra::B::size;
		const size_t nk = ket::A::size;

		double *data = data_.begin();

		size_t BStart = 3*ni*nj*nk; // + k*size;
		size_t CStart = 6*ni*nj*nk; // + k*size;
		for( size_t c = 0; c < 3*ni*nj*nk; c++){
		    data[CStart +c] = -(data[c] + data[BStart+c]);
		}//c

		//cblas_dcopy(ni*nj*nk,&data[c],1,&data[CStart +c]);
		//then to daxpy ???

	    }//operator()


	private:
	    Data data_;

	};



    } //namespace threec_transform
} //namespace rysq



#endif /* LIBCCHEM_RYSQ_SRC_RI_INT_TRANSFORM_HPP_ */
