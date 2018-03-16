#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include "basis/basis.hpp"
#include "basis/normalize.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace detail {
    struct wavefunction {
	typedef Basis basis_type;
    };
}
/**
*   @brief Wavefunction class - self-documenting
*
*   @author Andrey Asadchev
*
*   @date 12-23-14 Luke
*   	-Added set_single_occupied(), set_double_occupied(),
*   	set_single_and_virtuals(),single(), double_occ(),
*   	single_virtual(), single_occupied_, double_occupied_,
*   	and single_virtual_
*/
struct Wavefunction {
    typedef detail::wavefunction::basis_type Basis;
    typedef boost::numeric::ublas::matrix<
	double, boost::numeric::ublas::column_major> matrix_type;
    typedef boost::numeric::ublas::matrix_range<const matrix_type> matrix_range;
    typedef boost::numeric::ublas::vector<double> vector_type;
    typedef boost::numeric::ublas::vector_range<const vector_type> vector_range;

    struct Orbitals {
	Orbitals() : size_(0), start_(0), stop_(0) {}
	Orbitals(size_t start, size_t stop)
	    : size_(stop - start), start_(start), stop_(stop) {}
	size_t size() const { return size_; }
	size_t start()const {return start_; }
	size_t stop() const { return stop_; }
	operator boost::numeric::ublas::range() const {
	    return boost::numeric::ublas::range(start_, stop_);
	}
    private:
	size_t size_, start_, stop_;
    };

    Wavefunction(const Basis& basis) : basis_(basis) {}

    Wavefunction(const Basis& basis, const Basis& auxbasis) : basis_(basis), auxbasis_(auxbasis), do_spherical_(0) {}

    void set_occupied(size_t start, size_t stop, size_t core) {
	occupied_ = Orbitals(start, stop);
	active_ = Orbitals(start+core, stop);
    }

    void set_virtuals(size_t start, size_t stop) {
	virtuals_ = Orbitals(start, stop);
    }

    void set_single_occupied(size_t start, size_t stop) {
	single_occupied_ = Orbitals(start, stop);
    }

    void set_double_occupied(size_t start, size_t stop) {
	double_occupied_ = Orbitals(start, stop);
    }

    void set_single_and_virtuals(size_t start, size_t stop) {
	single_virtual_ = Orbitals(start, stop);
    }

    template<class E>
    void set_C(const boost::numeric::ublas::matrix_expression<E> &C) {
	C_ = C;
    }

    template<class E>
    void set_sperical_coeff(const boost::numeric::ublas::matrix_expression<E> &SC) {
        SphericalCoeff_.push_back(SC);
    }


    void set_spherical() {
        do_spherical_ = (int)1;
    }

    template<class E>
    void set_e(const boost::numeric::ublas::vector_expression<E> &e) {
	e_ = e;
    }

    template<class E>
    void set_F(const boost::numeric::ublas::matrix_expression<E> &F) {
	F_ = F;
    }

    const Basis& basis() const { return basis_; }
    const Basis& auxbasis() const { return auxbasis_; }
    const int& do_spherical() const { return do_spherical_; }

    const Orbitals& occupied() const { return occupied_; }
    const Orbitals& active() const { return active_; }

    const Orbitals& single() const { return single_occupied_; }

    const Orbitals& double_occ() const { return double_occupied_; }

    const Orbitals& single_virtual() const { return single_virtual_; }

    const Orbitals& virtuals() const { return virtuals_; }
    Orbitals orbitals() const {
	int start = std::min(occupied_.start(), virtuals_.start());
	int finish = std::max(occupied_.stop(), virtuals_.stop());
	return Orbitals(start, finish);
    }

    const vector_type& e() const { return e_; }
    vector_range e(const Orbitals& r) const {
	return boost::numeric::ublas::project(e_, r);
    }

    const matrix_type& F() const { return F_; }
    const matrix_type& C() const { return C_; }

    const matrix_type& spherical_coeff(int i) const {return SphericalCoeff_.at(i);}

    double * spherical_coeff_ptr(int i) {return SphericalCoeff_.at(i).data().begin();}

    matrix_range C(const Orbitals& r) const {
	namespace ublas = boost::numeric::ublas;
	return ublas::project(C(), ublas::range(0, C().size1()), r);
    }

    void normalize() {
	basis::normalize<0>(basis_.N(), C_);
    }

    void sort() {
	basis_.sort();
	reorder();
    }

    void auxsort() {
	auxbasis_.sort();
	//auxbasis_.reverse();
    }

    void add_spherical_info() {
	auxbasis_.add_spherical_info();
    }

    void reverse() {
	basis_.reverse();
	reorder();
    }

private:

    Basis basis_;
    Basis auxbasis_;
    int do_spherical_;
    Orbitals occupied_, active_, virtuals_,single_occupied_,double_occupied_,single_virtual_;
    matrix_type C_, F_, AuxC_;
    std::vector<matrix_type> SphericalCoeff_;
    vector_type e_;

    void reorder() {
	vector_type c(basis_.size());
	for (size_t j = 0; j < C_.size2(); ++j) {
	    for (size_t i = 0; i < basis_.size(); ++i) {
		c[i] = C_(basis_.index()[i], j);
	    }
	    column(C_,j).assign(c);
	}
    }

};

#endif // WAVEFUNCTION_HPP
