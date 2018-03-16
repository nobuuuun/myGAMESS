/*
 * twoc-rysq-core.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_TWOC_RYSQ_CORE_HPP_
#define LIBCCHEM_RYSQ_SRC_TWOC_RYSQ_CORE_HPP_

namespace rysq {

    typedef boost::array<int,2> Int2;

    namespace detail {

	//adapted for Doublets
	template<class F, typename Doublet>
	typename F::value_type twoc_add(const Doublet &q) {
	    return (F::value(q[0])+F::value(q[1]));
	}
	//adapted for Doublets
	template<class F, typename Doublet>
	typename F::value_type twoc_multiply(const Doublet &q) {
	    return (F::value(q[0])*F::value(q[1]));
	}
	//adapted for Doublets
	struct TwoC_L {
	    typedef size_t value_type;
	    static size_t value(const rysq::type &type) { return shell(type).L; }
	    static size_t value(const rysq::shell &shell) { return shell.L; }
	    template<typename Doublet>
	    static value_type twoc_add(const Doublet &q) {
		return detail::twoc_add<L>(q);
	    }
	};
	//adapted for Doublets
	struct TwoC_size {
	    typedef size_t value_type;
	    static size_t value(const rysq::type &type) { return shell(type).size; }
	    static size_t value(const rysq::shell &shell) { return shell.size; }
	    template<typename Doublet>
	    static value_type twoc_multiply(const Doublet &q) {
		return detail::twoc_multiply<size>(q);
	    }
	};
	//adapted for Doublets
	struct TwoC_K {
	    typedef size_t value_type;
	    static size_t value(const rysq::shell &shell) { return shell.K; }
	    template<typename Doublet>
	    static value_type twoc_multiply(const Doublet &q) {
		return detail::twoc_multiply<K>(q);
	    }
	};

    }//detail


    struct TwoC_State {
	TwoC_State(const Shell &A) : A(A) {}
	size_t K() const { return A.K; }
	size_t L() const { return A.L; }
	size_t size() const {return A.size; }
	size_t hybrid() const { return  A.is_hybrid(); }
    private:
	const Shell &A;
    };


    template<class  C>
    struct doublet_base : boost::array<C,2> {
	typedef boost::array<C,2> base;
	typedef C c_type[2];
	doublet_base() {}
	template<typename T>
	doublet_base(const T &t0, const T &t1)
	    : base(make(t0, t1)) {}
	template<typename T>
	doublet_base(const boost::array<T,2> &q)
	    : base(make(q[0], q[1])) {}
	operator const c_type&() const { return this->elems; }
    protected:
	template<typename T>
	static base make(const T &t0, const T &t1) {
	    boost::array<C,2> array = {{ t0, t1}};
	    return array;
	}
	size_t size() const;
    };

    template<class C>
    struct Doublet : doublet_base<C> {
	typedef doublet_base<C> base;
	Doublet() {}
	template<typename T>
	Doublet(const T &t0, const T &t1)
	    : base(t0, t1) {}
	template<typename T>
	Doublet(const T (&q)[2]) : base(q[0], q[1]) {}
	template<typename T>
  	Doublet(const boost::array<T,2> &q) : base(q[0], q[1]) {}
    };

    template<class C>
    inline std::ostream& operator<<(std::ostream &os, const Doublet<C> &doublet) {
	return os << "{"
		  << doublet[0] << "," << doublet[1]
		  << "}";
    }



    template<>
    struct Doublet<rysq::type> : doublet_base<rysq::type> {
	typedef doublet_base<rysq::type> base;
	Doublet(const type &a, const type &b)
	    : base(a,b) {}
	template<typename T>
  	Doublet(const boost::array<T,2> &doublet) : base(doublet) {}
	size_t L() const { return detail::L::add(*this); }
	size_t size() const { return detail::size::multiply(*this); }
    };



    template<>
    struct Doublet<rysq::shell> : doublet_base<rysq::shell> {
	typedef doublet_base<rysq::shell> base;
	Doublet(const shell &a, const shell &b)
	    : base(a,b) {}
	template<typename T>
  	Doublet(const boost::array<T,2> &doublet) : base(doublet) {}
	size_t L() const { return detail::L::add(*this); }
	size_t size() const { return detail::size::multiply(*this); }
	size_t K() const { return detail::K::multiply(*this); }
    };



    template<>
    class Doublet<Shell> {
    public:
	typedef TwoC_State Bra, Ket;

	static int size(int a, int b) {
	    return (shell(a).size*shell(b).size);
	}

	static size_t size(const rysq::type (&doublet)[2]) {
	    return size(doublet[0], doublet[1]);
	}

	Doublet(const Shell &a, const Shell &b)
	    : a_(a), b_(b)
	{
	    initialize();
	}
	Doublet(const Doublet &doublet)
	    : a_(doublet[0]), b_(doublet[1])
	{
	    initialize();
	}
	~Doublet() {}
	int L() const { return L_; }
	int size() const { return size_; }
	int hybrid() const { return hybrid_; }
	int nc() const { return nc_; }
	int K() const { return K_; }
	const Shell& operator[](size_t i) const { return *shells_[i]; }
	const Bra bra() const { return Bra((*this)[0]); }
	const Ket ket() const { return Ket((*this)[1]); }
	operator Doublet<rysq::type>() const { return cast<rysq::type>(); }
	operator Doublet<rysq::shell>() const  { return cast<rysq::shell>(); }
	boost::array<Shell*,2> data() const { return shells_; }
    private:
	boost::array<Shell*,2> shells_;
	Shell a_, b_;
	int L_, size_, hybrid_, nc_, K_;
	void initialize();
	void operator=(const Doublet&);
	template<typename T>
	Doublet<T> cast() const {
	    const Doublet &q = *this;
	    return Doublet<T>(T(q[0]), T(q[1]));
	}
    };


















    typedef boost::array<int,3> Int3;

    namespace detail {

	//adapted for Triplets
	template<class F, typename Triplet>
	typename F::value_type threec_add(const Triplet &q) {
	    return (F::value(q[0])+F::value(q[1])+F::value(q[2]));
	}
	//adapted for Triplets
	template<class F, typename Triplet>
	typename F::value_type threec_multiply(const Triplet &q) {
	    return (F::value(q[0])*F::value(q[1])*F::value(q[2]));
	}
	//adapted for Triplets
	struct ThreeC_L {
	    typedef size_t value_type;
	    static size_t value(const rysq::type &type) { return shell(type).L; }
	    static size_t value(const rysq::shell &shell) { return shell.L; }
	    template<typename Triplet>
	    static value_type threec_add(const Triplet &q) {
		return detail::threec_add<L>(q);
	    }
	};
	//adapted for Triplets
	struct ThreeC_size {
	    typedef size_t value_type;
	    static size_t value(const rysq::type &type) { return shell(type).size; }
	    static size_t value(const rysq::shell &shell) { return shell.size; }
	    template<typename Triplet>
	    static value_type threec_multiply(const Triplet &q) {
		return detail::threec_multiply<size>(q);
	    }
	};
	//adapted for Triplets
	struct ThreeC_K {
	    typedef size_t value_type;
	    static size_t value(const rysq::shell &shell) { return shell.K; }
	    template<typename Triplet>
	    static value_type threec_multiply(const Triplet &q) {
		return detail::threec_multiply<K>(q);
	    }
	};

    }//detail


    struct ThreeC_State {
	ThreeC_State(const Shell &A) : A(A) {}
	size_t K() const { return A.K; }
	size_t L() const { return A.L; }
	size_t size() const {return A.size; }
	size_t hybrid() const { return  A.is_hybrid(); }
    private:
	const Shell &A;
    };


    template<class  C>
    struct triplet_base : boost::array<C,3> {
	typedef boost::array<C,3> base;
	typedef C c_type[3];
	triplet_base() {}
	template<typename T>
	triplet_base(const T &t0, const T &t1, const T &t2)
	    : base(make(t0, t1, t2)) {}
	template<typename T>
	triplet_base(const boost::array<T,3> &q)
	    : base(make(q[0], q[1], q[2])) {}
	operator const c_type&() const { return this->elems; }
    protected:
	template<typename T>
	static base make(const T &t0, const T &t1, const T &t2) {
	    boost::array<C,3> array = {{ t0, t1, t2 }};
	    return array;
	}
	size_t size() const;
    };

    template<class C>
    struct Triplet : triplet_base<C> {
	typedef triplet_base<C> base;
	Triplet() {}
	template<typename T>
	Triplet(const T &t0, const T &t1, const T &t2)
	    : base(t0, t1, t2) {}
	template<typename T>
	Triplet(const T (&q)[3]) : base(q[0], q[1], q[2]) {}
	template<typename T>
 	Triplet(const boost::array<T,3> &q) : base(q[0], q[1], q[2]) {}
    };

    template<class C>
    inline std::ostream& operator<<(std::ostream &os, const Triplet<C> &triplet) {
	return os << "{"
		  << triplet[0] << "," << triplet[1] << "," << triplet[2]
		  << "}";
    }



    template<>
    struct Triplet<rysq::type> : triplet_base<rysq::type> {
	typedef triplet_base<rysq::type> base;
	Triplet(const type &a, const type &b, const type &c)
	    : base(a,b,c) {}
	template<typename T>
 	Triplet(const boost::array<T,3> &triplet) : base(triplet) {}
	size_t L() const { return detail::L::add(*this); }
	size_t size() const { return detail::size::multiply(*this); }
    };



    template<>
    struct Triplet<rysq::shell> : triplet_base<rysq::shell> {
	typedef triplet_base<rysq::shell> base;
	Triplet(const shell &a, const shell &b, const shell &c)
	    : base(a,b,c) {}
	template<typename T>
 	Triplet(const boost::array<T,3> &triplet) : base(triplet) {}
	size_t L() const { return detail::L::add(*this); }
	size_t size() const { return detail::size::multiply(*this); }
	size_t K() const { return detail::K::multiply(*this); }
    };



    template<>
    class Triplet<Shell> {
    public:
	typedef State Bra;
	typedef TwoC_State Ket;

	static int size(int a, int b, int c) {
	    return (shell(a).size*shell(b).size*shell(c).size);
	}

	static size_t size(const rysq::type (&triplet)[3]) {
	    return size(triplet[0], triplet[1], triplet[2]);
	}


	Triplet(const Shell &a, const Shell &b, const Shell &c)
	    : a_(a), b_(b),c_(c)
	{
	    initialize();
	}
	Triplet(const Triplet &triplet)
	    : a_(triplet[0]), b_(triplet[1]), c_(triplet[2])
	{
	    initialize();
	}
	~Triplet() {}
	int L() const { return L_; }
	int size() const { return size_; }
	int hybrid() const { return hybrid_; }
	int nc() const { return nc_; }
	int K() const { return K_; }
	const Shell& operator[](size_t i) const { return *shells_[i]; }
	const Bra bra() const { return Bra((*this)[0],(*this)[1]); }
	const Ket ket() const { return Ket((*this)[2]); }
	const Ket onec_ket() const { return Ket((*this)[2]); }
	operator Triplet<rysq::type>() const { return cast<rysq::type>(); }
	operator Triplet<rysq::shell>() const  { return cast<rysq::shell>(); }
	boost::array<Shell*,3> data() const { return shells_; }
    private:
	boost::array<Shell*,3> shells_;
	Shell a_, b_, c_;
	int L_, size_, hybrid_, nc_, K_;
	void initialize();
	void operator=(const Triplet&);
	template<typename T>
	Triplet<T> cast() const {
	    const Triplet &q = *this;
	    return Triplet<T>(T(q[0]), T(q[1]), T(q[2]));
	}
    };


} //namsepace rysq



#endif /* LIBCCHEM_RYSQ_SRC_TWOC_RYSQ_CORE_HPP_ */
