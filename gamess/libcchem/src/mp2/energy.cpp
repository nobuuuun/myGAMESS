#include "mp2/mp2.hpp"

#include "core/wavefunction.hpp"
#include "integrals/eri.hpp"

#include "runtime.hpp"
#include "parallel.hpp"
#include "thread.hpp"
#include "omp.hpp"
#include "exception.hpp"
#include "utility/progress.hpp"

#include "blas.hpp"
#if defined(HAVE_CUBLAS) && defined(CCHEM_INTEGRALS_ERI_CUDA)
#include "cublas.hpp"
#define CCHEM_MP2_CUDA
//#warning GPU MP2 disabled due to being slow
// #undef CCHEM_MP2_CUDA
#endif

#include "array/hdf5.hpp"
#include "utility/timer.hpp"

#include "cc/tensor.hpp"
#include "cc/utility.hpp"

#include <algorithm>
#include <iostream>
#include <memory>

#include <boost/noncopyable.hpp>
#include <boost/typeof/typeof.hpp>
//#include <boost/numeric/ublas/adaptors.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/thread.hpp>

#include "boost/utility/profiler.hpp"
#include <numeric>
namespace cchem {
namespace mp2 {
    typedef boost::numeric::ublas::matrix<
	double, boost::numeric::ublas::column_major> Matrix;
}
}

namespace cchem {
namespace mp2 {
namespace detail {

    namespace ublas = boost::numeric::ublas;
    namespace tensor =  cc::tensor;
    using cc::tensor_reference;
    using cc::tensor_const_reference;
    using cchem::Thread;

    struct Transform : boost::noncopyable {
    private:
	Thread thread_;

#ifdef CCHEM_MP2_CUDA
	struct Device : cchem::Device {
	    std::auto_ptr<integrals::Eri::Cuda> eri;
	    std::vector<cublas::handle_t> cublas_handle;
	    Device(size_t ns) : cchem::Device(ns) {
		foreach (Device::Stream &stream, this->streams) {
		    cublas::handle_t h = cublas::create();
		    cublas::set_stream(h, stream);
		    cublas_handle.push_back(h);
		}
	    }
	};
	std::auto_ptr<Device> device_;
#endif // CCHEM_MP2_CUDA

    private:
	boost::reference_wrapper<const Basis> basis_;
	boost::reference_wrapper<const Matrix> C_;
	integrals::Eri eri_;

    public:
	struct {
	    utility::timer::value_type eri, t1, t2;
	    utility::timer::value_type device[3];
	    utility::timer::value_type total;
	} profile;

	typedef tensor_reference<4> T4;

	static
	boost::array<size_t,4> memory(const Basis & basis, const Matrix &C) {
	    size_t maxs = basis.max().size();
	    size_t maxb = basis.max_block();
	    // size_t N = basis.size();
	    size_t M = maxs*maxs;
	    size_t n1 = C.size1();
	    boost::array<size_t,4> memory;
	    memory[0] = maxb*maxb*M;
	    memory[1] = n1*maxb*M;
	    memory[2] = n1*n1*M;
	    memory[3] = C.size1()*C.size2() + n1*maxs*M;

	    return memory;
	}

	Transform(boost::reference_wrapper<const Basis> basis,
		  boost::reference_wrapper<const Matrix> C,
		  boost::reference_wrapper<const integrals::Screening> screening,
		  const Runtime &rt)
	    : basis_(basis), C_(C),
	      eri_(basis_.get(), &screening.get()), profile()
	{

	    BOOST_AUTO(memory, this->memory(basis_.get(), C_.get()));
	    for (int i = 0; i < memory.size(); ++i) {
		thread_.data[i] = thread_.malloc(memory[i]);
	    }

#ifdef CCHEM_MP2_CUDA
	    std::vector<int> devices = rt.devices();
#pragma omp critical
	    if (omp::thread() < devices.size()) {
		try {

		    cuda::set_device(devices.at(omp::thread()));

		    BOOST_AUTO(const &basis, basis_.get());
		    BOOST_AUTO(const &C, C_.get());

		    device_.reset(new Device(16));
		    device_->eri.reset(new integrals::Eri::Cuda
				       (basis, &screening.get()));

		    size_t n = C.size1()*C.size2();
		    device_->S.set("C", device_->malloc(n),
				   boost::extents[1][1][C.size2()][C.size1()]);
		    cublas::set_vector(n, C.data().begin(),
				       device_->S["C"].data());

		    device_->S.set("t0\'", device_->malloc(memory[0]));
		    device_->S.set("t0\"", device_->malloc(memory[0]));
		    device_->S.set("t1", device_->malloc(memory[1]));
		    device_->S.set("t2", device_->malloc(memory[2]));
		    device_->S.set("h2", thread_.malloc(memory[2]));
			
		    foreach (Device::Stream &stream, device_->streams) {
			stream.S.set("u1", device_->malloc(memory[3]));
		    }
		}
		catch (const cuda::error &e) { // failed to allocate device/resources
		    std::cout << e.what() << std::endl;
		    device_.reset(NULL);
		}

	    }
#endif // CCHEM_MP2_CUDA

	}


	template<class A, class B>
	static void copy(A &a, ublas::range ra, const B &b, ublas::range rb,
			 blas::Context) {
	    ublas::range rows(0, a.size1());
	    noalias(project(a, rows, ra)) = project(b, rows, rb);
	}

#ifdef CCHEM_MP2_CUDA
	template<class A, class B>
	static void copy(A &a, ublas::range ra, const B &b, ublas::range rb,
			 cublas::handle_t handle) {
	    CCHEM_ASSERT(ra.size() == rb.size());
	    CCHEM_ASSERT(a.size1() == b.size1());
	    cublas::copy(handle,
			 (a.size1()*ra.size()),
			 b.data().begin() + b.size1()*rb.start(),
			 a.data().begin() + a.size1()*ra.start());
	}
#endif // CCHEM_MP2_CUDA


	/** (s,q,r,p) -> (i,s,q,r) */
	template<class C, class H>
	static
	void contract1(size_t np, size_t nq, size_t nr, size_t ns,
		       const std::vector<Basis::Shell> &shells,
		       const integrals::Eri::Quartets &quartets,
		       const C &C1, const double *G,
		       cc::tensor_reference<4> &u1,
		       std::map<size_t,size_t> &rs,
		       const std::vector< std::pair<H,double*> > &context) {
		    
	    BOOST_PROFILE_LINE;

	    using tensor::as_matrix;
	    using tensor::as_vector;

	    using ublas::make_matrix;
	    using ublas::column_major;

	    size_t n1 = C1.size1();
	    ublas::range r1(0,n1);

	    BOOST_AUTO(it, quartets.begin());
	    BOOST_AUTO(ctx, context.begin());

	    while (it != quartets.end()) {

		const H &handle = ctx->first;
		double *tmp = ctx->second;
		++ctx;
		if (ctx == context.end()) ctx = context.begin();

		size_t p = 0;
		size_t r = it->at(1);
		if (!rs.count(r)) {
		    size_t index = nr*rs.size();
		    rs[r] = index;
		}

		BOOST_AUTO(C_, (make_matrix<column_major>
				(C1.size1(), C1.size2(), tmp)));
		tmp += C1.size1()*C1.size2();
		    
		// pack C corresponding to same r shell
		{
		    //BOOST_PROFILE_LINE;
		    size_t nf = 0; // number of functions to pack
		    size_t f = shells[it->at(3)].start(); // first function
		    while (it->at(1) == r) {
			nf += np;
			++it;
			if (it == quartets.end()) break;
			size_t next = shells[it->at(3)].start();
			// new range to pack
			if (f+nf != next) {
			    copy(C_, ublas::range(p, p+nf),
				 C1, ublas::range(f, f+nf),
				 handle);
			    p += nf;
			    f = next;
			    nf = 0;
			}
		    }
		    copy(C_, ublas::range(p, p+nf),
			 C1, ublas::range(f, f+nf),
			 handle);
		    p += nf;
		}

		cc::tensor_const_reference<4> t0(G, boost::extents[p][nq][nr][ns]);
		G += p*nq*nr*ns;
		T4 t1(tmp, boost::extents[nq][nr][ns][n1]);
		BOOST_AUTO(t, (as_matrix<1,3>(t1)));
		{
		    //BOOST_PROFILE_LINE;
		    blas::gemm(1,
			       project(C_, r1, ublas::range(0,p)),
			       trans(as_matrix<3,1>(t0)),
			       0, t, handle);
		}

		{
		    // BOOST_PROFILE_LINE;
		    for (size_t k = 0, K = rs[r]; k < nr; ++k, ++K) {
			for (size_t j = 0; j < nq; ++j) {
			    BOOST_AUTO(t, as_vector(u1[K][j]));
			    blas::axpy(1, as_vector(t1[j][k]), t, handle);
			}
		    }
		}

	    }
	}


	/** (i,s,q,r) -> (j,i,s,q) */
	template<class C, class H>
	static 
	void contract2(size_t n1, size_t nq, size_t nr, size_t ns,
		       const std::vector<Basis::Shell> &shells,
		       std::map<size_t,size_t> &rs,
		       const C &C2,
		       const cc::tensor_reference<4> &t1,
		       cc::tensor_reference<4> &t2,
		       std::pair<H,double*> context) {

	    BOOST_PROFILE_LINE;

	    using ublas::make_matrix;
	    using ublas::column_major;
	    using tensor::as_matrix;

	    if (rs.empty()) return;

	    const H &handle = context.first;
	    double *tmp = context.second;

	    size_t n2 = C2.size1();
	    ublas::range r2(0,n2);
	    BOOST_AUTO(C_, (make_matrix<column_major>(C2.size1(), C2.size2(), tmp)));

	    BOOST_AUTO(it, rs.begin());
	    while (it != rs.end()) {
		size_t r = it->second;
		copy(C_, ublas::range(r, r+nr),
		     C2, shells[it->first].range(),
		     handle);
		++it;
	    }
	    ublas::range rr(0, nr*rs.size());
	    tensor_const_reference<4>
		u1(t1.data(), boost::extents[rr.size()][nq][ns][n1]);
	    BOOST_AUTO(t, (as_matrix<1,3>(t2)));
	    blas::gemm(1,
		       project(C_, r2, rr),
		       trans(as_matrix<3,1>(u1)),
		       1, t, handle);
	}

	/** evaluate (j,i,s,q) */
	tensor_const_reference<4> operator()(int q, int s) {   //declaration for tensor_const_reference<4> is set to double tensor_const_reference<4,double>
	    BOOST_PROFILE_LINE;

	    utility::timer timer;

	    namespace ublas = boost::numeric::ublas;
	    namespace tensor = cc::tensor;
	    using tensor::as_matrix;

	    typedef cc::tensor_reference<4> T;

	    BOOST_AUTO(const &basis, basis_.get());
	    BOOST_AUTO(const &C, C_.get());

	    const Basis::Shell &Q = basis.shells().at(q);
	    const Basis::Shell &S = basis.shells().at(s);

	    // size_t N = basis.size();
	    size_t n1 = C.size1();  //number of correlated orbitals (occupied - core)
	    ublas::range r1(0,n1);

	    struct {
		std::vector< std::pair<blas::Context, double*> > host;
#ifdef CCHEM_MP2_CUDA
		std::vector< std::pair<cublas::handle_t, double*> > device;
#endif // CCHEM_MP2_CUDA
	    } context;

	    {
		std::pair<blas::Context, double*> ctx;
		ctx.first = blas::Context();
		ctx.second = thread_.data[3];  //thread_.date[3] memory : C.size1()*C.size2() + n1*maxs*M
		context.host.push_back(ctx);
	    }

#ifdef CCHEM_MP2_CUDA
	    Device *device = this->device_.get();
	    if (device) {
		BOOST_PROFILE_LINE;
		device->S.reset("h2", boost::extents[Q.size()][S.size()][n1][n1]);
		device->S.reset("t2", boost::extents[Q.size()][S.size()][n1][n1]);
		cublas::clear(device->S["t2"].num_elements(),
			      device->S["t2"].data());
		for (size_t i = 0; i < device->streams.size(); ++i) {
		    std::pair<cublas::handle_t, double*> ctx;
		    ctx.first = device->cublas_handle.at(i);
		    ctx.second = device->streams.at(i).S["u1"].data();
		    context.device.push_back(ctx);
		}
	    }
#endif // CCHEM_MP2_CUDA

	    // (j,i,s,q)
	    T4 t2(thread_.data[2], boost::extents[Q.size()][S.size()][n1][n1]); //thread_.data[2] memory : n1*n1*M
	    std::fill_n(t2.data(), t2.num_elements(), 0);

	    integrals::Eri::Quartets quartets;

	    BOOST_PROFILE_LINE;

	    foreach (const Basis::Block &R, basis.blocks()) {

		size_t nr = R.shell().size();  //nr the number of shells in current loop
		size_t nq = Q.size();
		size_t ns = S.size();


		T4 u1(thread_.data[1], boost::extents[R.range().size()][nq][ns][n1]); // thread_.data[1] memory : n1*maxb*M : M -> maxshell**2
		std::fill_n(u1.data(), u1.num_elements(), 0);

#ifdef CCHEM_MP2_CUDA
		if (device) {
		    BOOST_PROFILE_LINE;
		    cuda::synchronize();
		    device->S.reset
			("t1", boost::extents[R.range().size()][nq][ns][n1]);
		    cublas::clear(device->S["t1"].num_elements(),
				  device->S["t1"].data());
		}
#endif // CCHEM_MP2_CUDA

		BOOST_AUTO(const &shells, basis.shells());

		// record non-zero r shell and its index
		struct {
		    std::map<size_t,size_t> host, device;
		} rs;

		BOOST_PROFILE_LINE;

		foreach (const Basis::Block &P, basis.blocks()) {

		    BOOST_PROFILE_LINE;

		    utility::timer timer;

		    size_t np = P.shell().size();

		    quartets.clear();
		    foreach (const Basis::Shell &r, R) {
			eri_.screen(S, r, Q, P, quartets);
		    }

#ifdef CCHEM_MP2_CUDA
		    if (device) {
			try {
			    BOOST_PROFILE_LINE;
			    BOOST_AUTO(&symbol, device->S);
			    (*device->eri)(S, R, Q, P, quartets,
					   symbol["t0\'"].data());
			    // wait for integrals to compute
			    // wait for previous transformation
			    cuda::synchronize();
			    symbol.swap("t0\'", "t0\"");
			    BOOST_PROFILE_LINE;
			    contract1(np, nq, nr, ns, shells, quartets,
				      as_matrix<1,3>(symbol["C"]),
				      symbol["t0\""].data(), symbol["t1"],
				      rs.device,
				      context.device);
			    continue;
			}
			catch (...) {}
		    }
#endif // CCHEM_MP2_CUDA

		    BOOST_PROFILE_LINE;

		    timer.reset();

		    double *G = thread_.data[0];    // thread_.data[0] memory : maxb*maxb*M
		    eri_(S, R, Q, P, quartets, G);
		    profile.eri += timer;

		    contract1(np, nq, nr, ns, shells, quartets,
		    	      C, G, u1, rs.host,
		    	      context.host);

		    profile.t1 += timer;
		} // P

		BOOST_PROFILE_LINE;
		    
		timer.reset();
 
#ifdef CCHEM_MP2_CUDA
		if (device) {
		    cuda::synchronize();
		    BOOST_AUTO(&symbol, device->S);
		    contract2(n1, nq, nr, ns, shells,
			      rs.device,
			      as_matrix<1,3>(symbol["C"]),
			      symbol["t1"], symbol["t2"],
			      context.device.front());
		    cuda::synchronize();
		}
#endif // CCHEM_MP2_CUDA

		{
		    BOOST_PROFILE_LINE;
		    contract2(n1, nq, nr, ns, shells, rs.host,
			      C, u1, t2,
			      context.host.front());
		}

		profile.t2 += timer;

	    } // R

#ifdef CCHEM_MP2_CUDA
	    if (device) {
		BOOST_PROFILE_LINE;
		cuda::synchronize();
		cuda::copy(device->S["t2"].num_elements(),
			   device->S["t2"].data(),
			   device->S["h2"].data(),
			   cuda::device_to_host);
		using tensor::as_vector;
		BOOST_AUTO(t, as_vector(t2));
		blas::axpy(1, as_vector(device->S["h2"]), t);
	    }
#endif // CCHEM_MP2_CUDA

	    return t2;

	}

    };

    struct energy {
    private:
    	boost::reference_wrapper<const Wavefunction> wf_;
	Matrix Cv_;
	double *buffer_[2]; 

    public:

    /**
    *   @brief Closed shell mp2 (called Closed shell-like term for ZAPT)
    *
    *   @author Andrey Asadchev
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param E energy component
    *
    *   @date 12-22-2014 Luke Roskop
    *   -added ns number of singly occupied orbitals (ns=0 for closed shell mp2)
    *
    */
	template<class T, class V>
	static double evaluate(const T &t, double dij, const V &e, const size_t &ns) {
	    double E = 0;
	    for (size_t b = ns; b < t.size2(); ++b) {
		for (size_t a = ns; a < t.size1(); ++a) {
		    double ab = t(a,b);
		    double ba = t(b,a);
		    double de = dij - (e[a] + e[b]);
		    // E += ab*(2*ab - ba)/de;
		    E += ab*(2*ab - ba)/de;
		}//a
	    }//b
	    return E;
	}
    /**
    *   @brief singly unoccupied term
    *   
    *   @details \sum_{\substack_ijsb}\frac{v_i_j^s^b(2v_i_j^s^b-v_j_i^s^b)}{\varepsilon_i+\varepsilon_j-\varepsilon_s^--\varepsilon_b}
    *
    *   @author Luke Roskop
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_m correction vector to singly occupied orbitals
    *   @param E energy component
    */
	template<class T, class V>
	static double evaluate2(const T &t, double dij, const V &e, size_t &ns, std::vector<double> &e_m) {
	    double E = 0;
	   	    for (size_t s = 0; s < ns; ++s) {
	   		for (size_t a = ns; a < t.size1(); ++a) {
			double as = t(a,s);
		    double sa = t(s,a);
		    double de = dij - (e[a] + e[s]+ e_m[s]);
		    E+= (as*as - as*sa + sa*sa)/de;
	   		}//a
	    }//s
	    return E; // /(double)2;
	}
    /**
    *   @brief singly occupied
    *
    *   @author Luke Roskop
    *
    *   @details \sum_{\substack_sjab}\frac{v_s_j^a^b(2v_s_j^a^b-v_j_s^a^b)}{\varepsilon_s^++\varepsilon_j-\varepsilon_a-\varepsilon_b}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_m correction vector to singly occupied orbitals
    *   @param E energy component
    */
	template<class T, class V>
	static double evaluate6(const T &t, double dij, const V &e, size_t &ns, double &e_m) {
	    double E = 0;
	   	    for (size_t b = ns; b < t.size1(); ++b) {
	   		for (size_t a = ns; a < t.size1(); ++a) {
			double ab = t(a,b);
		    double ba = t(b,a);
		    double de = dij - (e[a] + e[b])-e_m;
		    E += (ab*(2*ab - ba))/de;
	   		}//a
	    }//b
	    return E/(double)2;
	}
    /**
    *   @brief singly unoccupied/occupied
    *
    *   @author Luke Roskop
    *
    *   @details \sum_{\substack_sjbt}\frac{(v_s_j^b^tv_s_j^b^t)}{\varepsilon_s^++\varepsilon_j-\varepsilon_b-\varepsilon_t^-}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_s correction vector to singly occupied orbitals
    *   @param e_a energy correction to orbital singly occupied orbital
    *   @param E energy component
    */
	template<class T, class V>
	static double evaluate4(const T &t, double dij, const V &e, size_t &ns, double &e_a, std::vector<double> &e_s) {
	    double E = 0;
	   	    for (size_t s = 0; s < ns; ++s) {
	   		for (size_t a = ns; a < t.size1(); ++a) {
			double ab = t(s,a);
		    double de = dij - (e[a] + e[s]) - e_a - e_s[s];
			E += ab*ab/de;
	   		}//a
	    }//s
	    return E/(double)2;
	}
    /**
    *   @brief 2 singly occupied
    *
    *   @author Luke Roskop
    *
    *   @details \frac{1}{4} \sum_{\substack_stab}\frac{(v_s_t^a^b-v_t_s^a^b)(v_s_t^a^b-v_t_s^a^b)}{\varepsilon_s^++\varepsilon_t^+-\varepsilon_a-\varepsilon_b}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_s energy correction to orbital singly occupied orbital
    *   @param e_t energy correction to orbital singly occupied orbital
    *   @param E energy component
    */
	template<class T, class V>
	static double evaluate3(const T &t, double dij, const V &e, size_t &ns, double &e_s, double &e_t) {
	    double E = 0;
	    for (size_t b = ns; b < t.size1(); ++b) {
	    	for (size_t a = ns; a < t.size1(); ++a) {
	    		double ab = t(a,b);
	    		double ba = t(b,a);
	    		double de =dij -(e[a] + e[b]) - e_s - e_t;
	    		E += ((ab-ba)*(ab-ba))/de;
 	   		}//a
	    }//b
	    return E/(double)4;
	}

    /**
    *   @brief doubly unoccupied
    *
    *   @author Luke Roskop
    *
    *   @details \frac{1}{4} \sum_{\substack_ijst}\frac{(v_i_j^s^t-v_j_i^s^t)(v_i_j^s^t-v_j_i^s^t)}{\varepsilon_i+\varepsilon_j-\varepsilon_s^--\varepsilon_t^-}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_m correction vector to singly occupied orbitals
    *   @param E energy component
    */
	template<class T, class V>
	static double evaluate7(const T &t, double dij, const V &e, size_t &ns, std::vector<double> &e_m) {
	    double E = 0;
	    for (size_t b = 0; b < ns; ++b) {
	    	for (size_t a = 0; a < ns; ++a) {
	    		double ab = t(a,b);
	    		double ba = t(b,a);
	    		double de =dij -(e[a] + e[b]) -e_m[a] -e_m[b];
	    		E += ((ab-ba)*(ab-ba))/de;
	    	}//a
	    }//b
	    return E/(double)4;
	}

    /**
    *   @brief initialization
    *
    *   @author Andrey Asadchev
    *
    *   @param wf wave function
    *   @param memory memory
    *
    *   @date 12/22/2014 Luke Roskop
    *   -redimension Cv_ to include singly occupied orbitals (if any)
    */
	energy(boost::reference_wrapper<const Wavefunction> wf,
	       Runtime::Memory &memory)
	    : wf_(wf)
	{
	    size_t N =  wf_.get().basis().size();
	    Cv_ = wf.get().C(wf.get().single_virtual());
	    // thread buffers
	    buffer_[0] = memory.malloc<double>(N*N);
	    buffer_[1] = memory.malloc<double>(N*N);
	}

    /**
    *   @brief overload () operator to evaluate MBPT
    *
    *   @author Andrey Asadchev
    *
    *   @param ij current ij pair
    *   @param *data pointer to half-transformed integrals
    *   @param B block size
    *   @param b block index
    *	@param pfock array of points to need for ZAPT
    *	@param e_term store energy components
    *	@param e_m vector needed to modify singly occupied orbital energies
    *
    *   @date 12/22/2014 Luke Roskop
    *   -redefine no, nv, ev, and ev to reflect singly occupied orbitals
    *   -add ns, nd, es, and ed
    *   -add ZAPT energy terms (will only be computed if there are single occupied orbitals)
    */
	void operator()(std::pair<size_t,size_t> ij,
		  const double *data,
		  size_t B, size_t b, std::vector<double *> pfock, double (&e_term)[7],
		  std::vector<double> &e_m) {

    BOOST_PROFILE_LINE;

    namespace ublas = boost::numeric::ublas;
    using ublas::make_matrix;
    using ublas::column_major;

    BOOST_AUTO(const &wf, wf_.get());

    //wf.set_single_and_virtuals(check1,check2);

    BOOST_AUTO(const &ea, wf.e(wf.active()));
    BOOST_AUTO(const &ev, wf.e(wf.single_virtual()));
    BOOST_AUTO(const &es, wf.e(wf.single())); 	
    BOOST_AUTO(const &ed, wf.e(wf.double_occ()));
    BOOST_AUTO(const &shells, wf.basis().shells());

    size_t N =  wf.basis().size();

    size_t no = wf.double_occ().size();
    size_t nv = wf.single_virtual().size(); //includes the number of single occupied orbital (mp2 -> 0)

    size_t ns = wf.single().size();  
    size_t nd = wf.double_occ().size();

    double *t[5];
    // transformation buffers
    t[2] = buffer_[0];  //	    buffer_[0] = memory.malloc<double>(N*N);
    t[3] = buffer_[1];  //	    buffer_[1] = memory.malloc<double>(N*N);
    t[4] = buffer_[0];  //	    buffer_[0] = memory.malloc<double>(N*N);

    BOOST_AUTO(t2, make_matrix<column_major>(N, N, t[2]));
    BOOST_AUTO(t3, make_matrix<column_major>(nv, N, t[3]));
    BOOST_AUTO(t4, make_matrix<column_major>(nv, nv, t[4]));

    // restore (B,qs) to (Q,S,B)
    const double *ptr = data+b;
    for (size_t s = 0; s < shells.size(); ++s) {
	for (size_t q = 0; q <= s; ++q) {
	    const Basis::Shell &Q = shells.at(q);
	    const Basis::Shell &S = shells.at(s);
	    size_t size = Q.size()*S.size();
	    if (q == s) { // diagonal block
		copy(Q, S, ptr, B, t2);
		ptr += size*B;
	    }
	    else {
		copy(Q, S, ptr, B, t2);
		ptr += size*B;
		copy(S, Q, ptr, B, t2);
		ptr += size*B;
	    }
	}
    }

    blas::gemm(1, trans(Cv_), t2, 0, t3); //t3 = trans(Cv_) * t2
    blas::gemm(1, t3, Cv_, 0, t4); // t4 = t3 * Cv_

    size_t i = ij.first;
    size_t j = ij.second;
    double dij = (ea[i] + ea[j]);

//closed shell HF term (called closed shell-like term when performing ZAPT).
//	If closed shell MBPT is being performed, only the following function is performed,
// 	as the MBPT2 correlation energy is accumulated into e_term[0] (formally called E)
if(ij.second < nd)	e_term[0] += ((1+(i!=j))*evaluate(t4, dij, ev, ns));

//if this is closed shell MBPT then exit (i.e. ns == 0)
if(ns == 0) return;

//Everything below is for ZAPT
//closed shell-like term E1 / singly unoccupied E2 / two singly unoccupied E7
if(ij.second < nd)		{
	e_term[1] += ((1+(i!=j))*evaluate2(t4, dij, ev, ns, e_m));
    e_term[6] += ((1+(i!=j))*evaluate7(t4, dij, ev, ns, e_m));
	}

//singly occupied E6 / singly occupied/unoccupied E4
if(ij.second > nd-1 && ij.first < nd)	    {
	e_term[5] += ((1+(i!=j))*evaluate6(t4, dij, ev, ns, e_m[ij.second-nd]));
	e_term[3] += ((1+(i!=j))*evaluate4(t4, dij, ev, ns, e_m[ij.second-nd], e_m));
	}

//two singly occupied E3
if(ij.first > nd-1)	    e_term[2] += ((1+(i!=j))*evaluate3(t4, dij, ev, ns, e_m[ij.first-nd], e_m[ij.second-nd]));

//Fock matrix contribution E5 (one thread at a time)
#pragma omp critical
if(ij.second >= nd && ij.first < nd){
	BOOST_AUTO(const &w, row(t4, ij.second-nd));
	double *p = pfock.at(ij.first);
	for (int a = ns; a < w.size(); a++){*p += w[a]; p++;}
}


}//end

    private:
	template<class R, typename U, class T>
	void copy(const R &q, const R &s, U *data, size_t B, T &t) {
	    for (size_t j = s.start(); j < s.start()+s.size(); ++j) {
		for (size_t i = q.start(); i < q.start()+q.size(); ++i) {
		    t(i,j) = *data;
		    data += B;
		}
	    }
	}

    };

    struct Block {
	explicit Block(size_t s) : s(s), q(0,0), size(0), index(0) {}
	size_t s;  
	boost::numeric::ublas::range q;
	size_t size, index, progress;
	bool operator<(const Block &b) const {
	    return (this->size < b.size);
	}
    };

    struct async {
        /**
        *   @brief launch boost thread for asynchronous fetch of half-transformed integrals
        *
        *   @author Andrey Asadchev
        *
        *   @param k beginning of block
        *   @param m end of block
    	*   @param n length each ij
    	*   @param V reference to half-transformed integrals
    	*   @param *data pointer to fetched half-transformed integrals
    	*
    	*   @date 12-22-2014 Luke Roskop
    	*   -added k variable
        */
    	void get(size_t k, size_t m, size_t n, size_t b,
    		 const Array<double> &V, double *data) {
    	    thread_.reset(new boost::thread(launch, k, m, n, b, boost::cref(V), data));
    	}
    	void wait() {
	    thread_->join();
    	}
    	private:
    	std::auto_ptr<boost::thread> thread_;
    	/**
    	 *   @brief fetch half-transformed integrals
    	 *
    	 *   @author Andrey Asadchev
    	 *
    	 *   @param k beginning of block
    	 *   @param m end of block
    	 *   @param n length each ij
    	 *   @param V reference to half-transformed integrals
    	 *   @param *data pointer to fetched half-transformed integrals
    	 *
    	 *   @date 12-22-2014 Luke Roskop
    	 *   -added k variable
    	 */
    	static
		void launch(size_t k, size_t m, size_t n, size_t b,
				const Array<double> &V, double *data) {
    		BOOST_PROFILE_LINE;
    		utility::timer timer;
    		size_t start[] = { k, 0, b };
    		size_t finish[] = { m, n, b+1 };
    		V.get(data, start, finish);
    		//std::cout << "I/O: " << (N*N*B*8)/(timer*1e6) << std::endl;
    	}
    };

    struct openshell {
    private:
    	boost::reference_wrapper<const Wavefunction> wf_;
    	size_t nd, nv, ns, N;
    	size_t B_;
    	double *buffer_[2];
    	std::vector<double *> pfock_;
    	std::vector<double> e_m;

    public:

    	/**
    	 *   @brief openshell constructor
    	 *
    	 *   @author Luke Roskop
	 *
	 *   @details this routine creates arrays 
	 *       -e_m stores exchange integrals that are used to modify singly occupied orbital energies
         *       -pfock is an array of (nd= #doubly_occ) pointers. Each pointer designates a vector
         *          of length nv (# virtuals). These vectors store {-\sum_s v_i_a^s^s} type contractions
         *          in which the individual components are accumulated during evaluation of the PT energy 
         *          (due to the way in which integrals are stored).
    	 *
    	 *   @param wf wave function specifications
    	 *   @param memory adaptor
    	 *   @param B block size
    	 */
	openshell(boost::reference_wrapper<const Wavefunction> wf,
		  Runtime::Memory &memory,
		  size_t B): wf_(wf), B_(B)
    {
    	    ns = wf_.get().single().size();
    	    if(ns == 0)return;
    	    N = wf_.get().basis().size();
    		nd = wf_.get().double_occ().size();
    		nv = wf_.get().virtuals().size();
    		for (size_t i = 0; i < nd; i++){
    			pfock_.push_back(memory.malloc<double>(nv));
    			double *p = pfock_.back();
    			std::fill(p, p + nv, 0); //initialize vector
    		}
    		for (int i = 0; i < ns; ++i) e_m.push_back(0);

    };//openshell constructor

    	/**
    	 *   @brief get pfock vector (self-documenting)
    	 *
    	 *   @author Luke Roskop
    	 */
    	std::vector<double *> get_pfock() {return pfock_;};

    	/**
    	 *   @brief e_m vector (self-documenting)
    	 *
    	 *   @author Luke Roskop
    	 */
    	std::vector<double> & get_e_m() {return e_m;};

    	/**
    	 *   @brief compute fock-like energy
    	 *
    	 *   @author Luke Roskop
    	 *
         *   @details \frac{1}{2} \sum_{\substack_ia}\frac{K_i^aK_i^a}{\varepsilon_i-\varepsilon_a}
         *
    	 *   @param pe reference to parallel environment
    	 */
    	double compute_fock_like_term(Parallel &pe){
    		for (size_t i = 0; i < nd; i++)
    			pe.reduce("+",pfock_[i],nv);
    		//only the master process (rank 0) continues on
    		//(so Fock-like energy is summed correctly)
    		if(pe.rank() != 0) return 0;
    		BOOST_AUTO(const &wf, wf_.get());
    		BOOST_AUTO(const &ev, wf.e(wf.virtuals()));
    		BOOST_AUTO(const &ed, wf.e(wf.double_occ()));
    		double E = 0;
    		for (size_t d = 0; d < wf.double_occ().size(); d++){
    			double *p = pfock_.at(d); 
    			for (size_t a = 0; a < wf.virtuals().size(); a++){
    				E += pow((*p),2)/(ed[d]-ev[a]);
    				p++;
    			}//a
    		}//d
    		return E/(double)2;
    	};//compute_fock_like_term

    	/**
    	 *   @brief computes exchange integrals
    	 *
    	 *   @author Luke Roskop
    	 *
    	 *   @param V pointer to half-transformed integral array
    	 *   @param data2 pointer to scratch data
    	 *   @param index index for block data
    	 *   @param memory memory adaptor
    	 */
    	void build_exchange_integrals(const Array<double> *V, double *data2,
    			std::vector< std::pair<size_t,size_t> > index,
    			Runtime::Memory &memory){
//    		std::cout << memory.available() << " " << memory.used() << std::endl;
    		buffer_[0] = memory.malloc<double>(N*N);
    		buffer_[1] = memory.malloc<double>(N*N);
//    		std::cout << memory.available() << " " << memory.used() << std::endl;
    		BOOST_AUTO(const &wf, wf_.get());
    		detail::async async;
    		size_t nd = wf.double_occ().size();
    		Matrix Cv_ = wf.C(wf.single_virtual());
    		for(size_t s = nd+1; s <= nd+ns; s++){
    			size_t start = s*(s-1)/2;
    			for(size_t t = nd+1; t <= s ;t++){
    				size_t st = start+t-1; //find position of st pair
    				size_t position = st/B_; //find which block contains the st pair
    				size_t st_get = st - B_*position; //adjust to reflect position within the current block

    				async.get(st_get,st_get+1, N*N, position, *V, data2);
    				async.wait();
    				this->populate_e_m(index.at(st), data2, Cv_);
    			}//t
    		}//s
//    		std::cout << memory.available() << " " << memory.used() << std::endl;
    		memory.free(buffer_[1]); //free memory
    		memory.free(buffer_[0]); //free memory
//    		std::cout << memory.available() << " " << memory.used() << std::endl;
    	};//build_exchange_integrals

    	/**
    	 *   @brief build energy corrections for singly occupied orbitals (this is ZAPT)
    	 *
    	 *   @author Luke Roskop
    	 *
    	 *   @param index index for block data
    	 *   @param data pointer to scratch data
    	 *   @param Cv_ local copy of LCAO coefficients
    	 */
    	void populate_e_m(std::pair<size_t,size_t> ij,
    			const double *data,
				Matrix &Cv_) {

    		namespace ublas = boost::numeric::ublas;
    		using ublas::make_matrix;
    		using ublas::column_major;

    		BOOST_AUTO(const &wf, wf_.get());

    		BOOST_AUTO(const &shells, wf.basis().shells());

    		size_t ns = wf.single().size();
    		size_t nd = wf.double_occ().size();

    		double *t[5];
    		// transformation buffers
    		t[2] = buffer_[0];  //	    buffer_[0] = memory.malloc<double>(N*N);
    		t[3] = buffer_[1];  //	    buffer_[1] = memory.malloc<double>(N*N);
    		t[4] = buffer_[0];  //	    buffer_[0] = memory.malloc<double>(N*N);

    		BOOST_AUTO(t2, make_matrix<column_major>(N, N, t[2]));
    		BOOST_AUTO(t3, make_matrix<column_major>(ns, N, t[3]));
    		BOOST_AUTO(t4, make_matrix<column_major>(ns, ns, t[4]));

    		// restore (B,qs) to (Q,S,B)
    		const double *ptr = data;
    		for (size_t s = 0; s < shells.size(); ++s) {
    			for (size_t q = 0; q <= s; ++q) {
    				const Basis::Shell &Q = shells.at(q);
    				const Basis::Shell &S = shells.at(s);
    				size_t size = Q.size()*S.size();
    				if (q == s) { // diagonal block
    					copy(Q, S, ptr, 1, t2);
    					ptr += size;
    				}
    				else {
    					copy(Q, S, ptr, 1, t2);
    					ptr += size;
    					copy(S, Q, ptr, 1, t2);
    					ptr += size;
    				}
    			}
    		}

    		blas::gemm(1, trans(Cv_), t2, 0, t3); //t3 = trans(Cv_) * t2
    		blas::gemm(1, t3, Cv_, 0, t4);        // t4 = t3 * Cv_

    		size_t i = ij.first;
    		size_t j = ij.second;

    		e_m[j-nd] += t4(j-nd,i-nd)/2;
    		if(i != j)e_m[i-nd] += t4(j-nd,i-nd)/2;

    	};//populate_e_m

    	/**
    	 *   @brief sum and finish computing energy components
    	 *
    	 *   @author Luke Roskop
    	 *
    	 *   @param e_term thread public energy component vector
    	 *   @param private_e_term thread private energy component vector
    	 *	 @param pe
    	 */
    	void operator()(double (&e_term)[7],double (&private_e_term)[7], Parallel &pe){
#pragma omp critical // accumulate ZAPT energy components
    		for (int i = 0; i < 7; i++) e_term[i] += private_e_term[i];
#pragma omp barrier  //make sure threads are done contributing to pfock[i]->vector
#pragma omp master
    		e_term[4] = this->compute_fock_like_term(pe);
    	};//operator()

    private:
    	template<class R, typename U, class T>
    	void copy(const R &q, const R &s, U *data, size_t B, T &t) {
    		for (size_t j = s.start(); j < s.start()+s.size(); ++j) {
    			for (size_t i = q.start(); i < q.start()+q.size(); ++i) {
    				t(i,j) = *data;
    				data += B;
    			}
    		}
    	};//copy

    };//openshell
}//cchem
}//mp2

/**
*   @brief mp2/zapt energy evaluation
*
*   @author Andrey Asadchev
*
*   @param wf wavefunction
*   @param rt runtime
*
*/

//namespace ublas = boost::numeric::ublas;
double mp2::energy(Wavefunction wf, Runtime &rt) {

    BOOST_PROFILE_REGISTER_THREAD;

    namespace ublas = boost::numeric::ublas;
    using ublas::make_vector;
    using ublas::make_matrix;
    using ublas::column_major;

    using cc::Tensor;
    using cc::tensor_reference;
    namespace tensor = cc::tensor;
    using tensor::as_matrix;

    using detail::Transform;

    BOOST_AUTO(const &basis, wf.basis());
    BOOST_AUTO(const &shells, basis.shells());

    wf.sort();
    // wf.reverse();
    double cutoff = rt.get<double>("/mp2/integrals/cutoff", 1e-10);
    integrals::Screening screening(basis, cutoff);

    Matrix Ca = trans(wf.C(wf.active()));

    size_t N = basis.size();
    size_t no = wf.active().size();
    size_t nv = wf.virtuals().size();
    ublas::range ro(0,no), rv(0,nv), ra(0,N);

    Parallel pe;
    // suppress output
    if (pe.rank() != 0) rt.set_cout(0);
 
    {
	BOOST_AUTO(cout, rt.cout());
	cout << "active: " << no << std::endl;
	cout << "virtual: " << nv << std::endl;
	cout << "atomic: " << N << std::endl;
	cout << "eri cutoff: " << screening.value() << std::endl;
    }

    Runtime::Memory &memory = rt.memory();
    double used = 0;

    size_t max_threads = omp_get_max_threads();
    size_t nwords = (memory.available()/sizeof(double));

    // runtime may be different, ensure consistency
    pe.broadcast(&max_threads, 1, 0);
    pe.broadcast(&nwords, 1, 0);

    {
	BOOST_AUTO(const &memory, Transform::memory(basis, Ca));
	size_t required = 0;
	foreach (size_t n, memory) required += n;
	required = std::max(required, 2*N*N);
	required *= max_threads;
	if (required > nwords)
	    throw cchem::exception("not enough memory to run MP2");
	nwords -= required;
    }

    const size_t no2 = (no*no+no)/2;
    // block factor
    size_t B = 1;

    while (N*N*(B+1) < nwords/2 && B < no2) {++B;} // first dimension

    size_t nb = (no2+B-1)/B;

    while (B > 1 && nb%pe.size()) {
	--B; nb = (no2+B-1)/B;
    }

    while (B > 1 && (no2+B)/(B-1) == nb) --B; // shrink first dimension
    CCHEM_ASSERT(no2 <= B*nb);

    // qs shell block  MAXQS is the maximum block size
    const size_t MAXQS = std::min<size_t>(nwords/(nb*B*max_threads), N*basis.max().size());

    if (MAXQS < 2*basis.max().size()*basis.max().size())   //basis.max().size() size of largest shell
	throw cchem::exception("not enough memory to run MP2");

    {
	size_t dims[3] = { B, N*N, nb };
	size_t chunk[3] = { dims[0], dims[1], 1 };
	// 32 MB chunks limit  note: 2**25/1024/1024/1024 = 0.03125 GB (31.3 MB)
	chunk[0] = std::min((1LU << 25)/(B*sizeof(double)), chunk[0]);
	// pad to chunk size
	// dims[0] = ((dims[0] + chunk[0] - 1)/chunk[0])*chunk[0];
	rt.arrays().allocate<double,3>("mp2.v(b,qs,ij)", dims, pe, chunk);
    }
    Array<double> *V = rt.arrays().find< Array<double> >("mp2.v(b,qs,ij)");
    rt.cout() << *V << std::endl;

    struct {
	utility::timer time;
	utility::timer::value_type eri, t1, t2;
	utility::timer::value_type device[3];
	utility::timer::value_type transform;
	utility::timer::value_type input, output, reduce;
	size_t size;
    } profile = {};


    // create a vector of Block objects
    //	each Block contains information about how many basis functions pairs there are for the shells q and s
    //  the range of q shells being coupled to an s shell
    using detail::Block;
    std::vector<Block> blocks;
    for (size_t s = 0, index = 0; s < shells.size(); ++s) {      //loop over shells
	const Basis::Shell &S = shells.at(s);
	// QS blocks
	blocks.push_back(Block(s));
	blocks.back().index = index;
	for (size_t q = 0; q <= s; ++q) {
	    typedef ublas::range range;
	    const Basis::Shell &Q = shells.at(q);
	    size_t nqs = Q.size()*S.size()*((q == s) ? 1 : 2);  //if shells are differnt, size is multiplied by two
	    CCHEM_ASSERT(nqs <= MAXQS);      //MAXQS is the (number of basis functions * the largest shell) (ie for f - 10 basis functions or
	    								 // MAXQS = (nwords/(nb*B*max_threads)) -- whichever is smaller

	    if (blocks.back().size + nqs > MAXQS) {  //if QS is getting large, start a new block
		blocks.push_back(Block(s));
		blocks.back().q = range(q, q);
		blocks.back().index = index;

	    }
	    int start = blocks.back().q.start();
	    blocks.back().q = range(start, q+1);
	    blocks.back().size += nqs;
	    index += nqs;
	}
    }

    // permute blocks
    {
	std::reverse(blocks.begin(), blocks.end());
    	std::vector<Block> p;

	size_t N = std::max<size_t>((blocks.size()/pe.size()*max_threads), 1);  //pe.size is the number of procs

    	for (size_t i = 0; i < blocks.size(); ++i) {         // blocks.size() is the number of Block(s) in the blocks vector
    	    size_t M = blocks.size() - blocks.size()%N;      // blocks.size()%N is the remainder when blocks.size is divided by N
    	    size_t j = (i < M) ? (i/N + (i%N)*(M/N)) : i;
    	    p.push_back(blocks.at(j));
    	}
    	blocks = p;

    }

    //  progress reporting
    {
	size_t progress = 0;
	foreach (Block &b, blocks) {
	    progress += b.size;
	    b.progress = progress;
	}
    }

    utility::Progress progress;
    {
	BOOST_AUTO(cout, rt.cout());
	//cout << rt << std::endl;
	cout << "OpenMP threads: " << max_threads << std::endl;
	cout << "blocks: " << blocks.size() << std::endl;
    }
    if (pe.rank() == 0) progress.reset(N*N);

    pe.task().reset();

#pragma omp parallel
    if (pe.node().rank() == 0) {

	BOOST_PROFILE_LINE;

	detail::Thread::Task<Parallel::Task&> task(pe.task());

	double *buffer = memory.malloc<double>(MAXQS*nb*B);
#pragma omp barrier

	//initialize transformation, allocate memory
	Transform transform(boost::cref(basis), boost::cref(Ca),
			    boost::cref(screening), rt);

	while (++task < blocks.size()) {

	    Block block = blocks.at(task);
	    size_t s = block.s;
	    const Basis::Shell &S = shells.at(s);

	    boost::multi_array_ref<double,3>
	    	u(buffer, boost::extents[nb][block.size][B]);   //nb is related to blocking, block.size is the number of SQ pairs in the block (limit - MAXQS), B is block parameter
	    std::fill(u.data(), u.data()+B*block.size*nb, 0);   //initialize u

	    for (size_t q = block.q.start(), nqs = 0;
	    	 q < block.q.start()+block.q.size(); ++q) {

	    	// 1,2 index transform
	    	BOOST_AUTO(const &t2, transform(q, s)); //the trick is here to pick up on the parenthesis operator (eclipse (luna) misses it)
	    	size_t symm = 1 + (q != s);
	    	BOOST_AUTO(const &Q, shells.at(q));


	    	// pack t2
	    	size_t n = symm*Q.size()*S.size();
	    	for (size_t j = 0, ij = 0; j < no; ++j) {
	    	    for (size_t i = 0; i <= j; ++i, ++ij) {
	    		BOOST_AUTO(uij, u[ij/B]);
			size_t b = ij%B;  //b is the remainder when ij is divided by B
			size_t qs = nqs;
	    		for (size_t s = 0; s < S.size(); ++s) {
	    		    for (size_t q = 0; q < Q.size(); ++q) {
	    			uij[qs++][b] = t2[q][s][i][j];
	    		    }
	    		}
	    		if (symm == 1) continue;
	    		for (size_t q = 0; q < Q.size(); ++q) {
	    		    for (size_t s = 0; s < S.size(); ++s) {
	    			uij[qs++][b] = t2[q][s][j][i];
	    		    }
	    		}
	    	    }
		}
	    	nqs += n;

	    } // q in block

	    {
		BOOST_PROFILE_LINE;
		utility::timer timer;
		size_t start[] = { 0, block.index, 0 };
		size_t finish[] = { B, block.index+block.size, nb };

		V->put(u.data(), start, finish);
#pragma omp master
		{
		    profile.output += timer;
		    profile.size += block.size*B*nb;
		}
	    }

#pragma omp master
	    progress.jump(block.progress);

	} // blocks

#pragma omp master
	{
	    // += used in lieu of = causing ICE with some gcc 4.x versions
	    profile.eri += transform.profile.eri;  
	    profile.t1 += transform.profile.t1;
	    profile.t2 += transform.profile.t2;
	    std::copy(transform.profile.device,
		      transform.profile.device+3,
		      profile.device);
	    profile.transform += transform.profile.total;
	    progress.jump(blocks.back().progress);
	    used = memory.used();
	}
#pragma omp barrier

    } // parallel

    memory.clear();

    {
	BOOST_AUTO(cout, rt.cout());
	cout << std::endl;
	cout << "eri, transformations 1+2 (master thread):" << std::endl;
	cout << "    eri: " << profile.eri << std::endl;
	cout << "    trans 1: " << profile.t1 << std::endl;
	cout << "    trans 2: " << profile.t2 << std::endl;
	cout << "    I/O: "
	     << profile.output << ", "
	     << ((profile.size*sizeof(double)*1e-6)/
		 profile.output.total_seconds()) << " MB/s"
	     << std::endl;
	cout << "    time: " << profile.time << std::endl;
	cout << "    memory: " << used/(1<<20) << " MB" << std::endl;
    }

    {
	utility::timer timer;
	{
	    BOOST_PROFILE_LINE;
	    if (V) V->flush();
	    pe.barrier();
	}
	BOOST_AUTO(cout, rt.cout());
	cout << "synchronization: " << timer << std::endl;
	// cout << rt << std::endl;
	BOOST_PROFILE_DUMP(cout);
    }

    profile.time.reset();
    profile.size = 0;

    double E = 0;
	// e_terms stores energy components for open shell (note: for MP2 only e_term[0] is used)
	double e_term[7] = {0};

	//unite single+virtual orbitals 
	// (note: for closed shell, ns = 0 and nv is the number of virtuals)
	wf.set_single_and_virtuals(wf.double_occ().size(),wf.orbitals().size());
	size_t ns = wf.single().size();

    pe.task().reset();
    detail::Thread::Task<Parallel::Task&> task(pe.task());

    if (pe.node().rank() == 0) { // && pe.rank() == pe.size()-1) {

	BOOST_PROFILE_LINE;

	std::vector< std::pair<size_t,size_t> > index;
	for (size_t j = 0; j < no; ++j) {
	    for (size_t i = 0; i <= j; ++i) {
		index.push_back(std::pair<size_t,size_t>(i,j));
	    }
	}

	double *data1 = memory.malloc<double>(N*N*B); // get buffer
	double *data2 = memory.malloc<double>(N*N*B); // get buffer

	//Open Shell work (note: openshell object is instantiated regardless of wave function)
	cchem::mp2::detail::openshell OS_work(boost::cref(wf),memory, B);

	if(ns>0) OS_work.build_exchange_integrals(V,data2,index,memory); //exchange integrals are needed to modify singly occupied orbital energies

	task++;
	detail::async async;
	if (task < nb) async.get(0, B, N*N, task, *V, data2);
#pragma omp parallel reduction(+:E)
	{
	    detail::energy energy(boost::cref(wf), memory);		
	    double private_e_term[7] = {0}; //C++ OpenMP does not support array reduction so we have to do it ourselves
	    while (task < nb) {
		size_t b = task;

#pragma omp barrier
#pragma omp master
		{
		    BOOST_PROFILE_LINE;
		    utility::timer timer;
		    task++;
		    async.wait();
		    std::swap(data1, data2);
		    if (task < nb) {
			async.get(0, B, N*N, task, *V, data2);
		    }
		    profile.size += N*N*B;
		    profile.input += timer;
		}

#pragma omp barrier
#pragma omp for schedule(dynamic,1)
		for (int k = 0; k < (int)B; ++k) {
		    size_t ij = k + b*B;
		    if (ij >= index.size()) continue;
		    energy(index.at(ij), data1, B, k, OS_work.get_pfock() , private_e_term, OS_work.get_e_m());
		}
	    }

	    E = private_e_term[0]; //copy private closed shell mp2 energy to E

	    if(ns > 0)	OS_work(e_term, private_e_term, pe);// ZAPT energy terms only

	}  // omp parallel

	if(ns > 0) E = std::accumulate(e_term,e_term+7,0.0);//sum ZAPT energy components into E (for ROHF-PT, i.e. ns > 0)

    }
    used = memory.used();
    memory.clear();
    
    pe.reduce("+", &E, 1);

    if(ns > 0) for(int i = 0; i<7; i++)pe.reduce("+", &e_term[i], 1); //globally reduce ZAPT components

    {
	BOOST_AUTO(cout, rt.cout());
	cout << "transformations 3+4 (master thread):" << std::endl;
	cout << "    trans 3+4: " << profile.time << std::endl;
	cout << "    I/O: " << profile.input << ", "
	     << ((profile.size*sizeof(double)*1e-6)/
		 profile.input.total_seconds()) << " MB/s"
	     << std::endl;
	cout << "    memory: " << used/(1<<20) << " MB" << std::endl;
	cout << std::endl;
	BOOST_PROFILE_DUMP(cout);
    }

    //print out energy decomposition for ZAPT (i.e. ns > 0)
    if(ns > 0 && pe.rank()==0){
    std::cout << std::fixed <<std::setprecision(10) << std::showpos;
    std::cout <<" ZAPT ENERGY CONTRIBUTIONS" << std::endl;
    std::cout <<"       E1 = "<< e_term[0] <<"   CLOSED SHELL-LIKE TERM" << std::endl;
    std::cout <<"       E2 = "<< e_term[1] <<"   SINGLY UNOCCUPIED"<< std::endl;
    std::cout <<"       E3 = "<< e_term[2] <<"   2 SINGLY OCCUPIED"<< std::endl;
    std::cout <<"       E4 = "<< e_term[3] <<"   SINGLY UNOCCUPIED/OCCUPIED"<< std::endl;
    std::cout <<"       E5 = "<< e_term[4] <<"   \"FOCK MATRIX CONTRIBUTION\""<< std::endl;
    std::cout <<"       E6 = "<< e_term[5] <<"   SINGLY OCCUPIED"<< std::endl;
    std::cout <<"       E7 = "<< e_term[6] <<"   2 SINGLY UNOCCUPIED"<< std::endl;
    }
    rt.arrays().erase< Array<double> >("mp2.v(qsij)");

    return E;

}
} // namespace cchem
