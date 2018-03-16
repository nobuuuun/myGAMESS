/*
 * libcchem_libint.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: luke
 */



#include "utility/progress.hpp"
#include <libint2.hpp>


#include <vector>
//#include </opt/local/include/eigen3/Eigen/Core>
//#include </opt/local/include/eigen3/Eigen/Cholesky>
//#include </opt/local/include/eigen3/Eigen/Dense>
#include <Eigen/Dense>

#include "basis/basis.hpp"
#include "basis/molecule.hpp"

using libint2::Shell;
using libint2::Atom;
using libint2::BasisSet;

namespace libcchem_libint_interface{


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
Matrix;  // import dense, dynamically sized Matrix type from Eigen;
// this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
// to meet the layout of the integrals returned by the Libint integral library

std::vector<Atom> read_geometry(const std::string& filename);
Matrix compute_soad(const std::vector<Atom>& atoms);
// computes norm of shell-blocks of A
Matrix compute_shellblock_norm(const BasisSet& obs,
		const Matrix& A);
class GamessBasisSet;

template <libint2::OneBodyEngine::operator_type obtype>
std::array<Matrix, libint2::OneBodyEngine::operator_traits<obtype>::nopers>

compute_1body_ints(const GamessBasisSet& obs,
		const std::vector<Atom>& atoms = std::vector<Atom>());








template <libint2::OneBodyEngine::operator_type obtype>
std::vector<Matrix>
compute_1body_deriv_ints(unsigned deriv_order,
		const GamessBasisSet& obs,
		const std::vector<Atom>& atoms);

class GamessBasisSet : public std::vector<libint2::Shell> {


public:
	GamessBasisSet() : nbf_(-1), max_nprim_(0), max_l_(-1) {}
	GamessBasisSet(GamessBasisSet&& other) :
		std::vector<libint2::Shell>(std::move(other)),
		nbf_(other.nbf_),
		max_nprim_(other.max_nprim_),
		max_l_(other.max_l_),
		shell2bf_(std::move(other.shell2bf_)) {
	}
	~GamessBasisSet() = default;

	GamessBasisSet(::basis::Basis basis, Molecule molecule,const std::vector<Atom>& atoms){

		std::vector<libint2::Shell> ref_shells = GetBasisSetFromGamess(basis, molecule, atoms);

		//		int i = 0;
		for(auto s: ref_shells) {
			//			const Basis::Shell &Q = basis.shells().at(i);
			this->push_back(std::move(s));
			//			this->back().move({{atoms[Q.atom()].x, atoms[Q.atom()].y, atoms[Q.atom()].z}});
			//			i++;
		}//s

		init();

	}//GamessBasisSet

	std::vector<libint2::Shell> GetBasisSetFromGamess(::basis::Basis basis,
			Molecule molecule, const std::vector<Atom>& atoms){

		std::vector<libint2::Shell> new_ref;

		//BIG DEAL: inside of <libint install dir>/include/libint2/shell.h,
		// you must comment out renorm(); in order to compute the correct integrals.
		// you do not need to rebuild libint. The contraction coefficients are y normalized
		// and running those through the renormalization srews up. need to work this out later.
		for( auto ishell = 0; ishell < basis.shells().size(); ishell++){

			const Basis::Shell &Q = basis.shells().at(ishell);

			if(Q.size() != 4){

				auto l = (int)Q.L();
				std::vector<double> exps;
				std::vector<double> coeffs;

				for(int p = 0; p < Q.K(); p++){

					double e = Q.exp(p);
					double c = (Q.data()->C(0))[p];
					exps.emplace_back(e);
					coeffs.emplace_back(c);

				}//p

				new_ref.push_back(
						libint2::Shell{
					exps,
					{{l, 0, coeffs}},
					{{ atoms[Q.atom()].x, atoms[Q.atom()].y, atoms[Q.atom()].z }}
				}
				);

			}else{//if(Q.size() != 4)

				auto l = (int)Q.L();
				std::vector<double> exps;
				std::vector<double> coeffs_s, coeffs_p;

				for(int p = 0; p < Q.K(); p++){

					double e = Q.exp(p);
					double cs = (Q.data()->C(0))[p];
					double cp = (Q.data()->C(1))[p];

					exps.emplace_back(e);
					coeffs_s.emplace_back(cs);
					coeffs_p.emplace_back(cp);

				}//p

				new_ref.push_back(
						libint2::Shell{
					exps,
					{{0, 0, coeffs_s}},
					{{ atoms[Q.atom()].x, atoms[Q.atom()].y, atoms[Q.atom()].z }}
				}
				);

				new_ref.push_back(
						libint2::Shell{
					exps,
					{{1, 0, coeffs_p}},
					{{ atoms[Q.atom()].x, atoms[Q.atom()].y, atoms[Q.atom()].z }}
				}
				);

			}//if(Q.size() != 4)

		}//ishell

		return new_ref;
	}//GetBasisSetFromGamess


	/// @return the number of basis functions in the basis; -1 if uninitialized
	long nbf() const {
		return nbf_;
	}
	/// @return the maximum number of primitives in a contracted Shell, i.e. maximum contraction length; 0 if uninitialized
	size_t max_nprim() const {
		return max_nprim_;
	}
	/// @return the maximum angular momentum of a contraction; -1 if uninitialized
	size_t max_l() const {
		return max_l_;
	}
	/// @return the map from shell index to index of the first basis function from this shell
	/// \note basis functions are ordered as shells, i.e. shell2bf[i] >= shell2bf[j] iff i >= j
	const std::vector<size_t>& shell2bf() const {
		return shell2bf_;
	}

	/// Computes the map from shells to the corresponding atoms in \c atoms. If no atom matches the origin of a shell, it is mapped to -1.
	/// @note coordinates must match \em exactly , i.e. shell2atom[k] == l iff atoms[l].x == *this[k].O[0] && atoms[l].y == *this[k].O[1] &&  atoms[l].z == *this[k].O[2]
	/// @return the map from shell index to the atom in the list \c atoms that coincides with its origin;
	std::vector<long> shell2atom(const std::vector<Atom>& atoms) const {
		std::vector<long> result;
		result.reserve(size());
		for(const auto& s: *this) {
			auto a = std::find_if(atoms.begin(), atoms.end(), [&s](const Atom& a){ return s.O[0] == a.x && s.O[1] == a.y && s.O[2] == a.z; } );
			result.push_back( a != atoms.end() ? a - atoms.begin() : -1);
		}
		return result;
	}


private:
	long nbf_;
	size_t max_nprim_;
	int max_l_;
	std::vector<size_t> shell2bf_;

	void init() {
		nbf_ = nbf(*this);
		max_nprim_ = max_nprim(*this);
		max_l_ = max_l(*this);
		shell2bf_ = compute_shell2bf(*this);
	}


	static size_t nbf(const std::vector<libint2::Shell>& shells) {
		size_t n = 0;
		for (const auto& shell: shells)
			n += shell.size();
		return n;
	}

	static size_t max_nprim(const std::vector<libint2::Shell>& shells) {
		size_t n = 0;
		for (auto shell: shells)
			n = std::max(shell.nprim(), n);
		return n;
	}

	static int max_l(const std::vector<libint2::Shell>& shells) {
		int l = 0;
		for (auto shell: shells)
			for (auto c: shell.contr)
				l = std::max(c.l, l);
		return l;
	}

	static std::vector<size_t> compute_shell2bf(const std::vector<libint2::Shell>& shells) {
		std::vector<size_t> result;
		result.reserve(shells.size());

		size_t n = 0;
		for (auto shell: shells) {
			result.push_back(n);
			n += shell.size();
		}

		return result;
	}




};









void T_V_1body_deriv_contributions(size_t nbf, double * ptr_pmat, ::basis::Basis basis, Molecule molecule,
		double * ptr_wmat, Eigen::MatrixXd &TTerm, Eigen::MatrixXd &VTerm){

	//	int nthreads=1;




	std::vector<Atom> atoms;
	for (auto i = 0; i < molecule.size(); ++i){
		atoms.push_back( { molecule.get_atom(i).Z(),
			molecule(i)[0],
			molecule(i)[1],
			molecule(i)[2] } );
		//		std::cout << atoms[i].x  << " " <<
		//			atoms[i].y  << " " <<
		//			atoms[i].z  << " " << std::endl;

	}

	//    // set up thread pool
	//    {
	//      using libint2::nthreads;
	//      auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
	//      nthreads = 1;
	//      if (nthreads_cstr && strcmp(nthreads_cstr,"")) {
	//        std::istringstream iss(nthreads_cstr);
	//        iss >> nthreads;
	//        if (nthreads > 1<<16 || nthreads <= 0)
	//          nthreads = 1;
	//      }
	//#if defined(_OPENMP)
	//      omp_set_num_threads(nthreads);
	//#endif
	//      std::cout << "Will scale over " << nthreads
	//#if defined(_OPENMP)
	//                << " OpenMP"
	//#else
	//                << " C++11"
	//#endif
	//                << " threads" << std::endl;
	//    }

	GamessBasisSet obs2(basis,molecule,atoms);

	const auto n = obs2.nbf();
	const auto nshells = obs2.size();

	// count the number of electrons
	auto nelectron = 0;
	for (auto i = 0; i < atoms.size(); ++i)
		nelectron += atoms[i].atomic_number;
	const auto ndocc = nelectron / 2;

	// compute the nuclear repulsion energy
	auto enuc = 0.0;
	for (auto i = 0; i < atoms.size(); i++)
		for (auto j = i + 1; j < atoms.size(); j++) {
			auto xij = atoms[i].x - atoms[j].x;
			auto yij = atoms[i].y - atoms[j].y;
			auto zij = atoms[i].z - atoms[j].z;
			auto r2 = xij*xij + yij*yij + zij*zij;
			auto r = sqrt(r2);
			enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;

			//			std::cout << atoms[i].x << " " << atoms[j].x << std::endl;
			//			std::cout << atoms[i].y << " " << atoms[j].y	 << std::endl;
			//			std::cout << atoms[i].z << " " << atoms[j].z << std::endl;

		}

//	std::cout << "Nuclear repulsion energy = " << enuc << std::endl;

	libint2::init();

	//	auto T = compute_1body_ints<libint2::OneBodyEngine::kinetic>(obs2)[0];
	//    auto V = compute_1body_ints<libint2::OneBodyEngine::nuclear>(obs, atoms)[0];
	auto T1 = compute_1body_deriv_ints<libint2::OneBodyEngine::kinetic>(1, obs2, atoms);
	auto V1 = compute_1body_deriv_ints<libint2::OneBodyEngine::nuclear>(1, obs2, atoms);
//	auto O1 = compute_1body_deriv_ints<libint2::OneBodyEngine::overlap>(1, obs2, atoms);


	typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;

//	MapMatrixXd eg_global(ptr_eg_global,molecule.size(),3);

//	Eigen::MatrixXd eg(atoms.size(),3);
//	eg.setZero();
	MapMatrixXd pmat(ptr_pmat,n,n);
	MapMatrixXd wmat(ptr_wmat,n,n);

	for(auto atom=0, i=0; atom!=atoms.size(); ++atom) {
		for(auto xyz=0; xyz!=3; ++xyz, ++i) {
			//			auto force = (T1[i]).cwiseProduct(pmat).sum();
			auto force = T1[i].cwiseProduct(pmat).sum();
//			eg(atom, xyz) += force;
			TTerm(atom, xyz) += force;
		}//xyz
	}//atom
//	eg_global += eg;
//	std::cout << "kinetic energy gradient" << std::endl << eg << std::endl;

//	eg.setZero();

	for(auto atom=0, i=0; atom!=atoms.size(); ++atom) {
		for(auto xyz=0; xyz!=3; ++xyz, ++i) {
			auto force = V1[i].cwiseProduct(pmat).sum();
//			eg(atom, xyz) += force;
			VTerm(atom, xyz) += force;
		}//xyz
	}//atom
//	eg_global += eg;
//	std::cout << "nuclear (w/ Hellman) energy gradient" << std::endl << eg << std::endl;


//	eg.setZero();
//
//	for(auto atom=0, i=0; atom!=atoms.size(); ++atom) {
//		for(auto xyz=0; xyz!=3; ++xyz, ++i) {
//			auto force = O1[i].cwiseProduct(wmat).sum();
//			eg(atom, xyz) += force;
//		}//xyz
//	}//atom
//
//	std::cout << "overlap energy gradient" << std::endl << eg << std::endl;



	libint2::finalize();

	return;
}//T_V_1body_deriv_constributions















//template <libint2::OneBodyEngine::operator_type obtype>
//std::array<Matrix, libint2::OneBodyEngine::operator_traits<obtype>::nopers>
////compute_1body_ints(const BasisSet& obs,
////		const std::vector<Atom>& atoms)
//compute_1body_ints(const GamessBasisSet& obs,
//		const std::vector<Atom>& atoms)
//
//
//
//		{
//	const auto n = obs.nbf();
//	const auto nshells = obs.size();
////	std::cout << n <<  " " << nshells << std::endl;
//#ifdef _OPENMP
//	const auto nthreads = omp_get_max_threads();
//#else
//	const auto nthreads = 1;
//#endif
//	typedef std::array<Matrix, libint2::OneBodyEngine::operator_traits<obtype>::nopers> result_type;
//	const unsigned int nopers = libint2::OneBodyEngine::operator_traits<obtype>::nopers;
//	result_type result; for(auto& r: result) r = Matrix::Zero(n,n);
//
//	// construct the 1-body integrals engine
//	std::vector<libint2::OneBodyEngine> engines(nthreads);
//	engines[0] = libint2::OneBodyEngine(obtype, obs.max_nprim(), obs.max_l(), 0);
//
////	// construct the 2-electron repulsion integrals engine
////	typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
////	std::vector<coulomb_engine_type> tengines(nthreads);
////	tengines[0] = coulomb_engine_type(obs.max_nprim(), obs.max_l(), 0);
////	tengines[0].set_precision(0.); // !!! very important: cannot screen primitives in Schwartz computation !!!
////	for(size_t i=1; i!=nthreads; ++i) {
////		tengines[i] = tengines[0];
////	}
//
//
//	// nuclear attraction ints engine needs to know where the charges sit ...
//	// the nuclei are charges in this case; in QM/MM there will also be classical charges
//	if (obtype == libint2::OneBodyEngine::nuclear) {
//		std::vector<std::pair<double,std::array<double,3>>> q;
//		for(const auto& atom : atoms) {
//			q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
//		}
//		engines[0].set_params(q);
//	}
//	for(size_t i=1; i!=nthreads; ++i) {
//		engines[i] = engines[0];
//	}
//
//	auto shell2bf = obs.shell2bf();
//
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//	{
//#ifdef _OPENMP
//		auto thread_id = omp_get_thread_num();
//#else
//		auto thread_id = 0;
//#endif
////		std::cout << "nthreads "<< nthreads << std::endl;
//
//		// loop over unique shell pairs, {s1,s2} such that s1 >= s2
//		// this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
//		for(auto s1=0l, s12=0l; s1!=nshells; ++s1) {
//
//			auto bf1 = shell2bf[s1]; // first basis function in this shell
//			auto n1 = obs[s1].size();
//
//			for(auto s2=0; s2<=s1; ++s2) {
//
//				if (s12 % nthreads != thread_id)
//					continue;
//
//				auto bf2 = shell2bf[s2];
//				auto n2 = obs[s2].size();
//
//				auto n12 = n1 * n2;
//
//				// compute shell pair; return is the pointer to the buffer
//				const auto* buf = engines[thread_id].compute(obs[s1], obs[s2]);
////				std::cout << s1 << " " << s2 << std::endl;
////				const auto* tbuf = tengines[thread_id].compute(obs[s1], obs[s2],obs[s1], obs[s2]);
//				//        std::cout << n12 << std::endl;
////				        std::cout << buf[0] << std::endl;
//				for(unsigned int op=0; op!=nopers; ++op, buf+=n12) {
//
//					// "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
//					          Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
//					          result[op].block(bf1, bf2, n1, n2) = buf_mat;
////					          std::cout << buf_mat << std::endl;
//					          if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
//					            result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
//				}
//
//			}
//		}
//
//
//
//		//    std::cout << result[0] << std::endl;
//	} // omp parallel
//
//	return result;
//		}//compute_1body_ints












template <libint2::OneBodyEngine::operator_type obtype>
std::vector<Matrix>
compute_1body_deriv_ints(unsigned deriv_order,
		const GamessBasisSet& obs,
		const std::vector<Atom>& atoms)
		{

	const auto n = obs.nbf();
	const auto nshells = obs.size();
#ifdef _OPENMP
	const auto nthreads = omp_get_max_threads();
#else
	const auto nthreads = 1;
#endif
	//	  std::cout << libint2::OneBodyEngine::operator_traits<obtype>::nopers << std::endl;
	constexpr auto nopers = libint2::OneBodyEngine::operator_traits<obtype>::nopers;
	//	std::cout << nopers << std::endl;
	const auto nresults = nopers * libint2::num_geometrical_derivatives(atoms.size(),deriv_order);
	typedef std::vector<Matrix> result_type;
	result_type result(nresults); for(auto& r: result) r = Matrix::Zero(n,n);
	//	  std::cout << nopers << std::endl;

	// construct the 1-body integrals engine
	std::vector<libint2::OneBodyEngine> engines(nthreads);
	engines[0] = libint2::OneBodyEngine(obtype, obs.max_nprim(), obs.max_l(), deriv_order);
	// nuclear attraction ints engine needs to know where the charges sit ...
	// the nuclei are charges in this case; in QM/MM there will also be classical charges
	if (obtype == libint2::OneBodyEngine::nuclear) {
		std::vector<std::pair<double,std::array<double,3>>> q;
		for(const auto& atom : atoms) {
			q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
		}
		engines[0].set_params(q);
	}
	for(size_t i=1; i!=nthreads; ++i) {
		engines[i] = engines[0];
	}

	auto shell2bf = obs.shell2bf();
	auto shell2atom = obs.shell2atom(atoms);


#ifdef _OPENMP
#pragma omp parallel
#endif
	{
#ifdef _OPENMP
		auto thread_id = omp_get_thread_num();
#else
		auto thread_id = 0;
#endif
		//		std::cout << thread_id << std::endl;

		// loop over unique shell pairs, {s1,s2} such that s1 >= s2
		// this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
		for(auto s1=0l, s12=0l; s1!=nshells; ++s1) {

			auto bf1 = shell2bf[s1]; // first basis function in this shell
			auto n1 = obs[s1].size();
			auto atom1 = shell2atom[s1];
			assert(atom1 != -1);

			for(auto s2=0; s2<=s1; ++s2) {

				if (s12 % nthreads != thread_id)
					continue;

				auto bf2 = shell2bf[s2];
				auto n2 = obs[s2].size();
				auto atom2 = shell2atom[s2];
				auto n12 = n1 * n2;


				// compute shell pair; return is the pointer to the buffer
				const auto* buf = engines[thread_id].compute(obs[s1], obs[s2]);

				assert(deriv_order == 1); // the loop structure below needs to be generalized for higher-order derivatives
				// 1. process derivatives with respect to the Gaussian origins first ...
				for(unsigned int d=0; d!=6; ++d) { // 2 centers x 3 axes = 6 cartesian geometric derivatives
					auto atom = d < 3 ? atom1 : atom2;
					auto op_start = (3*atom+d%3) * nopers;
					auto op_fence = op_start + nopers;
					for(unsigned int op=op_start; op!=op_fence; ++op, buf+=n12) {
						Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
						result[op].block(bf1, bf2, n1, n2) += buf_mat;
						if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
							result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
					}//op
				}//d

				// 2. process derivatives of nuclear Coulomb operator, if needed
				if (obtype == libint2::OneBodyEngine::nuclear) {
					for(unsigned int atom=0; atom!=atoms.size(); ++atom) {
						for(unsigned int xyz=0; xyz!=3; ++xyz) {
							auto op_start = (3*atom+xyz) * nopers;
							auto op_fence = op_start + nopers;
							for(unsigned int op=op_start; op!=op_fence; ++op, buf+=n12) {
								Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
								result[op].block(bf1, bf2, n1, n2) += buf_mat;
								if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
									result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
							}//op
						}//xyz
					}//atom
				}//(obtype == libint2::OneBodyEngine::nuclear)

			}//s2
		}//s1
	} // omp parallel

	return result;
		}//compute_1body_deriv_ints
}//namespace libcchem_libint_interface



































//#include <iostream>
//#include <cmath>
//#include <sys/time.h>
//#include <cassert>
//
//#include <libint2.h>
//#include <libint2.hpp>
//#include <prep_libint2.h>
//#include <libint2/cgshell_ordering.h>
//
//libint2::FmEval_Chebyshev3<double> fmeval_chebyshev(std::max(LIBINT_MAX_AM,4)*4 + 2);
//libint2::FmEval_Taylor<double,6> fmeval_taylor(std::max(LIBINT_MAX_AM,4)*4 + 2, 1e-15);
//
//namespace libcchem_libint_interface{
//void test_4eri(int deriv_order,
//		int lmax_max) {
//
//	LIBINT2_PREFIXED_NAME(libint2_static_init)();
//
//	typedef unsigned int uint;
//
//	size_t veclen = LIBINT2_MAX_VECLEN;
//	size_t max_contrdepth = 3;
//	size_t max_contrdepth4 = max_contrdepth * max_contrdepth * max_contrdepth * max_contrdepth;
//	size_t contrdepth = max_contrdepth;
//
//	int lmax;
//	if (deriv_order == 0) lmax = LIBINT2_MAX_AM_ERI;
//	if (deriv_order == 1) lmax = LIBINT2_MAX_AM_ERI1;
//
//	Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth4);
//
//	//	void * aLibint_t = NULL;
//	//    std::cout << aLibint_t << std::endl;
//	//	posix_memalign(&aLibint_t, 8, max_contrdepth4*sizeof(Libint_t) );
//	//    std::cout << aLibint_t << std::endl;
//	//	Libint_t *inteval = new(aLibint_t) Libint_t[max_contrdepth4];
//
//
//
//	size_t ndubs;
//	if (deriv_order == 0) ndubs = libint2_need_memory_eri(lmax_max);
//	if (deriv_order == 1) ndubs = libint2_need_memory_eri1(lmax_max);
//	void * adubs = NULL;
//	posix_memalign(&adubs, 16, ndubs*sizeof(double) );
//	void *ptr_ndubs = new(adubs) double[ndubs];
//
//
//	//	std::cout << inteval << std::endl;
//	if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], // ptr to the first evaluator in the array
//			lmax,        // maximum angular momentum
//			ptr_ndubs);
//
//	if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_init_eri1)(&inteval[0], // ptr to the first evaluator in the array
//			lmax,        // maximum angular momentum
//			ptr_ndubs);
//
//
//
//
//	lmax = std::min(lmax_max, lmax);
//
//	//	for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
//	//		for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
//	//			for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
//	//				for (unsigned int l3 = 0; l3 <= lmax; ++l3) {
//
//	for (size_t l0 = 0; l0 <= lmax; ++l0) {
//		for (size_t l1 = 0; l1 <= lmax; ++l1) {
//			for (size_t l2 = 0; l2 <= lmax; ++l2) {
//				for (size_t l3 = 0; l3 <= lmax; ++l3) {
//
//					if(l0 != 0  ||  l1 != 0  ||  l2 != 0  ||  l3 != 0  )continue;
//
//					std::cout << l0 << " "<< l1 << " "<< l2 << " "<< l3 <<std::endl;
//					//	        	  for (unsigned int l0 = 0; l0 <= 0; ++l0) {
//					//	        	    for (unsigned int l1 = 0; l1 <= 0; ++l1) {
//					//	        	      for (unsigned int l2 = 0; l2 <= 0; ++l2) {
//					//	        	        for (unsigned int l3 = 0; l3 <= 0; ++l3) {
//
//					if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0){
//						std::cout << "   skipping" << std::endl;
//						continue;
//					}
//
//					if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_eri1)[l0][l1][l2][l3] == 0)
//						continue;
//
//
//					std::cout << " NOT skipping" << std::endl;
//					unsigned int am[4];
//					am[0] = l0;
//					am[1] = l1;
//					am[2] = l2;
//					am[3] = l3;
//					RandomShellSet<4u> rsqset(am, veclen, contrdepth);
//
//					//exp[4:shell][1:vectorization][contraction]
//
//					//suppose we have 4 shells, vectorization is 1, and contraction is 6
//					//for shell i    exp[i] : c0 c1 c2 c3 c4 c5
//
//					//suppose we have 4 shells, vectorization is 2, and contraction is 6
//					//for shell i    exp[i] : c0 c1
//					//                        c2 c3
//					//                        c4 c5
//
//
//					for (int i = 0; i < 4; i++ ){
//						std::cout << rsqset.R[i][0] <<  " "
//								<< rsqset.R[i][1] <<  " "
//								<< rsqset.R[i][2] <<  std::endl;
//						std::cout << "   " << rsqset.exp[i][0][0] << " " << rsqset.coef[i][0][0]<< std::endl << std::endl;
//					}//i
//
//
//					prep_libint2(inteval, rsqset, 0, deriv_order);
//
//					double scale_target = 1.0;
//
//					inteval[0].contrdepth = max_contrdepth4;
//					if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
//					if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
//
//
//					std::vector<LIBINT2_REALTYPE> new_eri;
//					if(deriv_order == 1){
//						for (unsigned int di = 0; di < 12; ++di) {
//							LIBINT2_REALTYPE value;
//							if (di<6)
//								value = scale_target * inteval[0].targets[di][0];
//							else if (di>=9)
//								value = scale_target * inteval[0].targets[di-3][0];
//							else // (di>=6 || di<=8)
//								value = -scale_target * (inteval[0].targets[di-6][0] +
//										inteval[0].targets[di-3][0] +
//										inteval[0].targets[di][0]);
//							new_eri.push_back(value);
//						}//di
//
//						for(unsigned int di = 0; di < new_eri.size(); di++){
//							std::cout << "look at this:  "<< new_eri[di] << std::endl;
//						}//di
//					}//(deriv_order = 1
//
//				}//l3
//			}//l2
//		}//l1
//	}//l0
//
//
//
//
//
//	//	LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
//	//	free(inteval);
//	delete inteval;
//	LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();
//
//
//}
//
//
//
//
//
//
//void test_T(int deriv_order,
//		int lmax_max) {
//
//	LIBINT2_PREFIXED_NAME(libint2_static_init)();
//
//	typedef unsigned int uint;
//
//	size_t veclen = LIBINT2_MAX_VECLEN;
//	size_t max_contrdepth = 1;
//	size_t max_contrdepth2 = max_contrdepth * max_contrdepth;
//	size_t contrdepth = max_contrdepth;
//
//	int lmax;
//	if (deriv_order == 0) lmax = LIBINT2_MAX_AM_kinetic;
//	if (deriv_order == 1) lmax = LIBINT2_MAX_AM_kinetic1;
//
//	Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth2);
//
//	//	void * aLibint_t = NULL;
//	//    std::cout << aLibint_t << std::endl;
//	//	posix_memalign(&aLibint_t, 8, max_contrdepth4*sizeof(Libint_t) );
//	//    std::cout << aLibint_t << std::endl;
//	//	Libint_t *inteval = new(aLibint_t) Libint_t[max_contrdepth4];
//
//
//
//	size_t ndubs;
//	if (deriv_order == 0) ndubs = libint2_need_memory_kinetic(lmax_max);
//	if (deriv_order == 1) ndubs = libint2_need_memory_kinetic1(lmax_max);
//	void * adubs = NULL;
//	posix_memalign(&adubs, 16, ndubs*sizeof(double) );
//	void *ptr_ndubs = new(adubs) double[ndubs];
//
//
//	//	std::cout << inteval << std::endl;
//	if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_init_kinetic)(&inteval[0], // ptr to the first evaluator in the array
//			lmax,        // maximum angular momentum
//			ptr_ndubs);
//
//	if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_init_kinetic1)(&inteval[0], // ptr to the first evaluator in the array
//			lmax,        // maximum angular momentum
//			ptr_ndubs);
//
//
//
//
//	lmax = std::min(lmax_max, lmax);
//
//	//	for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
//	//		for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
//	//			for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
//	//				for (unsigned int l3 = 0; l3 <= lmax; ++l3) {
//
//	for (size_t l0 = 0; l0 <= lmax; ++l0) {
//		for (size_t l1 = 0; l1 <= lmax; ++l1) {
//			//			for (size_t l2 = 0; l2 <= lmax; ++l2) {
//			//				for (size_t l3 = 0; l3 <= lmax; ++l3) {
//
//			//					if(l0 != 0  ||  l1 != 0  )continue;
//
//			std::cout << l0 << " "<< l1 <<std::endl;
//			//	        	  for (unsigned int l0 = 0; l0 <= 0; ++l0) {
//			//	        	    for (unsigned int l1 = 0; l1 <= 0; ++l1) {
//			//	        	      for (unsigned int l2 = 0; l2 <= 0; ++l2) {
//			//	        	        for (unsigned int l3 = 0; l3 <= 0; ++l3) {
//
//			if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_kinetic)[l0][l1] == 0){
//				std::cout << "   skipping" << std::endl;
//				continue;
//			}
//
//			if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_kinetic1)[l0][l1] == 0)
//				continue;
//
//
//			std::cout << " NOT skipping" << std::endl;
//			unsigned int am[2];
//			am[0] = l0;
//			am[1] = l1;
//			//					am[2] = l2;
//			//					am[3] = l3;
//			RandomShellSet<2u> rsqset(am, veclen, contrdepth);
//
//			//exp[4:shell][1:vectorization][contraction]
//
//			//suppose we have 4 shells, vectorization is 1, and contraction is 6
//			//for shell i    exp[i] : c0 c1 c2 c3 c4 c5
//
//			//suppose we have 4 shells, vectorization is 2, and contraction is 6
//			//for shell i    exp[i] : c0 c1
//			//                        c2 c3
//			//                        c4 c5
//
//
//			for (int i = 0; i < 2; i++ ){
//				std::cout << rsqset.R[i][0] <<  " "
//						<< rsqset.R[i][1] <<  " "
//						<< rsqset.R[i][2] <<  std::endl;
//				std::cout << "   " << rsqset.exp[i][0][0] << " " << rsqset.coef[i][0][0]<< std::endl << std::endl;
//			}//i
//
//
//			prep_libint2(inteval, rsqset, 0, deriv_order);
//
//			double scale_target = 1.0;
//
//			inteval[0].contrdepth = max_contrdepth2;
//			if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_build_kinetic)[am[0]][am[1]](&inteval[0]);
//			if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_build_kinetic1)[am[0]][am[1]](&inteval[0]);
//
//
//			std::vector<LIBINT2_REALTYPE> new_eri;
//			if(deriv_order == 1){
//				for (unsigned int di = 0; di < 3; ++di) {
//					LIBINT2_REALTYPE value;
//					//							if (di<3)
//					value = scale_target * inteval[0].targets[di][0];
//					std::cout << value << std::endl;
//					//							else if (di>=9)
//					//								value = scale_target * inteval[0].targets[di-3][0];
//					//							else // (di>=6 || di<=8)
//					//								value = -scale_target * (inteval[0].targets[di-6][0] +
//					//										inteval[0].targets[di-3][0] +
//					//										inteval[0].targets[di][0]);
//					new_eri.push_back(value);
//				}//di
//
//				for(unsigned int di = 0; di < new_eri.size(); di++){
//					std::cout << "look at this:  "<< new_eri[di] << std::endl;
//				}//di
//			}//(deriv_order = 1
//
//
//
//			if(deriv_order == 0){
//				int ij = 0;
//				int l00 = l0+1;
//				int l10 = l1+1;
//				for (int i = 0; i < l00*(l00+1)/2; i++){
//					for (int j = 0; j < l10*(l10+1)/2; j++,ij++){
//						LIBINT2_REALTYPE value;
//						//							if (di<3)
//						value = scale_target * inteval[0].targets[0][ij];
//						std::cout << value << std::endl;
//						//							else if (di>=9)
//						//								value = scale_target * inteval[0].targets[di-3][0];
//						//							else // (di>=6 || di<=8)
//						//								value = -scale_target * (inteval[0].targets[di-6][0] +
//						//										inteval[0].targets[di-3][0] +
//						//										inteval[0].targets[di][0]);
//						//									new_eri.push_back(value);
//
//
//						//								for(unsigned int di = 0; di < new_eri.size(); di++){
//						//									std::cout << "look at this:  "<< new_eri[di] << std::endl;
//						//								}//di
//
//					}
//				}
//
//			}//(deriv_order = 1
//
//
//
//			//				}//l3
//			//			}//l2
//		}//l1
//	}//l0
//
//
//
//
//
//	//	LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
//	//	free(inteval);
//	delete inteval;
//	LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();
//
//
//}
//
//
//
//
//
//}//namespace libcchem_libint_interface


//#include <boost/noncopyable.hpp>
#include <boost/typeof/typeof.hpp>
//#include <boost/thread.hpp>

typedef unsigned int uint;
template <unsigned int N>
struct ShellSet {
private:
	boost::reference_wrapper<const Basis> basis_;

public:
	ShellSet(boost::reference_wrapper<const Basis> basis) :
		basis_(basis){};



	void add_shell(size_t veclen,  //for vectorization
			uint shell_num, //which shell : 1->num_shells
			uint qpos,
			uint start) {  //0,1,2,3

		BOOST_AUTO(const &basis, basis_.get());

		const Basis::Shell &Shell = basis.shells().at(shell_num);

		l[qpos] = (uint)Shell.L();
		R[qpos].resize(3);
		R[qpos][0] = Shell.center(0);
		R[qpos][1] = Shell.center(1);
		R[qpos][2] = Shell.center(2);

		exp[qpos].resize(1);
		exp[qpos][0].resize(Shell.K());
		//		exp[qpos][0].resize(1);

		coef[qpos].resize(1);
		coef[qpos][0].resize(Shell.K());
		//		coef[qpos][0].resize(1);

		for (uint i = 0; i < Shell.K(); i++){
			//			for (uint i = 0; i < 1; i++){
			exp[qpos][0][i] = Shell.exp(i);
			//				exp[qpos][0][i] = Shell.exp(start);
			coef[qpos][0][i] = Shell.data()->C(0)[i];
			//			coef[qpos][0][i] = Shell.data()->C(0)[start];
		}

		//		coef[qpos][0][0] = 0.033498726390;
		//		coef[qpos][0][1] = 0.234800801174;
		//		coef[qpos][0][2] = 0.813682957883;

	}

	void print(){

		for (int qpos = 0; qpos < 4 ; qpos++){


			std::cout << l[qpos] << std::endl;
			std::cout << R[qpos][0] << " " << R[qpos][1] << " "  << R[qpos][2] << std::endl;

			for (int i = 0; i < exp[qpos][0].size(); i++){
				std::cout << exp[qpos][0][i] << " " << coef[qpos][0][i] << std::endl;
			}

			std::cout<<std::endl;
		}//i

	}//print

	uint l[N];                                  // angular momenta
	std::vector<double> R[N];                   // origins
	std::vector< std::vector<double> > exp[N];  // exponents
	std::vector< std::vector<double> > coef[N]; // coefficients
};//struct ShellSet



//#include <iostream>
//#include <cmath>
//#include <sys/time.h>
//#include <cassert>
//
//#include <libint2.h>
//#include <libint2.hpp>
//#include <prep_libint2.h>
//#include <libint2/cgshell_ordering.h>

libint2::FmEval_Chebyshev3<double> fmeval_chebyshev(std::max(LIBINT_MAX_AM,4)*4 + 2);
libint2::FmEval_Taylor<double,6> fmeval_taylor(std::max(LIBINT_MAX_AM,4)*4 + 2, 1e-15);



namespace libcchem_libint_interface
{


#include <prep_libint2.h>


void compute_schwartz_ints(Basis basis, Eigen::MatrixXd &schwartz_mat){

	LIBINT2_PREFIXED_NAME(libint2_static_init)();
	typedef unsigned int uint;
	size_t veclen = LIBINT2_MAX_VECLEN;
	//figure out the maximum contraction depth
	//figure out the maximum shell size
	size_t max_contrdepth = 0;
	size_t max_shell_size = 0;
	for(int ishell = 0; ishell < basis.shells().size(); ishell++){
		const Basis::Shell &Q = basis.shells().at(ishell);
		if(Q.K() > max_contrdepth) max_contrdepth = Q.K();
		if(Q.size() > max_shell_size) max_shell_size = Q.size();
//		std::cout << basis.max().L()<< " "<<max_contrdepth << std::endl;
	}//ishell
	size_t max_contrdepth4 = max_contrdepth*max_contrdepth*max_contrdepth*max_contrdepth;
	int libint_lmax = LIBINT2_MAX_AM_ERI1;
	int basis_lmax = basis.max().L();
	if(basis_lmax > libint_lmax) {
		std::cout << "this is not good; libint is not set up for this high of momentum" <<std::endl
		<< "max libint momentum l=" << libint_lmax <<std::endl
		<< "but you have functions with moment l=" << basis_lmax <<std::endl
		<< "rebuild libint with higher derivative momentum or" <<std::endl
		<< "select a different basis set" <<std::endl;
	}

	Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth4);

	LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], // ptr to the first evaluator in the array
			basis_lmax,        // maximum angular momentum
			0);


	uint am[4];

	ShellSet<4u> qset(boost::cref(basis) );
	//	(p q| r s)
	for(uint p = 0; p < basis.shells().size(); p++){
		const Basis::Shell &P1 = basis.shells().at(p);
		int pam = P1.L();
		int ram = pam;
		uint r = p;

		for(uint q = 0; q <= p; q++){
			const Basis::Shell &Q1 = basis.shells().at(q);
			int qam = Q1.L();
			int sam = qam;
			uint s = q;

			int pnew = p;
			int qnew = q;
			int rnew = r;
			int snew = s;

			//new to order shells special for libint
			if(pam < qam)std::swap(pnew,qnew);
			if(ram < sam)std::swap(snew,rnew);
			if( (pam+qam) > (ram+sam) ){
				std::swap(pnew,rnew);
				std::swap(qnew,snew);
			}



			qset.add_shell(veclen,pnew,0,0);
			qset.add_shell(veclen,qnew,1,0);
			qset.add_shell(veclen,rnew,2,0);
			qset.add_shell(veclen,snew,3,0);

			const Basis::Shell &P = basis.shells().at(pnew);
			am[0] = P.L();
			const Basis::Shell &Q = basis.shells().at(qnew);
			am[1] = Q.L();
			const Basis::Shell &R = basis.shells().at(rnew);
			am[2] = R.L();
			const Basis::Shell &S = basis.shells().at(snew);
			am[3] = S.L();

			inteval[0].contrdepth = P.K()*Q.K()*R.K()*S.K();

			if (LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]] == 0)continue;
			//			qset.print();

			prep_libint2(inteval, qset, 0, 1);

			// (p q | r s)
			// P.L() >= Q.L(), R.L() >= S.L(), R.L()+S.L() >= P.L()+Q.L()
			LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);

			typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
			MapMatrixXd buf_mat(&inteval[0].targets[0][0], P.size()*Q.size(), P.size()*Q.size() );


			double norm = buf_mat.diagonal().lpNorm<Eigen::Infinity>();
			//			schwartz_mat(p,q) = buf_mat.diagonal().norm();
			schwartz_mat(p,q) = std::sqrt(norm);
			schwartz_mat(q,p) = schwartz_mat(pnew,qnew);
		}//q

	}//p

	LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);

	free(inteval);
	LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();


}//compute_schwartz_ints



//void twobody_deriv_contributions(size_t nbf, double * ptr_pmp2, double * ptr_pscf,
//		Basis basis,
//		Molecule molecule,double *ptr_eg_global,
//		Eigen::MatrixXd &schwartz_mat){
//
//		typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
//
//		MapMatrixXd eg_global(ptr_eg_global,molecule.size(),3);
//		MapMatrixXd pmp2(ptr_pmp2,nbf,nbf);
//		MapMatrixXd pscf(ptr_pscf,nbf,nbf);
//
//
//		LIBINT2_PREFIXED_NAME(libint2_static_init)();
//
//		typedef unsigned int uint;
//
//		size_t veclen = LIBINT2_MAX_VECLEN;
//
//		//figure out the maximum contraction depth
//		//figure out the maximum shell size
//		size_t max_contrdepth = 0;
//		size_t max_shell_size = 0;
//		for(int ishell = 0; ishell < basis.shells().size(); ishell++){
//			const Basis::Shell &Q = basis.shells().at(ishell);
//			if(Q.K() > max_contrdepth) max_contrdepth = Q.K();
//			if(Q.size() > max_shell_size) max_shell_size = Q.size();
//			std::cout << basis.max().L()<< " "<<max_contrdepth << std::endl;
//		}//ishell
//		size_t max_contrdepth4 = max_contrdepth*max_contrdepth*max_contrdepth*max_contrdepth;
////		size_t contrdepth = max_contrdepth;
//
////		//enough space for the larget p,q,r,s shell combination of density
////		size_t ndubs = max_shell_size*max_shell_size*max_shell_size*max_shell_size;
////		void * adubs = NULL;
////		posix_memalign(&adubs, 16, ndubs*sizeof(double) );
////		double *ptr_den = new(adubs) double[ndubs];
//
//
//
//		int libint_lmax = LIBINT2_MAX_AM_ERI1;
//		int basis_lmax = basis.max().L();
//		if(basis_lmax > libint_lmax) {
//			std::cout << "this is not good; libint is not set up for this high of momentum" <<std::endl;
//			std::cout << "max libint momentum l=" << libint_lmax <<std::endl;
//			std::cout << "but you have functions with moment l=" << basis_lmax <<std::endl;
//			std::cout << "rebuild libint with higher derivative momentum or" <<std::endl;
//			std::cout << "select a different basis set" <<std::endl;
//		}
//
//		uint global_skipped = 0;
//#pragma omp parallel
//		{
//
//			Eigen::MatrixXd eg(molecule.size(),3);
//			eg.setZero();
////#pragma omp critical
////			std::cout << eg << std::endl;
//
////#pragma omp barrier
//			Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth4);
//
////		size_t ndubs;
////		ndubs = libint2_need_memory_eri1(basis_lmax);
////		void * adubs = NULL;
////		posix_memalign(&adubs, 16, ndubs*sizeof(double) );
////		double *ptr_ndubs = new(adubs) double[ndubs];
//
//		LIBINT2_PREFIXED_NAME(libint2_init_eri1)(&inteval[0], // ptr to the first evaluator in the array
//				basis_lmax,        // maximum angular momentum
//				0);
//
//		//        exp[c].resize(veclen);
//		//        coef[c].resize(veclen);
//		//        for(uint v=0; v<veclen; ++v) {
//		//          exp[c][v].resize(contrdepth); generate(exp[c][v].begin(), exp[c][v].end(), die);
//		//          coef[c][v].resize(contrdepth); generate(coef[c][v].begin(), coef[c][v].end(), die);
//		//        }
//
//		uint am[4];
////	    std::vector<double> Center[4];                   // origins
////	    std::vector< std::vector<double> > exp[4];  // exponents
////	    std::vector< std::vector<double> > coef[4]; // coefficients
////	    for(uint i = 0; i < 4; i++)Center[i].resize(3);
////		for(uint i = 0; i < 4; i++)
////			{exp[i].resize(1); //right now veclen is hardwired to 1
////			exp[i][0].resize(max_contrdepth);
////			}//i
////		for(uint i = 0; i < 4; i++)
////			{coef[i].resize(1); //right now veclen is hardwired to 1
////			coef[i][0].resize(max_contrdepth);
////			}//i
//
//		ShellSet<4u> qset(boost::cref(basis) );
//		uint skipped = 0;
////		(p q | r s)
//#pragma omp for schedule(dynamic)
//		for(uint p = 0; p < basis.shells().size(); p++){
//			const Basis::Shell &P1 = basis.shells().at(p);
//			am[0] = P1.L();
//			int pam = P1.L();
//			std::cout << omp_get_thread_num() << " " << p << std::endl;
//
////			for(uint q = 0; q < basis.shells().size(); q++){
//				for(uint q = 0; q <= p; q++){
//				const Basis::Shell &Q1 = basis.shells().at(q);
//				am[1] = Q1.L();
//				int qam = Q1.L();
//
////				if(omp_get_thread_num() == 0)std::cout << schwartz_mat(p,q)*schwartz_mat(p,q)<< std::endl;
////				for(uint r = 0; r < basis.shells().size(); r++){  //noddy
////				for(uint r = 0; r <= p; r++){   //not noddy
//					for(uint r = 0; r <= p; r++){   //not noddy
//					const Basis::Shell &R1 = basis.shells().at(r);
//					am[2] = R1.L();
//					int ram = R1.L();
//
////					for(uint s = 0; s < basis.shells().size(); s++){
//			         const auto s4_max = (p == r) ? q : r;
//						for(uint s = 0; s <= s4_max; s++){
//						const Basis::Shell &S1 = basis.shells().at(s);
//						am[3] = S1.L();
//						int sam = S1.L();
//
//
//						bool debug = 0;
//
//						if(P1.atom() ==  Q1.atom() && R1.atom() == S1.atom() && P1.atom() == R1.atom()) continue;
//
//						if( schwartz_mat(p,q)*schwartz_mat(r,s) < (double)1.0E-9 ){
//							skipped++;
//							continue;
//						}
//
//
////						if(schwartz_mat(p,q)*schwartz_mat(r,s) < 0.000000000001 )std::cout << "skip" << std::endl;
//
//
////			            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
////			            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
////			            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
////			            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;
//
//			            auto s12_deg = (p == q) ? 1.0 : 2.0;
//			            auto s34_deg = (r == s) ? 1.0 : 2.0;
//			            auto s12_34_deg = (p == r) ? (q == s ? 1.0 : 2.0) : 2.0;
//			            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;
//
//
//
//						int pnew = p;
//						int qnew = q;
//						int rnew = r;
//						int snew = s;
//
//						if(pam < qam){
//							pnew = q; qnew = p;
//						}
//						if(ram < sam){
//							rnew = s; snew = r;
//						}
//						if( (pam+qam) > (ram+sam) ){
//							int temp = rnew;
//							rnew = pnew;
//							pnew = temp;
//							temp = snew;
//							snew = qnew;
//							qnew = temp;
//						}
//
//
//
//
////						double sfac =1.0;
//////						if(p != r || q != s)sfac /= 2.0;
////						if(p != q)sfac /= 2.0;
////						if(r != s)sfac /= 2.0;
////						if(p != r)sfac *= 2.0;
////						if(std::max(p,q) != std::max(r,s) ||
////								std::min(p,q) != std::min(r,s) ) sfac /=2;
//
////										             IF(IIJJ.NE.KKLL) SFAC = HALF
//						qset.add_shell(veclen,pnew,0,0);
//						qset.add_shell(veclen,qnew,1,0);
//						qset.add_shell(veclen,rnew,2,0);
//						qset.add_shell(veclen,snew,3,0);
//
//						const Basis::Shell &P = basis.shells().at(pnew);
//						am[0] = P.L();
//						const Basis::Shell &Q = basis.shells().at(qnew);
//						am[1] = Q.L();
//						const Basis::Shell &R = basis.shells().at(rnew);
//						am[2] = R.L();
//						const Basis::Shell &S = basis.shells().at(snew);
//						am[3] = S.L();
//
//
//						inteval[0].contrdepth = P.K()*Q.K()*R.K()*S.K();
//
//
//			if (LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]] == 0)continue;
////			qset.print();
//
////			continue;
//			prep_libint2(inteval, qset, 0, 1);
//
////			continue;
//
//			// (p q | r s)
//			// P.L() >= Q.L(), R.L() >= S.L(), R.L()+S.L() >= P.L()+Q.L()
//			LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);
//
//
//
//
//
//
//			for (int pf = 0; pf < P.size(); pf++){
//
//				for (int qf = 0; qf < Q.size(); qf++){ //noddy
////				int maxqf =Q.size()-1; //-1 for the '<=' below
////				if (pnew == qnew) maxqf = pf;
////					for (int qf = 0; qf <= maxqf; qf++){
//
//
//
//						for (int rf = 0; rf < R.size(); rf++){ //noddy
////						int maxrf = R.size()-1; //-1 for the '<=' below
////						if(pnew == rnew && qnew == snew)maxrf = pf;
////						for (int rf = 0; rf <=maxrf; rf++){
//
//
//						for (int sf = 0; sf < S.size(); sf++){ //noddy
////						int maxsf =S.size()-1; //-1 for the '<=' below
////						if (rnew == snew) maxsf = rf;
////						for (int sf = 0; sf <= maxsf; sf++){
//
//							int mu = P.start() +pf;
//							int nu = Q.start() +qf;
//							int la = R.start() +rf;
//							int si = S.start() +sf;
//
//							double dhf=0.0;
//							dhf = pscf(mu,nu)*pscf(la,si)*4 -
//									pscf(mu,si)*pscf(la,nu) -
//									pscf(mu,la)*pscf(nu,si);
////							if(debug)std::cout << "dhf       "<< dhf << std::endl;
//
//
//							double dsep=0.0;
//							dsep = pmp2(mu,nu)*pscf(la,si)*4 +
//									pmp2(la,si)*pscf(mu,nu)*4 -
//									pmp2(mu,si)*pscf(la,nu) -
//									pmp2(mu,la)*pscf(nu,si) -
//									pmp2(nu,si)*pscf(mu,la) -
//									pmp2(la,nu)*pscf(mu,si);
////							if(debug)std::cout << "separable "<< dsep << std::endl;
//
//							double dtotal = dhf + dsep;
////							double dtotal = dhf;
////							double dtotal = dsep;
//
//
////							if(mu == nu)dtotal = dtotal/2.0; //not noddy
////							if(la == si)dtotal = dtotal/2.0; //not noddy
////							if(std::max(mu,nu) == std::max(la,si) && //not noddy
////								std::min(mu,nu) == std::min(la,si) ) dtotal /=2; //not noddy
////							dtotal *=sfac; //not noddy
////							if(p != r)dtotal *= 2;
//							dtotal *= s1234_deg;
//
//							int ijkl = ((pf*Q.size()+qf)*R.size()+rf)*S.size()+sf;
//
//
//							for (unsigned int di = 0; di < 12; ++di) {
//							  LIBINT2_REALTYPE value;
//							  if (di<3)
//							  {
//							    value = inteval[0].targets[di][ijkl];
////							    deri(di%3,0 ) += value;
//							    eg(P.atom(),di%3)+= value*dtotal;
////							    eg(P.atom(),di%3)+= value;
////							    eg(P.atom(),di%3)+= dtotal;
////								 std::cout << "from binding " << di<< " " << value << std::endl;
//
//							  }
//							  else if (di<6)
//							  {
//							    value = inteval[0].targets[di][ijkl];
////							    deri(di%3,1 ) += value;
//
//							    eg(Q.atom(),di%3)+= value*dtotal;
////							    eg(Q.atom(),di%3)+= value;
////							    eg(Q.atom(),di%3)+= dtotal;
////								 std::cout << "from binding " << di<< " " << value << std::endl;
//							  }
//							  else if (di>=9)
//							  {
//								  value = inteval[0].targets[di-3][ijkl];
////								  deri(di%3,3) += value;
//								   eg(S.atom(),di%3)+= value*dtotal;
////								  eg(S.atom(),di%3)+= value;
////								  eg(S.atom(),di%3)+= dtotal;
//							  }
//							  else // (di>=6 || di<=8)
//							  {
//								  value = -(inteval[0].targets[di-6][ijkl] +
//							        inteval[0].targets[di-3][ijkl] +
//							        inteval[0].targets[di][ijkl]);
////								  deri(di%3,2) += value;
//								  eg(R.atom(),di%3)+= value*dtotal;
////								  eg(R.atom(),di%3)+= value;
////								  eg(R.atom(),di%3)-= 3*dtotal;
////									 std::cout << "from binding " << di<< " " << value << std::endl;
//							  }
//
////							  new_eri.push_back(value);
////								 std::cout << "from binding " << di<< " " << value << std::endl;
//							}//di
//
////							ijkl++;
//
////							std::cout << std::endl;
//
//									}//sf
//								}//rf
//							}//qf
//						}//pf
//
////									std::cout << deri << std::endl;
////						std::cout <<std::endl << eg << std::endl;
////				}//logic
//
//
////				-----------------
////			std::cout << "sfac " << sfac << std::endl;
////			std::cout <<std::endl << eg << std::endl;
////				eg.setZero();
//					}//s
//				}//r
//			}//q
//		}//p
//
////#pragma omp critical
////		std::cout << omp_get_thread_num() << std::endl<<std::endl;
//
//
//		eg /= 8.0; //noddy
//
////#pragma omp critical
////		std::cout << "deri" << std::endl;
////		std::cout <<std::endl << eg << std::endl;
//
//#pragma omp critical
//			{eg_global += eg;
//			global_skipped += skipped;}
//
//
//
//
////#pragma omp barrier
//
//
//
//		LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
//
////		 if(omp_get_thread_num() == 0) std::cout << eg_global << std::endl;
//
////		delete ptr_den;
////		delete inteval;
//		free(inteval);
//
//
//		}//#pragma omp parallel
//
//		std::cout << "quartets skipped " << global_skipped << std::endl;
//		LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();
//
//
//
//}//twobody_deriv_contributions



void twobody_deriv_contributions(size_t nbf, double * ptr_pmp2, double * ptr_pscf,
		Basis basis,
		Molecule molecule,double *ptr_eg_global,
		Eigen::MatrixXd &schwartz_mat){

	typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;


	Eigen::MatrixXd g_coulomb(molecule.size(),3);
	g_coulomb.setZero();
	Eigen::MatrixXd g_exchange(molecule.size(),3);
	g_exchange.setZero();


	MapMatrixXd eg_global(ptr_eg_global,molecule.size(),3);
	MapMatrixXd pmp2(ptr_pmp2,nbf,nbf);
	MapMatrixXd pscf(ptr_pscf,nbf,nbf);

	double cutoff = 1.0e-9;
	double cutoff2 = cutoff/2.0;

	LIBINT2_PREFIXED_NAME(libint2_static_init)();

	typedef unsigned int uint;

	size_t veclen = LIBINT2_MAX_VECLEN;

	//figure out the maximum contraction depth
	//figure out the maximum shell size
	size_t max_contrdepth = 0;
	size_t max_shell_size = 0;
	for(int ishell = 0; ishell < basis.shells().size(); ishell++){
		const Basis::Shell &Q = basis.shells().at(ishell);
		if(Q.K() > max_contrdepth) max_contrdepth = Q.K();
		if(Q.size() > max_shell_size) max_shell_size = Q.size();
//		std::cout << basis.max().L()<< " "<<max_contrdepth << std::endl;
	}//ishell
	size_t max_contrdepth4 = max_contrdepth*max_contrdepth*max_contrdepth*max_contrdepth;

	int libint_lmax = LIBINT2_MAX_AM_ERI1;
	int basis_lmax = basis.max().L();
	if(basis_lmax > libint_lmax) {
		std::cout << "this is not good; Libint is not set up for this high of momentum" <<std::endl
		<< "max Libint momentum l=" << libint_lmax <<std::endl
		<< "but you have functions with moment l=" << basis_lmax <<std::endl
		<< "rebuild libint with higher derivative momentum or" <<std::endl
		<< "select a different basis set" <<std::endl;
	}


	std::cout  <<std::endl
			<< "Contracting separable densities with 4-center Derivative ERIs "
			<< std::endl << std::endl;
	utility::Progress progress;
	progress.reset(basis.shells().size());

	uint global_skipped = 0;
#pragma omp parallel
	{

		Eigen::VectorXd den_vec(basis.max().size()*
				basis.max().size()*
				basis.max().size()*
				basis.max().size() );

		Eigen::VectorXd coulomb_vec(basis.max().size()*
				basis.max().size()*
				basis.max().size()*
				basis.max().size() );
		Eigen::VectorXd exchange_vec(basis.max().size()*
				basis.max().size()*
				basis.max().size()*
				basis.max().size() );

		Eigen::MatrixXd eg(molecule.size(),3);
		eg.setZero();

		Eigen::MatrixXd coulomb(molecule.size(),3);
		coulomb.setZero();
		Eigen::MatrixXd exchange(molecule.size(),3);
		exchange.setZero();


		Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth4);


		LIBINT2_PREFIXED_NAME(libint2_init_eri1)(&inteval[0], // ptr to the first evaluator in the array
				basis_lmax,        // maximum angular momentum
				0);

		uint am[4];

		ShellSet<4u> qset(boost::cref(basis) );
		uint skipped = 0;

		int thread_id = omp_get_thread_num();
		//		(p q | r s)
#pragma omp for schedule(dynamic)
		for(uint p = 0; p < basis.shells().size(); p++){
			const Basis::Shell &P1 = basis.shells().at(p);
			int pam = P1.L();
//			std::cout << omp_get_thread_num() << " " << p << std::endl;

			for(uint q = 0; q <= p; q++){
				const Basis::Shell &Q1 = basis.shells().at(q);
				int qam = Q1.L();

				for(uint r = 0; r <= p; r++){   //not noddy
					const Basis::Shell &R1 = basis.shells().at(r);
					int ram = R1.L();

					const auto s4_max = (p == r) ? q : r;
					for(uint s = 0; s <= s4_max; s++){
						const Basis::Shell &S1 = basis.shells().at(s);
						int sam = S1.L();

						bool debug = 0;

						if(P1.atom() ==  Q1.atom() && R1.atom() == S1.atom() && P1.atom() == R1.atom()) continue;

						auto s12_deg = (p == q) ? 1.0 : 2.0;
						auto s34_deg = (r == s) ? 1.0 : 2.0;
						auto s12_34_deg = (p == r) ? (q == s ? 1.0 : 2.0) : 2.0;
						auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

						if( schwartz_mat(p,q)*schwartz_mat(r,s) < cutoff ){
							skipped++;
							continue;
						}

						int pnew = p;
						int qnew = q;
						int rnew = r;
						int snew = s;

						//new to order shells special for libint
						if(pam < qam)std::swap(pnew,qnew);
						if(ram < sam)std::swap(snew,rnew);
						if( (pam+qam) > (ram+sam) ){
							std::swap(pnew,rnew);
							std::swap(qnew,snew);
						}

						const Basis::Shell &P = basis.shells().at(pnew);
						am[0] = P.L();
						const Basis::Shell &Q = basis.shells().at(qnew);
						am[1] = Q.L();
						const Basis::Shell &R = basis.shells().at(rnew);
						am[2] = R.L();
						const Basis::Shell &S = basis.shells().at(snew);
						am[3] = S.L();

						double max_den = (double)0.0;
						for (int mu = P.start(), ijkl = 0; mu < P.stop(); mu++){
							for (int nu = Q.start(); nu < Q.stop(); nu++){
								for (int la = R.start(); la < R.stop(); la++){
									for (int si = S.start(); si < S.stop(); si++,ijkl++){


										//										int mu = P.start() +pf;
										//										int nu = Q.start() +qf;
										//										int la = R.start() +rf;
										//										int si = S.start() +sf;

										double dhf=0.0;
										dhf = pscf(mu,nu)*pscf(la,si)*4 -
												pscf(mu,si)*pscf(la,nu) -
												pscf(mu,la)*pscf(nu,si);

										double dhf_c = pscf(mu,nu)*pscf(la,si)*4;
										double dhf_e = -pscf(mu,si)*pscf(la,nu) -
												pscf(mu,la)*pscf(nu,si);

										double dsep=0.0;
										dsep = pmp2(mu,nu)*pscf(la,si)*4 +
												pmp2(la,si)*pscf(mu,nu)*4 -
												pmp2(mu,si)*pscf(la,nu) -
												pmp2(mu,la)*pscf(nu,si) -
												pmp2(nu,si)*pscf(mu,la) -
												pmp2(la,nu)*pscf(mu,si);

										double dsep_c = pmp2(mu,nu)*pscf(la,si)*4 +
												pmp2(la,si)*pscf(mu,nu)*4;
										double dsep_e = -pmp2(mu,si)*pscf(la,nu) -
												pmp2(mu,la)*pscf(nu,si) -
												pmp2(nu,si)*pscf(mu,la) -
												pmp2(la,nu)*pscf(mu,si);


										double dtotal = dhf + dsep;

										double dtotal_c = dhf_c + dsep_c;
										double dtotal_e = dhf_e + dsep_e;

										dtotal *= s1234_deg;
										den_vec(ijkl) = dtotal;


										coulomb_vec(ijkl) = dtotal_c*s1234_deg;
										exchange_vec(ijkl) = dtotal_e*s1234_deg;

										if( std::abs(dtotal) > max_den ) max_den = std::abs(dtotal);
									}//s
								}//r
							}//q
						}//p

						if( max_den*schwartz_mat(p,q)*schwartz_mat(r,s) < cutoff2 ){
							skipped++;
							continue;
						}

						//construct quartet set
						qset.add_shell(veclen,pnew,0,0);
						qset.add_shell(veclen,qnew,1,0);
						qset.add_shell(veclen,rnew,2,0);
						qset.add_shell(veclen,snew,3,0);

						inteval[0].contrdepth = P.K()*Q.K()*R.K()*S.K();


						if (LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]] == 0)continue;
						//			qset.print();

						prep_libint2(inteval, qset, 0, 1);

						// (p q | r s)
						// P.L() >= Q.L(), R.L() >= S.L(), R.L()+S.L() >= P.L()+Q.L()
						LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);


						for (int pf = 0, den_pos = 0; pf < P.size(); pf++){

							for (int qf = 0; qf < Q.size(); qf++){ //noddy

								for (int rf = 0; rf < R.size(); rf++){ //noddy

									for (int sf = 0; sf < S.size(); sf++,den_pos++){ //noddy


										//							int mu = P.start() +pf;
										//							int nu = Q.start() +qf;
										//							int la = R.start() +rf;
										//							int si = S.start() +sf;
										//
										//							double dhf=0.0;
										//							dhf = pscf(mu,nu)*pscf(la,si)*4 -
										//									pscf(mu,si)*pscf(la,nu) -
										//									pscf(mu,la)*pscf(nu,si);
										//
										//							double dsep=0.0;
										//							dsep = pmp2(mu,nu)*pscf(la,si)*4 +
										//									pmp2(la,si)*pscf(mu,nu)*4 -
										//									pmp2(mu,si)*pscf(la,nu) -
										//									pmp2(mu,la)*pscf(nu,si) -
										//									pmp2(nu,si)*pscf(mu,la) -
										//									pmp2(la,nu)*pscf(mu,si);
										//
										//							double dtotal = dhf + dsep;
										//							dtotal *= s1234_deg;

										double dtotal = den_vec(den_pos);

										int ijkl = ((pf*Q.size()+qf)*R.size()+rf)*S.size()+sf;


										for (unsigned int di = 0; di < 12; ++di) {
											LIBINT2_REALTYPE value;
											if (di<3)
											{
												value = inteval[0].targets[di][ijkl];
												eg(P.atom(),di%3)+= value*dtotal;
												coulomb(P.atom(),di%3) += value*coulomb_vec(den_pos);
												exchange(P.atom(),di%3) += value*exchange_vec(den_pos);
											}
											else if (di<6)
											{
												value = inteval[0].targets[di][ijkl];
												eg(Q.atom(),di%3)+= value*dtotal;
												coulomb(Q.atom(),di%3) += value*coulomb_vec(den_pos);
												exchange(Q.atom(),di%3) += value*exchange_vec(den_pos);
											}
											else if (di>=9)
											{
												value = inteval[0].targets[di-3][ijkl];
												eg(S.atom(),di%3)+= value*dtotal;
												coulomb(S.atom(),di%3) += value*coulomb_vec(den_pos);
												exchange(S.atom(),di%3) += value*exchange_vec(den_pos);
											}
											else // (di>=6 || di<=8)
											{
												value = -(inteval[0].targets[di-6][ijkl] +
														inteval[0].targets[di-3][ijkl] +
														inteval[0].targets[di][ijkl]);
												eg(R.atom(),di%3)+= value*dtotal;

												coulomb(R.atom(),di%3) += value*coulomb_vec(den_pos);
												exchange(R.atom(),di%3) += value*exchange_vec(den_pos);
											}

										}//di

									}//sf
								}//rf
							}//qf
						}//pf

					}//s
				}//r
			}//q

			if(thread_id == 0)progress.jump(p);

		}//p

		if(thread_id == 0)progress.jump(basis.shells().size());
		//#pragma omp critical
		//		std::cout << omp_get_thread_num() << std::endl<<std::endl;


		eg /= 8.0; //noddy

		coulomb /= 8.0;
		exchange /= 8.0;
		//#pragma omp critical
		//		std::cout << "deri" << std::endl;
		//		std::cout <<std::endl << eg << std::endl;

#pragma omp critical
		{eg_global += eg;
		g_coulomb += coulomb;
		g_exchange += exchange;
		global_skipped += skipped;}




		//#pragma omp barrier



		LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);

		//		 if(omp_get_thread_num() == 0) std::cout << eg_global << std::endl;

		free(inteval);


	}//#pragma omp parallel

//#pragma omp barrier
	std::cout << std::endl << "  skipped " << global_skipped << " blocks" << std::endl;
	LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();

//	std::cout << g_coulomb << std::endl << std::endl;
//	std::cout << g_exchange << std::endl << std::endl;

}//twobody_deriv_contributions

















void twobody_eri(size_t nbf,
		Basis basis,
		Molecule molecule,
		Eigen::MatrixXd &schwartz_mat){

	typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;


	double cutoff = 1.0e-9;
	double cutoff2 = cutoff/2.0;

	LIBINT2_PREFIXED_NAME(libint2_static_init)();

	typedef unsigned int uint;

	size_t veclen = LIBINT2_MAX_VECLEN;

	//figure out the maximum contraction depth
	//figure out the maximum shell size
	size_t max_contrdepth = 0;
	size_t max_shell_size = 0;
	for(int ishell = 0; ishell < basis.shells().size(); ishell++){
		const Basis::Shell &Q = basis.shells().at(ishell);
		if(Q.K() > max_contrdepth) max_contrdepth = Q.K();
		if(Q.size() > max_shell_size) max_shell_size = Q.size();
//		std::cout << basis.max().L()<< " "<<max_contrdepth << std::endl;
	}//ishell
	size_t max_contrdepth4 = max_contrdepth*max_contrdepth*max_contrdepth*max_contrdepth;

	int libint_lmax = LIBINT2_MAX_AM_ERI1;
	int basis_lmax = basis.max().L();
	if(basis_lmax > libint_lmax) {
		std::cout << "this is not good; Libint is not set up for this high of momentum" <<std::endl
		<< "max Libint momentum l=" << libint_lmax <<std::endl
		<< "but you have functions with moment l=" << basis_lmax <<std::endl
		<< "rebuild libint with higher derivative momentum or" <<std::endl
		<< "select a different basis set" <<std::endl;
	}


	std::cout  <<std::endl
			<< "4-center ERIs "
			<< std::endl << std::endl;
	utility::Progress progress;
	progress.reset(basis.shells().size());

	uint global_skipped = 0;
#pragma omp parallel
	{



		Libint_t* inteval = libint2::malloc<Libint_t>(max_contrdepth4);


		LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], // ptr to the first evaluator in the array
				basis_lmax,        // maximum angular momentum
				0);

		uint am[4];

		ShellSet<4u> qset(boost::cref(basis) );
		uint skipped = 0;

		int thread_id = omp_get_thread_num();
		//		(p q | r s)
#pragma omp for schedule(dynamic)
		for(uint p = 0; p < basis.shells().size(); p++){
			const Basis::Shell &P1 = basis.shells().at(p);
			int pam = P1.L();
//			std::cout << omp_get_thread_num() << " " << p << std::endl;

			for(uint q = 0; q <= p; q++){
				const Basis::Shell &Q1 = basis.shells().at(q);
				int qam = Q1.L();

				for(uint r = 0; r <= p; r++){   //not noddy
					const Basis::Shell &R1 = basis.shells().at(r);
					int ram = R1.L();

					const auto s4_max = (p == r) ? q : r;
					for(uint s = 0; s <= s4_max; s++){
						const Basis::Shell &S1 = basis.shells().at(s);
						int sam = S1.L();

						bool debug = 0;

						if(P1.atom() ==  Q1.atom() && R1.atom() == S1.atom() && P1.atom() == R1.atom()) continue;

						auto s12_deg = (p == q) ? 1.0 : 2.0;
						auto s34_deg = (r == s) ? 1.0 : 2.0;
						auto s12_34_deg = (p == r) ? (q == s ? 1.0 : 2.0) : 2.0;
						auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

						if( schwartz_mat(p,q)*schwartz_mat(r,s) < cutoff ){
							skipped++;
							continue;
						}

						int pnew = p;
						int qnew = q;
						int rnew = r;
						int snew = s;

						//new to order shells special for libint
						if(pam < qam)std::swap(pnew,qnew);
						if(ram < sam)std::swap(snew,rnew);
						if( (pam+qam) > (ram+sam) ){
							std::swap(pnew,rnew);
							std::swap(qnew,snew);
						}

						const Basis::Shell &P = basis.shells().at(pnew);
						am[0] = P.L();
						const Basis::Shell &Q = basis.shells().at(qnew);
						am[1] = Q.L();
						const Basis::Shell &R = basis.shells().at(rnew);
						am[2] = R.L();
						const Basis::Shell &S = basis.shells().at(snew);
						am[3] = S.L();





						//construct quartet set
						qset.add_shell(veclen,pnew,0,0);
						qset.add_shell(veclen,qnew,1,0);
						qset.add_shell(veclen,rnew,2,0);
						qset.add_shell(veclen,snew,3,0);

						inteval[0].contrdepth = P.K()*Q.K()*R.K()*S.K();


						if (LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]] == 0)continue;
						//			qset.print();

						prep_libint2(inteval, qset, 0, 1);

						// (p q | r s)
						// P.L() >= Q.L(), R.L() >= S.L(), R.L()+S.L() >= P.L()+Q.L()
						LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](&inteval[0]);

					}//s
				}//r
			}//q

		}//p

#pragma omp critical
		{
		global_skipped += skipped;
		}

		LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);

		//		 if(omp_get_thread_num() == 0) std::cout << eg_global << std::endl;

		free(inteval);


	}//#pragma omp parallel

//#pragma omp barrier
	std::cout << std::endl << "  skipped " << global_skipped << " blocks" << std::endl;
	LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();



}//twobody_deriv_contributions



}//namespace libcchem_libint_interface


