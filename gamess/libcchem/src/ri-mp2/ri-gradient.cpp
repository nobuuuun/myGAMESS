/*
 * ri-gradient.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: luke
 */


#include "god.hpp"

#include "ri-gradient.hpp"
#include <metric.hpp>

//these are routines needed to perform rimp2
#include "ri-mp2/ri-integrals.hpp"
#include "ri-async.hpp"
#include "ri-energy-terms.hpp" //temp
#include "array/hdf5.hpp"
#include "utility/timer.hpp"
#include <iostream>


#include <Eigen/Dense>


//#include <cmath> //std::abs


#include <math.hpp>
//for temporary libint interface
#if HAVE_LIBINT
#include "libint_interface/libcchem_libint.h"
#endif

#include <gradient-helpers.hpp>
#include <three-center-work.hpp>
#include <ri-zvector.hpp>
#include <ri-lagrangian.hpp>

#if HAVE_CUBLAS
#include <deviceGradientEngine.hpp>
#endif

//for reading files
#include<fstream>

#include <JK.hpp>
#include <DFIntGradient.hpp>

namespace cchem{
namespace rimp2_gradient{

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> AlignedMapMatrixXd;
typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef const Eigen::Map<const Eigen::MatrixXd,Eigen::AutoAlign> ConstMapMatrixXd;
//typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,Eigen::AutoAlign> MapRowMatrixXd;


typedef Eigen::Map<Eigen::VectorXd,Eigen::AutoAlign> MapVectorXd;

typedef ::rysq::TwoCenterEri TWO_ERI;
//typedef ::rysq::TwoCenterDerivativeEri INT_EXP;

typedef ::rysq::ThreeCenterEri THREE_ERI;
//typedef ::rysq::ThreeCenterDerivativeEri THREE_ERI;


}
}



// potential pitfalls:
//
//    1) the wmat[I] matrix component runs over occ-occ and not act-occ
//				see EQ [181] in Aitkens paper and Wik(2) in head-gordon paper


namespace cchem{

template<typename T>
void print_matrix(T & mat){

std::cout << std::right << std::setprecision(10);
for (int i = 0; i < mat.cols() ; i+=5){

	int ij_max = std::min( (int)mat.cols(), (int)(i+5) );
	std::cout << std::endl<< std::endl<< std::endl<< std::endl;
	for (int j = i; j < ij_max; j++){
		if(j == i)std::cout << std::setw(8) << j ;
		if(j > i)std::cout << std::setw(14) << j ;
	}
	std::cout<<std::endl << std::endl;
	std::cout << mat.block(0,i,mat.rows(),ij_max-i);
	std::cout<<std::endl << std::endl;
}//i

}//print_matrix


double rimp2_gradient::gradient(Wavefunction wf, Runtime &rt, Molecule molecule, double *EG) {



    //lukedef (3-locations)
    //    #undef HAVE_CUBLAS

	//	things to try to fix problems:
	//	  1) I turn off modification to the frozen/active blocks of wmat (maybe turn that back on?)
	//    2) is there a sorting problem
	//    3) are you passing in hardwired orbitals?

	utility::timer globalTimer;
	globalTimer.reset();


	BOOST_AUTO(const &basis, wf.basis());
	BOOST_AUTO(const &shells, basis.shells());

	BOOST_AUTO(const &auxbasis, wf.auxbasis());
	BOOST_AUTO(const &auxshells, auxbasis.shells());

	wf.sort();
	wf.auxsort();

//	Basis LibintBasis = wf.basis();
//	LibintBasis.sort();


	MapMatrixXd egGlobal(EG, molecule.size(), 3);
	egGlobal.setZero();

	//debugging
	bool gradientDebug = 0;

	const double pseudoSquareRootCutoff = rt.get<double>("/rimp2/sqrt/cutoff", 1e-10);

	const double pseudoTol = rt.get<double>("/rimp2/inversion/cutoff", 1e-10);

	double E = 0;
	double e_term[7] = {0};

//	double cutoff = rt.get<double>("/mp2/integrals/cutoff", 1e-10);
//	integrals::Screening screening(basis, cutoff);



	//size of orbitals basis set
	size_t N = basis.size();



	//number of occupied orbitals
	size_t no = wf.occupied().size();


	//number of active orbitals (doubly and singly occupied)
	//    size_t no = wf.active().size();
	size_t na = wf.active().size();
	size_t nv = wf.virtuals().size();

	//number of frozen orbitals
	size_t nf = no - na;

	//number of singly occupied orbitals
	size_t ns = wf.single().size();

	//number of singly+virtual orbitals
	size_t nvs = ns+nv;

	//number of doubly occupied orbitals
	size_t nd = wf.double_occ().size();

	//set up singly+virtual occupied orbitals
	wf.set_single_and_virtuals(wf.double_occ().size(),wf.orbitals().size());

	//singly occupied and virtual orbital coefficients
	Matrix Cv = trans(wf.C(wf.single_virtual()));  //zapt

	//active orbitals coefficients (doubly occupied and singly occupied)
	//    Matrix Ca = trans(wf.C(wf.active()));
	Matrix Ca = trans(wf.C(wf.occupied()));



	int spherical = wf.do_spherical();

	//parallel environment
	Parallel pe;
	// suppress output
	if (pe.rank() != 0) rt.set_cout(0);
	BOOST_AUTO(cout, rt.cout());

	//setup the cartesian->spherical ERI transformation matrix
	std::vector<MapMatrixXd> sph_c;
	std::vector<double *> pSphTrans;
	int max_auxbasis_shell_size = auxbasis.max().size();
	if(spherical){
		cout << std::endl<< "spherical auxiliary basis set will be used to construct (L|M)" << std::endl;
		for(int lang = 0; lang<7; lang++){
			int dim = (lang+1)*(lang+2)/2;
			MapMatrixXd temp( wf.spherical_coeff_ptr(lang),dim,dim);
			sph_c.push_back(temp);
			pSphTrans.push_back( wf.spherical_coeff_ptr(lang) );
		}//lang
		max_auxbasis_shell_size = 13; //I : (2*6+1) <- (2L+1)
	}//(spherical)

	//determine size of auxiliary basis set
	size_t nl = auxbasis.size();
	if(spherical){
		nl = auxbasis.spherical_size();
		cout << auxbasis.size() << " cartesian auxiliary functions ---> "
				<< auxbasis.spherical_size() << " spherical auxiliary functions"
				<< std::endl << std::endl;
	}//(spherical)


		cout << "frozen:     " << nf << std::endl;
		cout << "double:     " << nd << std::endl;
		cout << "occupied:   " << no << std::endl;
		cout << "single:     " << ns << std::endl;
		cout << "active:     " << na << std::endl;
		cout << "virtual:    " << nv << std::endl;
		cout << "atomic:     " << N  << std::endl;
		cout << "auxiliary:  " << nl << std::endl;
		cout << "pseudo square root cutoff : " << pseudoSquareRootCutoff << std::endl;
		cout << "pseudo inversion tolerance : "	<< pseudoTol << std::endl << std::endl;
//		cout << "eri cutoff: " << screening.value() << std::endl<< std::endl;


	Runtime::Memory &memory = rt.memory();
	double used = 0;

	size_t max_threads = omp_get_max_threads();
	size_t nwords = (memory.available()/sizeof(double));

	// runtime may be different, ensure consistency
	pe.broadcast(&max_threads, 1, 0);
	pe.broadcast(&nwords, 1, 0);


	// for (int i = 0; i < pe.size(); i++){
	//     if(pe.rank() == i)
	// 	std::cout <<  "this broadcast worked!!! "  << max_threads << " " << nwords << std::endl;
	//     pe.barrier();
	// }



	{

		BOOST_AUTO(const &mem_twoc, cchem::ri::AuxiliaryTwoCenterInt<TWO_ERI>::memory(auxbasis));
		BOOST_AUTO(const &mem_threec, cchem::ri::AuxiliaryThreeCenterInt<THREE_ERI>::memory(auxbasis, basis, Ca, Cv));

		size_t required = 0;
		foreach (size_t n, mem_twoc) required += n;
		required *= max_threads;
		required +=nl*nl; //for L^-1
		if (required > nwords)
			throw cchem::exception("not enough memory to run RI-MP2");

		required = 0;
		foreach (size_t n, mem_threec) required += n;
		required *= max_threads;
		required +=nl*nl; //for L^-1 - carried over from two-center Eri
		if (required > nwords)
			throw cchem::exception("not enough memory to run RI-MP2");

		nwords -= required;  //more thought is need here

	}


	struct {
	    utility::timer::value_type two_center,cholesky,inversion;
	    utility::timer::value_type three_center;
	    utility::timer::value_type t2_read,t3_trans,t3_write;
	    utility::timer::value_type eri_formation,energy_accum,eri_formationii,eri_formationij,eri;

	    utility::timer::value_type zvector, deri_4center;
	    utility::timer::value_type pmatGamma;
	    utility::timer::value_type recompute;
	    utility::timer::value_type remaining_density;
	    utility::timer::value_type lag_part;
	    utility::timer::value_type overlap_derivative;
	    utility::timer::value_type one_derivative;
	    utility::timer::value_type twoc_derivative;
	    utility::timer::value_type threec_derivative;

	} profile = {};

	utility::timer timer;
	utility::timer timerTest;


	::cchem::rimp2::device * pGPU = NULL;


	#if HAVE_CUBLAS


 /////////////////////////////////////////////
 // set up GPU devices
 /////////////////////////////////////////////
		
 //blocks per cuda kernel
 int cudaBlocks = 8; //do not mess with yet
 //threads per block
 int cudaThreadsPB = 64; //do not mess with yet

 int nstreams = ::cchem::rimp2::device::getNumberOfStreams();
 int ndevice = ::cchem::rimp2::device::getNumberOfDevices();

 ::cchem::rimp2::device GPU(ndevice,nstreams, cudaBlocks, cudaThreadsPB, (pe.rank()==0) );

 pGPU = &GPU;

 int cudaThreadsPerDevice = cudaBlocks*cudaThreadsPB;

 const size_t maxGPUWords = pGPU->maxWords();

 size_t freeWords = maxGPUWords;

 size_t NMBDevice = freeWords*sizeof(double)/1000/1000;
 size_t NMBDeviceTemp = NMBDevice;

 for(size_t ip = 0; ip < pe.size(); ip++){

     if(ip == pe.rank()){
 	 pe.broadcast(&NMBDevice, 1, ip);
     }else{
  	 pe.broadcast(&NMBDeviceTemp, 1, ip);
     }
     if(NMBDeviceTemp < NMBDevice) NMBDevice = NMBDeviceTemp;

 } //ip



#endif


	/*
	       Coulomb metric work

	 1) build two-center ERIs
	 2) pseudo invert two-center ERI matrix (P|Q)^-0.5
	 3) form (P|Q)^-1

	 */

	//(P|Q)^-1 storage - no initialization necessary
	double *pPqm1 = new double[nl*nl];

	//(P|Q)^-1/2 storage - no initialization necessary
	double *pPqm12 = new double[nl*nl];

	{  //scope: coulomb metric
		timer.reset();

		//only one copy need for all threads --so allocate here--
		//for ERI storage - no initialization necessary
		double *pTwoC = new double[nl*nl];

		::cchem::rimp2_gradient::detail::Metric
		 metric(boost::cref(auxbasis), pe, spherical, nl, pTwoC,pPqm1);

		//build (A|B) matrix
		metric.BuildTwoCenterEri(sph_c);

		profile.two_center += timer;

		timer.reset();

		//form (P|Q)^-0.5
		metric.pseudoInvertedSqrt(pseudoTol, pPqm12, pTwoC, (pe.rank()==0) );

		profile.inversion += timer;

		//-----use LDLT--- this does everything to form inverse metric
		//metric.DecomposeLDLT(); profile.inversion += timer; timer.reset();

		//-----use LLT--- need a couple additional steps to form inverse metric
		//metric.DecomposeLLT();     profile.cholesky += timer; timer.reset();
		//metric.InvertL();
		//metric.Lm1_to_pqm1();      profile.inversion += timer;



		//form (PQ)^-1 from (P|Q)^-0.5
		//  nb: A = (P|Q)^-0.5 * (Q|P)^-0.5
		//      A^-1 = (Q|P)^-0.5 * (P|Q)^-0.5  (nb swapped transpose)
		cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
				nl,nl,nl,
				1.0, pPqm12, nl,
				pPqm12, nl,
				0.0, pPqm1, nl);

		delete [] pTwoC;

	} //scope: coulomb metric


	cout << "time for two-center eri            "<< profile.two_center << std::endl;
	//	std::cout << "time for Cholesky Decomposition    "<< profile.cholesky <<std::endl;
	cout << "time to form   (P|Q)^-1            "<< profile.inversion <<std::endl<< std::endl;








	//set up storage
	//figure out the largest amount of memory allocated at anyone time for arrays:

	size_t initialArrays = 3*no*nvs*nl + no*N*nl + 2*N*(N+1)/2*nl;
	size_t secondArrays = no*N*nl + N*(N+1)/2*nl + (no+nv)*nl + N*N*nl;


	cout << "InitialArrays(MB) = " << initialArrays*sizeof(double)/1e6 << std::endl
	     << "     Array 1: "  << no*nvs*nl*sizeof(double)/1e6 << std::endl
	     << "     Array 2: "  << no*nvs*nl*sizeof(double)/1e6 << std::endl
	     << "     Array 3: "  << no*nvs*nl*sizeof(double)/1e6 << std::endl
	     << "     Array 4: "  << no*N*nl*sizeof(double)/1e6 << std::endl
	     << "     Array 5: "  << N*(N+1)/2*nl*sizeof(double)/1e6 << std::endl
	     << "     Array 6: "  << N*(N+1)/2*nl*sizeof(double)/1e6 << std::endl

	     << "SecondArrays(MB) = " << secondArrays*sizeof(double)/1e6 << std::endl
             << "     Array 1: "  << no*N*nl*sizeof(double)/1e6 << std::endl
             << "     Array 2: "  << N*(N+1)/2*nl*sizeof(double)/1e6 << std::endl
             << "     Array 3: "  << no*(no+nv)*nl*sizeof(double)/1e6 << std::endl
             << "     Array 4: "  << N*N*nl*sizeof(double)/1e6 << std::endl


	     << "Largest amount of memory allocated for arrays at one time: "
	     << std::max(initialArrays,secondArrays)*sizeof(double)/1e6
	     << " MB" <<std::endl;

	cout << "Creating arrays: " << std::endl;
	
	//lm1:  electronic hessian (zvector)
	//          vooo and vvvo contributions to the lagrangian
	//          energy weighted density matrix construction
	//bare:     JK gradients 
	size_t dimsMuNuSym[2] = {  N*(N+1)/2, nl };
	size_t chunkMuNuSym[2] = {  dimsMuNuSym[0], 1 };
	rt.arrays().allocate<double,2>("rimp2.lm1MuNuSym", dimsMuNuSym, pe, chunkMuNuSym );
	rt.arrays().allocate<double,2>("rimp2.bareMuNuSym", dimsMuNuSym, pe, chunkMuNuSym );
	
	
	size_t dims[2] = {  no*(nv+ns), nl };
	size_t chunk[2] = {  dims[0], 1 };
	//		size_t chunk[2] = {  (nv+ns), nl }; //idea?
	rt.arrays().allocate<double,2>("rimp2.gamma_ia_q", dims, pe, chunk); //forming with (P|Q)-1
	rt.arrays().allocate<double,2>("rimp2.bare_ia", dims, pe, chunk);
	rt.arrays().allocate<double,2>("rimp2.coeff_ia", dims, pe, chunk);
	
	//for mixed ao-mo gamma (GAMMA_i-nu_P)
	size_t dims_gamma_inu_P[2] = {  no*N, nl };
	size_t chunk_gamma_inu_P[2] = {  dims_gamma_inu_P[0], 1 };
	rt.arrays().allocate<double,2>("rimp2.gamma_inu_P", dims_gamma_inu_P, pe, chunk_gamma_inu_P);

	//(ia|Q)
	Array<double> *BARE_IA = rt.arrays().find< Array<double> >("rimp2.bare_ia");
	rt.cout() << "     " << *BARE_IA << std::endl;

	//C_ia^Q coefficients [ i.e. (ia|Q)(Q|P)-1 ]
	Array<double> *CIAQ = rt.arrays().find< Array<double> >("rimp2.coeff_ia");
	rt.cout() << "     " << *CIAQ << std::endl;

	//GAMMA_ia^Q = t_ia^ab * C_ia^Q
	Array<double> *GAMMA_IA_Q = rt.arrays().find< Array<double> >("rimp2.gamma_ia_q");
	rt.cout() << "     " << *GAMMA_IA_Q << std::endl;

	//gamma_i_nu_P = gamma_ia^Q*virt(nv,N)
	Array<double> *GAMMA_INU_P = rt.arrays().find< Array<double> >("rimp2.gamma_inu_P");
	rt.cout() << "     " << *GAMMA_INU_P << std::endl;

	//BARE_MN_SYM (mu nu | Q) --> (N_*(N_+1)/2 | nl)
	Array<double> *BARE_MN_SYM = rt.arrays().find< Array<double> >("rimp2.bareMuNuSym");
	rt.cout() << "     " << *BARE_MN_SYM << std::endl;

	//LM1_MN_SYM (mu nu | Q)(Q|P)^-1/2  OR  (mu nu | Q)L^-1
	Array<double> *LM1_MN_SYM = rt.arrays().find< Array<double> >("rimp2.lm1MuNuSym");
	rt.cout() << "     " << *LM1_MN_SYM << std::endl << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////
	// 1) evaluate three-center ERIs
	// 2) perform 1/3 and 2/3 AO to MO transformations (AO->MO transformation occurs here)
	////////////////////////////////////////////////////////////////////////////////////////

	pe.barrier();
	
	pe.task().reset(); // if(pe.rank() != 0)goto skip0;

	
	timer.reset();
	cout << "Starting Three-center ERIs with " 
	     << omp_get_max_threads() 
	     << " OpenMP threads. " << std::endl;
	
	{   //scope: three-center Eri and transform

	    ::cchem::rimp2_gradient::detail::ThreeCenterWork
		three_center_eri(boost::cref(basis),boost::cref(auxbasis),
				 boost::cref(Ca),boost::cref(Cv),
				 pe, spherical, max_auxbasis_shell_size,no,ns,nv,N,
				 BARE_IA, NULL, NULL, BARE_MN_SYM);
	    
	    //build (L|IA) integrals
	    //build (L| mu nu) integrals
	    three_center_eri.BuildMOInts(sph_c);

	    pe.barrier();
	    
	    profile.three_center += timer;
	    cout<< "     Time for three-center ERI and 1/3+2/3 Transformation: " 
		<< profile.three_center  
		<< std::endl << std::endl;
	    
	} //scope: three-center Eri and transform
	
	//skip0:;


	{//scope: gradient work starts here

	    ::cchem::rimp2_gradient::detail::gradientTerms egTerms( molecule.size() );

	    //////////////////////////////////////////
	    //work out global memory storage
	    //  1) MO-Lagrangian (lag_mo)
	    //  2) relaxed density (pmat)
	    //  3) energy weighted density (wmat)
	    //////////////////////////////////////////
	    
	    //Lagrangian(nv,no)
	    double *ptr_lag_mo = new double[nv*no];
	    MapMatrixXd lag_mo(ptr_lag_mo,nv,no);
	    memset(ptr_lag_mo, 0, nv*no*sizeof(double) );
	    
	    //pmat(basis.size(),basis.size())
	    double *ptr_pmat = new double[N*N];
	    MapMatrixXd pmat(ptr_pmat, N, N );
	    memset(ptr_pmat, 0, N*N*sizeof(double) );
	    
	    
	    //wmat(basis.size(),basis.size())
	    double *ptr_wmat = new double[N*N];
	    MapMatrixXd wmat(ptr_wmat, N, N );
	    memset(ptr_wmat, 0, N*N*sizeof(double) );
	    
	    //storage for gamma_rs (two-centered force non-seperable 'density')
	    double *ptr_gamma_rs2 = new double[nl*nl];
	    MapMatrixXd gamma_rs2(ptr_gamma_rs2,nl,nl);
	    memset(ptr_gamma_rs2, 0, nl*nl*sizeof(double) );
	    
	    //////////////////////////////////////////////////////////////
	    //build RI coefficient matrix C_ia^Q - store on disk or in DM.
	    //build RI B_ia^Q - store on disk or in DM.
	    //   nb. C_ia^Q = (ia|P)(P|Q)^(-1)
	    //       B_ia^Q = (ia|P)L_PQ^(-1)
	    //////////////////////////////////////////////////////////////
	    
	    {//scope:  build C_ia^Q and B_ia^Q
		
		size_t NMB = 5000;
		bool read_restrict = 1; bool async = 1;
		bool ReadRow = 1;
		bool MagicSwap = 0;
		bool loopRestrict = 0;
		
		//build C_ia^Q coefficient matrix here (i.e. BARE_IA -> CIAQ ):
		//  (ia|Q) -> (ia|Q)(Q|R)-1
		{//scope: C_ia^Q 
		    timer.reset();
		    bool printDebug = 0;//(pe.rank() == 0);
		    
		    ::cchem::rimp2_gradient::detail::DFAsync 
			  MakeCIAQ(NMB,
				   BARE_IA, no, 0, 1, 0,
				   CIAQ, no, 0, 1, 0,
				   read_restrict, loopRestrict, ReadRow, 
				   MagicSwap, pe, printDebug);
		    
		    ::cchem::rimp2_gradient::detail::transform_functor_new2
			  Functor( pPqm1 );
		    
		    MakeCIAQ.DoOpAsync_R_W_NODE_PARALLEL( Functor );
		    
		    cout << "(ia|P) -> (ia|P)(P|Q)^-1 3/3 Transformation: " 
			 << timer << std::endl;
		}//scope: C_ia^Q
		
		
		//build  B_mn^Q (i.e. BARE_MN_SYM -> LM1_MN_SYM):
		//  (mu nu| Q) -> (mu nu|Q)L-1
		{//scope: B_mn^Q
		    timer.reset();
		    bool printDebug = 0;//(pe.rank() == 0);
		    
		    ::cchem::rimp2_gradient::detail::DFAsync 
			  MakeBMNQ2(NMB,
				    BARE_MN_SYM, N*(N+1)/2, 0, 1, 0,
				    LM1_MN_SYM, N*(N+1)/2, 0, 1, 0,
				    read_restrict, loopRestrict,  ReadRow, 
				    MagicSwap, pe, printDebug);
		    
		    ::cchem::rimp2_gradient::detail::transform_functor_new2
			  Functor( pPqm12 );
		    
		    MakeBMNQ2.DoOpAsync_R_W_NODE_PARALLEL( Functor );
		    
		    cout << "(mu nu|P) -> (mu nu|P)Lm1 3/3 Transformation: " 
			 << timer << std::endl << std::endl;
		}//scope: B_mn^Q
		
		
	    } //scope: build C_ia^Q and B_mn^Q
	    
	    pe.barrier();







	    timer.reset();
	    //////////////////////////////////////////
	    //build gamma_ia_P
	    //back-transform gamma_ia_Q -> gamma_inu_Q
	    //build pvv block of pmat 'P(2)'
	    //build pij block of pmat 'P(2)'
	    //build gamma^RS
	    //compute energy
	    //////////////////////////////////////////
	    
#if !HAVE_CUBLAS
	    
	    pe.task().reset();
	    
	    pe.barrier();
	    
	    
	    {//scope: Pmat, gamma_ia_P, gamma_inu_P
		
		timer.reset();
		
		BOOST_AUTO(const &ea, wf.e(wf.double_occ()));
		BOOST_AUTO(const &ev, wf.e(wf.single_virtual()));
		
		
		//storage for GAMMA_ia^P_long
		double *ptr_gamma_ia_P_long = new double[na*nv*nl];
		std::vector<double *> pGammaLong;
		for (int i = 0; i < na; i++){
		    pGammaLong.push_back(&ptr_gamma_ia_P_long[i*nv*nl]);
		    memset(&ptr_gamma_ia_P_long[i*nv*nl],0,nv*nl*sizeof(double));
		}//i
		
		size_t reqMem = 
		    ::cchem::rimp2_gradient::detail::CPUPmatGammaFunctor::CPUMem(nv,na,nl,no);
		
		cout << "CPUPmatGammaFunctor: Each node allocating "
		     << reqMem/1000/1000 
		     << " MB for internal structures" << std::endl;
		
		
		double energy = 0.0;
		size_t NMB = 17000; bool read_restrict = 0; bool ReadRow = 1;bool MagicSwap = 1;
		
		//			NMB = 100;
		//			NMB = 200;
		//			NMB=5;
		//			NMB=20;
		//this is the total minimum size of buffers needed for asynchronous I/O
		size_t minNMB = 2*(na+1)*nv*nl*sizeof(double)/1000/1000 +1;
		
		cout << std::endl
		     << "Need at least " 
		     << minNMB
		     << " MB for PMAT and GAMMA asynchronous I/O" << std::endl;
		
		if( NMB < minNMB ){ cout << "Problem!!! Adjusting buffer size "
					 << "for PMAT and GAMMA work" << std::endl;
		    cout << "NMB is increased from " << NMB
			 << " to " << minNMB << " MWORDS" << std::endl;
		    NMB = minNMB;
		}//( NMB < minNMB )
		
		
		{//scope
		    
		    bool printDebug = 0;//(pe.rank() == 0);
		    bool loopRestrict = 0;
		    ::cchem::rimp2_gradient::detail::DFAsync 
			  GAMMA(NMB,
				BARE_IA, no, nf, 1, 0,
				CIAQ,    no, nf, 1, 0,
				read_restrict, loopRestrict, ReadRow, MagicSwap, pe, printDebug);
		    size_t maxNA2B = GAMMA.getMaxNAB(1);
		    
		    
		    
		    ::cchem::rimp2_gradient::detail::CPUPmatGammaFunctor 
			  GammaPmatFunctor(
					   nf,no,nv,nl,
					   &(ea.data()[0]),&(ev.data()[nd]),
					   energy,
					   pmat,
					   pGammaLong,
					   ptr_gamma_rs2,
					   pPqm1);
		    
		    
		    GAMMA.DoOpAsync_RR(
				       GammaPmatFunctor
				       ); //GAMMA.DoOpAsync_RR
		    
		}//scope
     
		
		pe.reduce("+",&energy, (size_t)1 ); 
		cout << std::setprecision(10) << std::fixed;
		if(pe.rank() == 0)
		    cout << "energy outside" << energy << std::endl << std::flush;
		e_term[0] = energy;
		pe.reduce("+",ptr_gamma_rs2, (size_t)nl*nl); 
		pe.reduce("+",pmat.data(), (size_t)N*N); 
		pe.reduce("+",ptr_gamma_ia_P_long, (size_t)na*nv*nl); 
		
		
		profile.pmatGamma += timer;
		cout << "CPU time for (occ-occ) and (virt-virt) " 
		     << "relaxed density blocks and Gamma_ia^Q: " 
		     << profile.pmatGamma 
		     << std::endl;
		
		
		timer.reset();
		//storage for GAMMA_inu^P : '1/2' back transformed
		double *ptr_gamma_inu_P = new double[N*nl];
		
		for (int iocc = nf+pe.rank(); iocc < no; iocc+=pe.size() ){
		    
		    //dump gamma_ia_P here
		    size_t start[] = { iocc*nv, 0 };
		    size_t finish[] = { (iocc+1)*nv, nl };
		    GAMMA_IA_Q->put(pGammaLong[iocc-nf], start, finish);
		    
		    //back tranform gamma_ia_P --> gamma_i-nu_P //(N,nl) = (N,nv) (nv,nl)
		    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
				N,nl,nv,
				1.0, Cv.data().begin(), nv,
				pGammaLong[iocc-nf],nv,
				0.0, ptr_gamma_inu_P, N);
		    
		    //dump 1/2 back-transformed gamma_inu_P here
		    size_t start_inu[] = { iocc*N, 0};
		    size_t finish_inu[] = { (iocc+1)*N,nl};
		    GAMMA_INU_P->put(ptr_gamma_inu_P, start_inu, finish_inu);
		}//iocc
		delete [] ptr_gamma_inu_P;
		delete [] ptr_gamma_ia_P_long;
		
		pe.barrier();
		
	    }//scope: Pmat, gamma_ia_P, gamma_inu_P
	    
	    
	    profile.pmatGamma += timer;
	    cout << "time for write and trans+write: " 
		 << timer << std::endl<< std::endl;
	    
	    pe.barrier();
	    
#endif //!HAVE_CUBLAS
	    
#if HAVE_CUBLAS

	    { //new scope
		timer.reset();
		
		BOOST_AUTO(const &hostea, wf.e(wf.double_occ()));
		BOOST_AUTO(const &hostev, wf.e(wf.single_virtual()));
		
		//storage for GAMMA_ia^P_long
		double *ptr_gamma_ia_P_long = new double[na*nv*nl];
		std::vector<double *> pHostGammaLong;
		for (int i = 0; i < na; i++){
		    pHostGammaLong.push_back(&ptr_gamma_ia_P_long[i*nv*nl]);
		    memset(&ptr_gamma_ia_P_long[i*nv*nl],0,nv*nl*sizeof(double));
		}//i
		
		
		int nThreads = omp_get_max_threads();
		omp_set_num_threads(1);
		
		size_t reqMem = 
		    ::cchem::rimp2_gradient::detail::GPUPmatGammaFunctor::
		    GPUMem(nl,nv,no,pGPU);
		
		cout << "GPU: Using " 
		     << (size_t)(reqMem/1e6 +1)
		     << " MB per device"
		     << " out of "
		     << NMBDevice
		     << " MB avaible (max per device)" << std::flush << std::endl;
		
		if( (size_t)(reqMem/1e6) > NMBDevice)
		    cout << "Not enough memory on GPU device. "
			 << " Try reducing the number of CUDA streams."
			 << std::flush << std::endl;
		
		
		
		size_t NMB = 17000; bool read_restrict = 0; bool ReadRow = 1;bool MagicSwap = 1;
		// NMB=14000;
		bool printDebug= 0;//(pe.rank() == 0);
		//NMB=5;
		bool loopRestrict = 0;
		::cchem::rimp2_gradient::detail::DFAsync 
		      GammaDevice(NMB,
				  BARE_IA, no, nf, 1, 0,
				  CIAQ,    no, nf, 1, 0,
				  read_restrict, loopRestrict, ReadRow, MagicSwap, pe, printDebug);
		size_t maxNA2B = GammaDevice.getMaxNAB(1);
		
		double energy;
		
		{//scope
		    
		    ::cchem::rimp2_gradient::detail::GPUPmatGammaFunctor
			GPUFunctor(nv, no, nl, nf,
				   &(hostea.data()[0]),&(hostev.data()[no]),
				   pGPU,
				   pPqm1,
				   pmat,
				   ptr_gamma_rs2,
				   energy,
				   pHostGammaLong);
		    
		    GammaDevice.DoOpAsync_RR ( 
		     			      GPUFunctor 
		     			       );
		    
		}//end scope
		
		delete [] pPqm12;		
		
		pe.reduce("+",&energy, (size_t)1 ); 
		cout << std::setprecision(10);
		if(pe.rank() == 0)
		    cout << "Correlation Energy: " << energy << std::endl << std::endl;
		e_term[0] = energy;
		pe.reduce("+",ptr_gamma_rs2, (size_t)nl*nl); 
		pe.reduce("+",pmat.data(), (size_t)N*N); 
		pe.reduce("+",ptr_gamma_ia_P_long, (size_t)na*nv*nl); 
		
		
		profile.pmatGamma += timer;
		cout << "GPU time for (occ-occ) and (virt-virt) "
		     << "relaxed density blocks and Gamma_ia^Q: " 
		     << profile.pmatGamma 
		     << std::endl;
		
		
		timer.reset();
		//storage for GAMMA_inu^P : '1/2' back transformed
		double *ptr_gamma_inu_P = new double[N*nl];
		
		for (int iocc = nf+pe.rank(); iocc < no; iocc+=pe.size() ){     
		    
		    //dump gamma_ia_P here
		    size_t start[] = { iocc*nv, 0 };
		    size_t finish[] = { (iocc+1)*nv, nl };
		    GAMMA_IA_Q->put(pHostGammaLong[iocc-nf], start, finish);
		    
		    //back tranform gamma_ia_P --> gamma_i-nu_P //(N,nl) = (N,nv) (nv,nl)
		    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
				N,nl,nv,
				1.0, Cv.data().begin(), nv,
				pHostGammaLong[iocc-nf],nv,
				0.0, ptr_gamma_inu_P, N);
		    
		    //dump 1/2 back-transformed gamma_inu_P here
		    size_t start_inu[] = { iocc*N, 0};
		    size_t finish_inu[] = { (iocc+1)*N,nl};
		    GAMMA_INU_P->put(ptr_gamma_inu_P, start_inu, finish_inu);
		}//iocc
		
		delete [] ptr_gamma_inu_P;
		delete [] ptr_gamma_ia_P_long;
		
		profile.pmatGamma += timer;
		cout << "     Time to write Gamma_ia^Q and trans+write Gamma_inu^P: " 
		     << timer << std::endl<< std::endl;
		
		omp_set_num_threads(nThreads);
		
	    }//new scope
	    
	    pe.barrier();
	    
#endif //HAVE_CUBLAS



	    timer.reset();
	    cout << "Starting mixed Lagrangian." << std::endl;
	    /////////////////////////////////////////////////
	    //recompute three-center eri (in AO basis. There is NO AO->MO transformation)
	    //build AO/MO Lagrangian L_mu_i(1)  and L_a_nu(2)
	    //construct wvo block of W(2)
	    //construct wao block of W(2)
	    //construct wvv block of W(2)
	    //construct paf block of P(2)
	    //transform L_mu_i(1) -> L_ai(1)
	    //transform L_a_nu(2) -> L_ai(2)
	    /////////////////////////////////////////////////
	    
	    
	    {//scope: build vvoo components of lagrangian
		
		BOOST_AUTO(const &ea, wf.e(wf.double_occ()));
		
		MapMatrixXd virt_coeff_mat(  Cv.data().begin(), Cv.size1(), Cv.size2() );
		MapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() );
		
		//storage of 'mixed' AO/MO lagrangian
		double *ptr_lag_a_nu = new double[nv*N];
		MapMatrixXd lag_a_nu(ptr_lag_a_nu, nv, N);
		memset(ptr_lag_a_nu, 0, nv*N*sizeof(double) );
		
		//storage of 'mixed' AO/MO lagrangian
		double *ptr_lag_mu_i = new double[N*no];
		MapMatrixXd lag_mu_i(ptr_lag_mu_i, N, no);
		memset(ptr_lag_mu_i, 0, N*no*sizeof(double) );
		
		//storage for GAMMA_inu^P : '1/2' back transformed
		double *ptr_gamma_inu_P = new double[N*no];
		MapMatrixXd gamma_inu_P(ptr_gamma_inu_P, N, no);
		
		//storage for GAMMA_ia^P
		double *ptr_gamma_ia_P = new double[nv*no];
		MapMatrixXd gamma_ia_P(ptr_gamma_ia_P, nv, no);
		
		//////////////////////////////////////////////////////////////////////////////////
		// 1) re-evaluate three-center ERIs (to build mixed 'AO/MO' lag_mu_i and lag_a_nu)
		//		*these are re-evaluated since the GAMMA_inu^P and the GAMMA_ia^P  
		//               are needed
		//	      to compute the 'mixed' lag_mu_i/lag_a_nu Lagrangians
		// 2) perform 1/3 AO to occ MO transformations (this is only done to form lag_a_nu)
		//////////////////////////////////////////////////////////////////////////////////
    
		pe.task().reset(); 
    
		{   //scope: recompute three-center ERI --> mixed Lag, W[I], P[core,active], MO lag
		    
		    ::cchem::rimp2_gradient::detail::ThreeCenterWork
			lag_work(boost::cref(basis),boost::cref(auxbasis),
				 boost::cref(Ca),boost::cref(Cv),
				 pe, spherical, max_auxbasis_shell_size,no,ns,nv,N,
				 BARE_IA,
				 GAMMA_INU_P,GAMMA_IA_Q,NULL);
		    
		    lag_work.BuildMixedLagrangian(ptr_gamma_ia_P,ptr_gamma_inu_P,
						  sph_c,lag_mu_i,lag_a_nu);
		    
		    //build wvo block of W(2) 
		    //     (this is the wmat (NOT wmat(II) and NOT wmat(III))
		    wmat.block(no,0,nv,no) = lag_a_nu * occ_coeff_mat.transpose();
		    wmat.block(0,no,no,nv) = wmat.block(no,0,nv,no).transpose();
		    
		    //build wao block of W(2) (note the minus sign)
		    wmat.block(0,0,no,no) = -occ_coeff_mat * lag_mu_i;
		    
		    //I am not sure about the next bit
		    // it does not seems to change the result if it is on or off
		    // why did i have it? maybe to have wmat matrices match up between gamess and this?
		    //lets keep it around but 'off' for now
		    //	wmat.block(0,nf,nf,na) /= (double)2.0;
		    //  wmat.block(nf,0,na,nf) = wmat.block(0,nf,nf,na).transpose();
		    
		    //build wvv block of W(2) !this works for w vv-block
		    wmat.block(no,no,nv,nv) = lag_a_nu * virt_coeff_mat.transpose();
		    
		    
		    //build p(active,frozen) block of P(2) (pmat)
		    for(int i = nf; i < no; i++){
			for(int K = 0; K < nf; K++){ //note capital K
			    pmat(i,K) = occ_coeff_mat.row(K) * lag_mu_i.col(i);
			    pmat(i,K) /= ea[i]-ea[K];
			    pmat(K,i) = pmat(i,K);
			}//K
		    }//i
		    
		    //transform L_mu_i(1) -> L_ai(1) (AO->MO)
		    lag_mo = lag_a_nu*occ_coeff_mat.transpose();
		    
		    //transform L_a_nu(2) -> L_ai(2) (AO->MO)
		    lag_mo += virt_coeff_mat * lag_mu_i;
		    
		    
	
		    
		}//scope: recompute three-center ERI --> mixed Lag, W[I], P[core,active], MO lag
		
		
		delete[] ptr_gamma_ia_P;
		delete[] ptr_gamma_inu_P;
		
		delete[] ptr_lag_mu_i;
		delete[] ptr_lag_a_nu;
		
	    }//scope: build vvoo components of lagrangian
	    profile.recompute += timer;
	    
	    cout << "     Time for mixed Lagrangians L(occ,mu) and L(virtual,nu): "
		 << profile.recompute 
		 << std::endl<< std::endl;
	    
	    
	    //at this stage, Pia blocks are zero (this is referred to as the un-relaxed density)
	    // note: pmat itself has not been multiplied by two
	    double *ptr_zai = new double[(nv+no)*(nv+no)];
	    MapMatrixXd zai(ptr_zai,no+nv,no+nv);
	    memset(ptr_zai, 0, (no+nv)*(no+nv)*sizeof(double) );
	    
	    double *pFAO = new double[N*N];
	    MapMatrixXd FAO(pFAO,N,N);
	    memset(pFAO, 0, N*N*sizeof(double) );
	    
	    
	    {//scope: build vooo and vvvo contributions to the Lagrangian, compute z-vector (i.e. determine active/virtual block of pmat)
		
		
		BOOST_AUTO(const &ea, wf.e(wf.double_occ()));
		BOOST_AUTO(const &ev, wf.e(wf.single_virtual()));
		
		
		// The lagrangian matrix, needed for the zvector equations to get Zai (ie occ-virt density - Pai),
		//  can be constructed from terms that involve density terms (not the Pai part!) and terms that
		//  involve amplitudes. see equation (18) in the J.Comp.Chem vol.:28 page:839 (2007).
		//
		//       The density matrix terms involve the orbital hessian matrix (A_pqrs).
		//       The amplituded terms involved the (ab|Q) and (ij|Q) Eris and the 'gamma' matrix
		{//scope: Lagrangian formation / P(o,v) formation
		    
		    //build lagranigan (vvvo and vooo contributions) with J/K matrices
		    
		    cout << "Starting Lagrangian" 
			 << std::endl;
		    timer.reset();
		    ::cchem::rimp2_gradient::detail::matrix_factor 
			  lagrangian(N,no,nv,nl,
				     Ca.data().begin(), Cv.data().begin(),
				     LM1_MN_SYM,pmat,pe,pGPU);
		    
		    timerTest.reset();

		    //get pos/neg pseudosquare root of pmat
		    // pmat so far:
		    //   f a v
		    //  -------
		    // f|0 x 0|
		    // a|x x 0|
		    // v|0 0 x|   
		    //  -------
		    lagrangian.pseudoSquareRoot( pseudoSquareRootCutoff );
		    cout << "     Time for pseudoSquareRoot: " 
			 << timerTest << std::flush << std::endl;
		    
		    pe.barrier();
		    
		    timerTest.reset();

		    //build JK matrices
		    lagrangian.buildJK();

		    //for MI lagrangian
		    lagrangian.buildLagMO(lag_mo,FAO);

		    cout << "     Time to build MO Lagrangian: " << timerTest << std::endl;
		    
		    profile.lag_part += timer;
		    cout << "     Total time for VOOO and VVVO contributions to Lagrangian: "
			 << profile.lag_part << std::endl<< std::endl;
		    
		}//scope: Lagrangian formation / P(o,v) formation 

		
		/////////////////////////////////////////////////////
		//-----solve Z-Vector equations---------
		//   i.e. compute Pav block of pmat (relaxed density)
		/////////////////////////////////////////////////////

		pe.barrier();
		
		{//scope: zvector
		    int maxitc=50;          //this should be a control variable
		    double uconv = 1.0e-10; //this should be a control variable
		    double small = 1.0e-13; //this should be a control variable
		    size_t NMB = 500;       //this should be a control variable
		    
		    timer.reset();
		    
		    //solve Z-Vector equations with J/K matrices
		    ::cchem::rimp2_gradient::detail::JK_RIZVector 
			  zvector(N,nd,ns,nv,nl,maxitc,uconv,small,
				  boost::cref(ea),boost::cref(ev),
				  NMB,
				  boost::cref(Ca),boost::cref(Cv),
				  LM1_MN_SYM,
				  Ca.data().begin(),
				  pe,
				  pGPU
				  );

		    zvector.solve(lag_mo, pmat, zai);

		    profile.zvector += timer;
		    
		}//scope: zvector
		
		
		
		
		
		cout << std::endl << "     Time for solving the ZVector equations: " 
		     << profile.zvector 
		     << std::endl<< std::endl;
		
	    }//scope: build vooo and vvvo contributions to the Lagrangian, compute z-vector


	    timer.reset();
	    
	    //build remaining energy weighted density terms
	    //   wb_mat [ie w_baror W(II)] and wbb_mat [eg w_bar_bar or W(III)]
	    
	    pe.barrier();
	    
	    cout << "Starting remianing density terms" << std::endl;
	    
	    
	    {//scope: remaining density terms
		
		BOOST_AUTO(const &ea, wf.e(wf.double_occ()));
		BOOST_AUTO(const &ev, wf.e(wf.single_virtual()));
		
		// oo block of w_bar ( wmat[II] )
		for(int i = 0; i < no; i++){
		    for(int k = 0; k < no; k++){
			//wb_mat(i,k) = -0.5*pmat(i,k)*(ea[i] + ea[k]);
			wmat(i,k) -= 0.5*pmat(i,k)*(ea[i] + ea[k]);
		    }//k
		}//i
		
		// vv block of w_bar ( wmat[II] )
		for(int a = 0; a < nv; a++){
		    for(int c = 0; c < nv; c++){
			//wb_mat(a+no,c+no) = -0.5*pmat(a+no,c+no)*(ev[a]+ev[c]);
			wmat(a+no,c+no) -= 0.5*pmat(a+no,c+no)*(ev[a]+ev[c]);
		    }//c
		}//a
		
		// vo block of w_bar ( wmat[II] )
		for(int a = 0; a < nv; a++){
		    for(int i = 0; i < no; i++){
			//wb_mat(a+no,i) = -pmat(a+no,i)*ea[i];
			//wb_mat(i,a+no) = -pmat(a+no,i)*ea[i];
			wmat(a+no,i) -= pmat(a+no,i)*ea[i];
			wmat(i,a+no) -= pmat(a+no,i)*ea[i];
		    }//i
		}//a
		
		
		
		
		//use J/K matrices to compute wmat[III]
		zai *=(double)2;
		FAO *=(double)2;
		//FAO *=(double)-1;
		
		//build lagranigan (vvvo and vooo contributions) with J/K matrices
		//			omp_set_num_threads(4);
		timer.reset();
		//    size_t NMBtemp = 500;
		::cchem::rimp2_gradient::detail::matrix_factor 
		      wmat3(N,no,nv,nl,
			    Ca.data().begin(), Cv.data().begin(),
			    LM1_MN_SYM,zai,pe,pGPU);

		//if(pe.rank()==0)std::cout << zai << std::endl;
		//get pos/neg pseudosquare roots of zai
		wmat3.pseudoSquareRoot( pseudoSquareRootCutoff );

		//build JK matrices
		wmat3.buildJK();
		//form wmat[III]
		wmat3.buildWmat3( wmat, FAO );

		
	    }//scope: remaining density terms
		
	    //delete unneeded arrays
	    delete BARE_IA;
	    delete CIAQ;
	    delete GAMMA_IA_Q;
	    delete LM1_MN_SYM;
	    
	    pe.barrier();
	    
	    profile.remaining_density += timer;
	    cout << "     Time for remaining density terms: "
		 << profile.remaining_density << std::endl<< std::endl;
	    
	    delete[] pFAO;
	    delete[] ptr_zai;

	    //scale wmat by two
	    wmat *=(double)2;
	    
	    //  ADD wmat(SCF) TERM TO FORM wmat(MP2) IN THE MO BASIS
	    BOOST_AUTO(const &ea, wf.e(wf.double_occ()));
	    for(int i=0; i< no; i++) { wmat(i,i) -= (double)2*ea[i]; } //i


	    //scale pmat by two
	    pmat *=(double)2;

	    //build the relaxed density (in MO basis)
	    Eigen::MatrixXd Pmo(no+nv,no+nv);
	    Pmo = 2*pmat.block(0,0,no+nv,no+nv);
	    //Pmat is how the denisty changes from HF to relaxed 1PDM, 
	    //   so we need to add the Pscf denisty back in to get the real relaxed density
	    //    (which we know for RHF, 2 on diagonal of MO density matrix)
	    for(int i=0; i< no; i++) { Pmo(i,i) += double(2); }

	    
	    //back transform wmat from MO->AO basis
	    //set up btrans object
	    cchem::rimp2_gradient::detail::backtransform 
		btrans(N,
		       boost::cref(Ca),
		       boost::cref(Cv));
	    
	    //carry out back-transformation on energy weight density
	    //wmat is now in AO basis
	    btrans(wmat);

	    //carry out backtransformation on relaxed density
	    // pmat is now AO basis
	    btrans(pmat);
	    
	    //build P(SCF) in AO basis --> scf HF AO density
	    MapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() );
	    double *ptr_pscf = new double[N*N];
	    MapMatrixXd Pscf(ptr_pscf, N, N );
	    Pscf = 2.0*occ_coeff_mat.transpose()*occ_coeff_mat;

	    timer.reset();
	    
	    int natoms = molecule.size();
	    MapMatrixXd pqm1(pPqm1,nl,nl);
	    
	    ::cchem::rimp2_gradient::DFIntGradient::gradient dfintgradient
		  (boost::cref(Ca),boost::cref(Cv),
		   Pscf, Pmo, N, no, nv , nl, BARE_MN_SYM, pqm1,
		   boost::cref(auxbasis), boost::cref(basis), sph_c, pSphTrans,
		   egTerms.getMatrix(::cchem::rimp2_gradient::detail::gradientTerms::DFIntGradient),
		   natoms, spherical, pe, rt);
	    
	    dfintgradient.compute_DFIntGradient();
		
	    
	    profile.deri_4center +=timer;
	    cout << std::endl << "Done with DFIntegral gradients: " 
		 << profile.deri_4center
		 << std::endl << std::endl;

	    //adjust relaxed density to incorperate scf density (in AO basis)
	    pmat += Pscf;

	    delete [] ptr_pscf;
	    delete [] pPqm1;

	    cout <<std::setprecision(10);
	    cout << "Correlation Energy: " << e_term[0] << std::endl;
	    
	    E = std::accumulate(e_term,e_term+7,0.0);
	    
	    timer.reset();
	
	    {//scope: derivative overlap
		
		pe.task().reset(); 
		if(pe.rank() != 0) goto skip2;
		
		cout << "MAKE SURE NORMALIZATION IS CORRECT" << std::endl;
#pragma omp parallel
		if (pe.node().rank() == 0) {
		    
		    Eigen::MatrixXd eg(molecule.size(),3);
		    eg.setZero();
		    detail::Thread::Task<Parallel::Task&> task(pe.task());
#pragma omp barrier
		    
		    cchem::ri::AuxiliaryTwoCenterInt< ::rysq::TwoCenterDerivativeOverlap > 
			ds(boost::cref(basis));
		    
		    while (++task < shells.size()) {
			const Basis::Shell &S = basis.shells().at(task);
			const size_t s = task;
			
			for(size_t q= 0; q < task; ++q){
			    const Basis::Shell &Q = basis.shells().at(q);
			    
			    if( S.atom() == Q.atom() )continue;
			    
			    int satom = S.atom();
			    int qatom = Q.atom();
			    
			    //in order to compare to gamess integrals, you may have to mess with the
			    //normalization (right now it should be off)
			    //normalization taken care of inside of wmat
			    MapMatrixXd twoc_batch(ds(Q,S), Q.size()*S.size(),6);
			    
			    for(int qf = 0; qf<Q.size(); qf++){
				for(int sf = 0; sf<S.size(); sf++){
				    eg(satom,0) += 2*twoc_batch(sf+qf*S.size(),0)*wmat(Q.start()+qf,S.start()+sf);
				    eg(satom,1) += 2*twoc_batch(sf+qf*S.size(),1)*wmat(Q.start()+qf,S.start()+sf);
				    eg(satom,2) += 2*twoc_batch(sf+qf*S.size(),2)*wmat(Q.start()+qf,S.start()+sf);
				    
				    eg(qatom,0) -= 2*twoc_batch(sf+qf*S.size(),0)*wmat(Q.start()+qf,S.start()+sf);
				    eg(qatom,1) -= 2*twoc_batch(sf+qf*S.size(),1)*wmat(Q.start()+qf,S.start()+sf);
				    eg(qatom,2) -= 2*twoc_batch(sf+qf*S.size(),2)*wmat(Q.start()+qf,S.start()+sf);
				}//sf
			    }//qf
////-----debugging-----
//						std::cout << task << " " << q << std::endl;
//						for(int qf = 0; qf<Q.size(); qf++){
//							for(int sf = 0; sf<S.size(); sf++){
//								std::cout << "      "<< twoc_batch(sf+qf*S.size(),0) << " "
//										<< twoc_batch(sf+qf*S.size(),1) << " "
//										<< twoc_batch(sf+qf*S.size(),2) << std::endl;
//							}//sf
//						}//qf
////-----debugging-----
			} //q


		    } //task

#pragma omp critical
		    egTerms.getMatrix(::cchem::rimp2_gradient::detail::gradientTerms::densityForce ) += eg;
		    
		}//(pe.node().rank() == 0)
		
	    skip2:
		;
	    }//scope: derivative overlap
	    
	    profile.overlap_derivative += timer;
	    cout << "time for overlap derivative contributions: "
		 << profile.overlap_derivative << std::endl<< std::endl;
	    
	    timer.reset();
	    {//scope: do pmat terms (fock derivatve) with libint
		
		

#ifdef HAVE_LIBINT
		::libcchem_libint_interface
		    ::T_V_1body_deriv_contributions(N,ptr_pmat, basis, molecule ,
						    ptr_wmat,
						    egTerms.getMatrix(::cchem
								      ::rimp2_gradient
								      ::detail
								      ::gradientTerms::T),
						    egTerms.getMatrix(::cchem
								      ::rimp2_gradient
								      ::detail
								      ::gradientTerms::V)
						    );
#else
		cout << "WARNING: LIBINT is not present, "
		     << "T&V derivative contributions are not "
		     << "included in gradient." << std::endl;
#endif
		
	    }//scope: do pmat terms (fock derivatve) with libint
	    
	    profile.overlap_derivative += timer;
	    cout << "time for one-electron derivative contributions: "
		 << profile.overlap_derivative << std::endl<< std::endl;
	    


	    timer.reset();

	    {//scope: two-center derivative eri
		
		
		pe.task().reset(); 
		if(pe.rank() != 0) goto skip3;
		
		cout << "MAKE SURE NORMALIZATION IS CORRECT" << std::endl;
#pragma omp parallel
		if (pe.node().rank() == 0) {
		    
		    Eigen::MatrixXd eg(molecule.size(),3);
		    eg.setZero();
		    
		    detail::Thread::Task<Parallel::Task&> task(pe.task());
#pragma omp barrier
		    
		    cchem::ri::AuxiliaryTwoCenterInt< ::rysq::TwoCenterDerivativeEri > auxiliary_eri(boost::cref(auxbasis));
		    
		    while (++task < auxshells.size()) {
			const Basis::Shell &S = auxbasis.shells().at(task);
			const size_t s = task;
			const int satom = S.atom();
			
			for(size_t q= 0; q <= task; ++q){
			    const Basis::Shell &Q = auxbasis.shells().at(q);
			    const int qatom = Q.atom();
			    
			    if( satom == qatom )continue;
			    
			    //since we are computing the lower triangular, use translation invarience to compute
			    //   the corresponding ket derivatives


					//Q.shell() <= S.shell()     -->   < S | 1/r | Q >  !!! when sorted !!!

					// so S becomes the <bra| and Q is the |ket> -- we do root contraction over <bra| states
					//				std::cout << Q.size() << " " << S.size() << std::endl;

			    MapMatrixXd twoc_batch(auxiliary_eri(Q,S), Q.size()*S.size(),6);

			    //storage layout for twoc_batch e.g. ( "d-shell" | "P-shell" )
			    //
			    //   n.b.  dAi is the derivative with respect the ith coordinate of the bra center (NOT ket)
			    //         "A" means bra center
			    //
			    //        dAx  dAy  dAz  dBx  dBy  dBz
			    // (xx|:  |x)  |x)  |x)  |x)  |x)  |x)
			    // (yy|:  |x)  |x)  |x)  |x)  |x)  |x)
			    // (zz|:  |x)  |x)  |x)  |x)  |x)  |x)
			    // (xy|:  |x)  |x)  |x)  |x)  |x)  |x)
			    // (xz|:  |x)  |x)  |x)  |x)  |x)  |x)
			    // (yz|:  |x)  |x)  |x)  |x)  |x)  |x)
			    // (xx|:  |y)  |y)  |y)  |y)  |y)  |y)
			    // (yy|:  |y)  |y)  |y)  |y)  |y)  |y)
			    // (zz|:  |y)  |y)  |y)  |y)  |y)  |y)
			    // (xy|:  |y)  |y)  |y)  |y)  |y)  |y)
			    // (xz|:  |y)  |y)  |y)  |y)  |y)  |y)
			    // (yz|:  |y)  |y)  |y)  |y)  |y)  |y)
			    // (xx|:  |z)  |z)  |z)  |z)  |z)  |z)
			    // (yy|:  |z)  |z)  |z)  |z)  |z)  |z)
			    // (zz|:  |z)  |z)  |z)  |z)  |z)  |z)
			    // (xy|:  |z)  |z)  |z)  |z)  |z)  |z)
			    // (xz|:  |z)  |z)  |z)  |z)  |z)  |z)
			    // (yz|:  |z)  |z)  |z)  |z)  |z)  |z)
			    //
			    //
			    //  i.e. twoc_batch(0,0)  =  d/dAx ( dxx | px ) = 2*alpha*( fxxx | px ) - 2*( px | px )
			    //						  => d/dBx ( dxx | px ) = -d/dAx ( dxx | px )
			    //
			    //       twoc_batch(5,2)  =  d/dAy ( dyz | px ) = 2*alpha*( fyyz| px ) - ( pz | px )
			    //				          => d/dBy ( dyz | px ) = -d/dAy ( dyz | px )
			    //
			    
			    if(spherical){
				Eigen::MatrixXd transformed_batch(Q.sphsize(),S.sphsize());
				//    				transformed_batch = sph_c[S.Lmin()].block(0,0,S.size(),S.sphsize()).transpose()*twoc_batch*sph_c[Q.Lmin()].block(0,0,Q.size(),Q.sphsize());
				
				//transform over each coordinate (e.g. d/dAx)
				double *ptr_deriv = twoc_batch.data();
				
				
				for(int nder = 0; nder < 3; nder++){
				    MapMatrixXd coord_batch(ptr_deriv,S.size(),Q.size());
				    transformed_batch = sph_c[S.Lmin()].block(0,0,S.size(),S.sphsize()).transpose()*coord_batch*sph_c[Q.Lmin()].block(0,0,Q.size(),Q.sphsize());
				    ptr_deriv += Q.size()*S.size();
				    
				    
				    for(int st = 0; st < S.sphsize(); st++){
					for(int qt = 0; qt < Q.sphsize(); qt++){
					    eg(satom,nder) -= 2*gamma_rs2(Q.sphstart()+qt,S.sphstart()+st)*transformed_batch(st,qt);
					    eg(qatom,nder) += 2*gamma_rs2(Q.sphstart()+qt,S.sphstart()+st)*transformed_batch(st,qt);
					    
					}//q
				    }//s
				    
				}///ic
				
				
			    }else{ //(spherical)
				
				double *ptr_deriv = twoc_batch.data();
				
				
				for(int nder = 0; nder < 3; nder++){
				    MapMatrixXd cart_batch(ptr_deriv,S.size(),Q.size());
				    
				    
				    for(int st = 0; st < S.size(); st++){
					for(int qt = 0; qt < Q.size(); qt++){
					    eg(satom,nder) -= 2*gamma_rs2(Q.start()+qt,S.start()+st)*cart_batch(st,qt);
					    
					    eg(qatom,nder) += 2*gamma_rs2(Q.start()+qt,S.start()+st)*cart_batch(st,qt);
					    
					}//q
				    }//s
				    
				    ptr_deriv += Q.size()*S.size();
				    
				}///ic
				
			    }//(spherical)
			    
			} //q
			
		    } //task
		    
#pragma omp critical
		    egTerms.getMatrix(::cchem
				      ::rimp2_gradient
				      ::detail
				      ::gradientTerms::twoCenterNonSeparable ) += eg;
		    
		    
		}//(pe.node().rank() == 0)
		
	    skip3:;
		
	    } //scope: two-center derivative eri
	    profile.twoc_derivative += timer;
	    cout << "time for 2-e two-center derivative contributions: "
		 << profile.twoc_derivative << std::endl<< std::endl;
	    
	    
	    ////////////////////////
	    //clean up global memory
	    ////////////////////////
	    delete[] ptr_gamma_rs2;
	    delete[] ptr_wmat;
	    delete[] ptr_pmat;
	    delete[] ptr_lag_mo;
	    





	    //
	    //storage layout for threec_batch e.g. ( "d-shell" "p-shell" | "p shell" )
	    //
	    //   n.b.  "A" is the first  center ( "A" B  |  C )
	    //         "B" is the second center (  A "B" |  C )
	    //         "C" is the third  center (  A  B  | "C" )
	    //         d/dAi is the derivative with respect the ith coordinate of the first  bra center (NOT ket)
	    //         d/dBi is the derivative with respect the ith coordinate of the second bra center (NOT ket)
	    //         d/dCi is the derivative with respect the ith coordinate of the ket center (NOT bra)
	    //		   alpha is the gaussian exponent on center A
	    //		   beta is the gaussian exponent on center B
	    //
	    //
	    //         dAx  dAy  dAz  dBx  dBy  dBz  dCx  dCy  dCz
	    // (xx x|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (yy x|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (zz x|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xy x|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xz x|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (yz x|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xx y|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (yy y|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (zz y|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xy y|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xz y|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (yz y|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xx z|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (yy z|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (zz z|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xy z|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xz z|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (yz z|: |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)  |x)
	    // (xx x|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (yy x|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (zz x|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xy x|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xz x|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (yz x|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xx y|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (yy y|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (zz y|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xy y|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xz y|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (yz y|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xx z|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (yy z|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (zz z|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xy z|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xz z|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (yz z|: |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)  |y)
	    // (xx x|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (yy x|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (zz x|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xy x|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xz x|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (yz x|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xx y|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (yy y|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (zz y|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xy y|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xz y|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (yz y|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xx z|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (yy z|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (zz z|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xy z|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (xz z|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    // (yz z|: |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)  |z)
	    //
	    //  i.e. threec_batch(0,0)  =  d/dAx (xx x|x) = 2*alpha*(xxx x|x) - 2*(x x|x)
	    //                          => d/dBx (xx x|x) = 2*beta*(xx xx|x) - (xx s|x)
	    //                          => d/dCx = -d/dAx -d/dBx
	    //
	    
	    //auxbasis.max().size() : cartesian
	    //max_auxbasis_shell_size : cartesian/spherical (overide)
	    
	    timer.reset();
	    //	omp_set_num_threads(4);
	    
	    { //scope: three-center derivative eri





		pe.task().reset(); 
		
		Eigen::MatrixXd egNode(molecule.size(),3);
		egNode.setZero();
		

#pragma omp parallel
		if (pe.node().rank() == 0) { //i think rank is the process in the global space
		    
		    Eigen::MatrixXd eg(molecule.size(),3);
		    eg.setZero();
		    detail::Thread::Task<Parallel::Task&> task(pe.task());
		    
		    //storage for GAMMA_inu^P : '1/2' back transformed
		    double *ptr_gamma_inu_P = new double[max_auxbasis_shell_size*N*(no-nf)];
		    
		    //			std::vector<MapMatrixXd>gamma_inu_P_vec;
		    double *ptr_temp = ptr_gamma_inu_P;
		    std::vector<double *>pGammaINuPVec;
		    for(int l = 0; l < max_auxbasis_shell_size; l++){
			//				MapMatrixXd temp(ptr_temp,N, (no-nf) );
			//				gamma_inu_P_vec.push_back(temp);
			pGammaINuPVec.push_back(ptr_temp);
			ptr_temp += N*(no-nf);
		    }//l
		    
		    
		    
		    //storage for GAMMA_munu^P : '2/2' back transformed
		    double *ptr_gamma_munu_P = new double[max_auxbasis_shell_size*N*N];
		    
		    //			std::vector<MapMatrixXd>gamma_munu_P_vec;
		    double *ptr_temp2 = ptr_gamma_munu_P;
		    std::vector<double *> pGammaMuNuPVec;
		    for(int l = 0; l < max_auxbasis_shell_size; l++){
			//				MapMatrixXd temp(ptr_temp2,N,N);
			//				gamma_munu_P_vec.push_back(temp);
			pGammaMuNuPVec.push_back(ptr_temp2);
			ptr_temp2 += N*N;
		    }//l
		    
		    
		    
		    const size_t MaxObsShell = basis.max().size();    // #functions in largest obs shell
		    const size_t MaxAuxShell = auxbasis.max().size(); // #functions in largest aux shell
		    double *pTransBatch = NULL;
		    if(spherical)pTransBatch = new double[MaxObsShell*MaxObsShell*MaxAuxShell];
		    
#pragma omp barrier
		    
		    //			MapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() ); //(no,N)
		    
		    cchem::ri::AuxiliaryThreeCenterInt < ::rysq::ThreeCenterDerivativeEri > 
			auxiliary_eri(
				      boost::cref(auxbasis),
				      boost::cref(basis),
				      boost::cref(Ca),
				      boost::cref(Cv) 
				      );
		    
		    while (++task < auxbasis.shells().size() ) { //this is the L shell
			const Basis::Shell &L = auxbasis.shells().at(task);
			const int latom = L.atom();
			
			const size_t LSize = L.size();
			const size_t LStart = L.start();
			
			size_t gammaLSize = LSize;
			size_t gammaLStart = LStart;
			if(spherical){
			    gammaLSize = L.sphsize();
			    gammaLStart = L.sphstart();
			}
			for(int l = 0; l < gammaLSize; l++){
			    double *pGammaMuNuP = pGammaMuNuPVec[l];
			    double *pGammaINuP = pGammaINuPVec[l];
			    
			    size_t start[] ={ nf*N, gammaLStart +l };
			    size_t finish[] ={ no*N, gammaLStart +l +1 };
			    GAMMA_INU_P->get(pGammaINuP,start,finish);
			    
			    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
					N, N, na,
					1.0, pGammaINuP, N,
					&Ca.data().begin()[nf], no,
					0.0, pGammaMuNuP, N );
			    
	 				//symmetrize result
			    for(int i = 0; i < N; i++){
				for(int j = 0; j < i; j++){
				    pGammaMuNuP[i +j*N] += pGammaMuNuP[j +i*N];
				    pGammaMuNuP[j +i*N] = pGammaMuNuP[i +j*N];
				} //j
				pGammaMuNuP[i +i*N] *= 2;
			    } //i
			    
			    cblas_dscal(N*N, 0.5,  pGammaMuNuP, 1);
			    
			}//l
			
			
			for(size_t s = 0; s < basis.shells().size(); s++){
			    const Basis::Shell &S = basis.shells().at(s);
			    const int satom = S.atom();
			    const size_t SStart = S.start();
			    const size_t SSize = S.size();
			    
			    for (size_t q = 0; q <= s; ++q) {
				const Basis::Shell &Q = basis.shells().at(q);
				const int qatom = Q.atom();
				const size_t QStart = Q.start();
				const size_t QSize = Q.size();
				const size_t SQSize = SSize*QSize;
				
				
				if(latom == satom && latom == qatom)continue;
				
				double perm = 4.0;
				if( q!=s ) perm = 8.0;
				
				
				//when sorted: s >= q
				MapMatrixXd threec_batch( auxiliary_eri(Q,S,L), Q.size()*S.size()*L.size(), 9);
				
				double * ptr_deriv = threec_batch.data();
				
				if(spherical){
				    
				    const size_t LSphSize = L.sphsize();
				    const size_t LSphStart = L.sphstart();
				    double *pLmin = pSphTrans[L.Lmin()];
				    
				    //if something is not working, do we have sorted shells???
				    for(int nder = 0; nder < 9; nder++){
					int center = L.atom();
					if(nder < 6)center = Q.atom();
					if(nder < 3)center = S.atom();
					
					const int ic = (nder % 3);
					
					//MapMatrixXd deriv_batch(ptr_deriv,S.size()*Q.size(),L.size());
					
					//Eigen::MatrixXd transformed_batch( S.size()*Q.size(),L.sphsize() );
					//transformed_batch = deriv_batch*sph_c[L.Lmin()].block(0,0,L.size(),L.sphsize());
					//ptr_deriv += Q.size()*S.size()*L.size();
					
					cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
						    SQSize, LSphSize, LSize,
						    1.0, ptr_deriv, SQSize,
						    pLmin, LSize,
						    0.0, pTransBatch, SQSize );
					
					
					double val = 0.0;
					for (int lt = 0; lt < LSphSize; lt++) {
					    
					    const double *pInt = &pTransBatch[lt*SQSize];
					    
					    for(int st = SStart, qst=0; st < SStart+SSize; st++){
						const double *pGammaMuNuP = &pGammaMuNuPVec[lt][st*N];
						
						for(int qt = QStart; qt < QStart+QSize; qt++,qst++){

						    //eg(center,ic) += perm*transformed_batch(qst ,lt )*gamma_munu_P(qt,st);
						    //eg(center,ic) += perm*pInt[qst]*gamma_munu_P(qt,st);
						    //eg(center,ic) += perm*pInt[qst]*pGammaMuNuP[ qt ]; // gamma_munu_P(qt,st);
						    val += pInt[qst]*pGammaMuNuP[ qt ]; // gamma_munu_P(qt,st);

						}//qt
					    }//st
					}//lt
					
					eg(center,ic) += perm*val;
					ptr_deriv += SQSize*LSize;
					
				    }//nder
				    
				}else{ //(spherical)
				    
				    
				    for(int nder = 0; nder < 9; nder++){
					int center = L.atom();
					if(nder < 6)center = Q.atom();
					if(nder < 3)center = S.atom();
					
					const int ic = (nder % 3);
					
					double val = 0.0;
					for (int lt = 0; lt < L.size(); lt++) {
					    
					    const double *pInt = &ptr_deriv[lt*SQSize];
					    //									MapMatrixXd & gamma_munu_P = gamma_munu_P_vec.at(lt);
					    
					    for(int st = SStart, qst = 0; st < SStart+SSize; st++){
						const double *pGammaMuNuP = &pGammaMuNuPVec[lt][st*N];
						
						for(int qt = QStart; qt < QStart+QSize; qt++,qst++){
						    
						    //eg(center,ic) += perm*threec_batch(qsl ,nder )*gamma_munu_P(qt,st);
						    val += pInt[qst]*pGammaMuNuP[ qt ]; //*gamma_munu_P(qt,st);
						    
						}//st
					    }//qt
					}//lt
					
					eg(center,ic) += perm*val;
					ptr_deriv += SQSize*LSize;
					
				    }//nder
				    
				} //(spherical)
				
				
			    }// q
			    
			}// s
			
		    }//++task : auxiliary shells
		    
		    
#pragma omp critical
		    egNode += eg;
		    
		    delete [] pTransBatch;
		    
		    delete [] ptr_gamma_munu_P;
		    delete [] ptr_gamma_inu_P;

		}//pe.node().rank() == 0
		
		
		pe.reduce("+",egNode.data(), (size_t)(molecule.size()*3)); 
		
		egTerms.getMatrix(::cchem::rimp2_gradient::detail::gradientTerms::threeCenterNonSeparable ) += egNode;
		
		
	    } //scope: three-center derivative eri
	    profile.threec_derivative += timer;
	    cout << "time for 2-e three-center derivative contributions: "
		 << profile.threec_derivative << std::endl<< std::endl;
	    
	    //compute classical nuclear-nuclear gradient
	    ::cchem::rimp2_gradient::detail::nuclear_deriv(molecule, egTerms.getPointer(::cchem::rimp2_gradient::detail::gradientTerms::nuclearForce));
	    
	    if(pe.rank() == 0)egTerms.printGradientTerms();
	    
	    egTerms.gradientSum(egGlobal);
	    
	    if(gradientDebug && pe.rank() == 0){
		egTerms.printGradientTerms();
	    }//(egGlobal)
	    
	    cout << std::endl << "total gradient" << std::endl;
	    cout << "if this looks really wrong, make sure you have " <<
		"renorm() commented out in libint" << std::endl;
	    cout << egGlobal << std::endl;

	}//scope: gradient work starts here

	used = memory.used();
	memory.clear();
	
	{
	    BOOST_AUTO(cout, rt.cout());
	    
	    cout << "    memory: " << used/(1<<20) << " MB" << std::endl;
	    cout << std::endl;
	    BOOST_PROFILE_DUMP(cout);
	}
	
	
	//print out energy decomposition for ZAPT (i.e. ns > 0)
	if(ns > 0 && pe.rank()==0){
	    cout << std::fixed <<std::setprecision(10) << std::showpos;
	    cout <<"   RI-ZAPT ENERGY CONTRIBUTIONS" << std::endl;
	    cout <<"       E1 = "<< e_term[0] <<"   CLOSED SHELL-LIKE TERM" << std::endl;
	    cout <<"       E2 = "<< e_term[1] <<"   SINGLY UNOCCUPIED"<< std::endl;
	    cout <<"       E3 = "<< e_term[2] <<"   2 SINGLY OCCUPIED"<< std::endl;
	    cout <<"       E4 = "<< e_term[3] <<"   SINGLY UNOCCUPIED/OCCUPIED"<< std::endl;
	    cout <<"       E5 = "<< e_term[4] <<"   \"FOCK MATRIX CONTRIBUTION\""<< std::endl;
	    cout <<"       E6 = "<< e_term[5] <<"   SINGLY OCCUPIED"<< std::endl;
	    cout <<"       E7 = "<< e_term[6] <<"   2 SINGLY UNOCCUPIED"<< std::endl;
	}

	//delete unneeded arrays
	delete GAMMA_INU_P;
	delete BARE_MN_SYM;

	pe.barrier();

	if(pe.rank() == 0)std::cout << "global rimp2 gradient time: "<< globalTimer << std::endl;
	
	return e_term[0];

}

}//namespace cchem



