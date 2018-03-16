/*
 * ri-energy.cpp
 *
 *  Created on: May 7, 2015
 *      Author: luke
 */



#include "ri-mp2/rimp2.hpp"

#include "core/wavefunction.hpp"


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
//#include <cuda_profiler_api.h>
//#include "/share/apps/cuda7/include/cuda_profiler_api.h"
#if HAVE_CUBLAS
#include <cuda_runtime.h>
#include <cuda.h>
#endif



//#pragma GCC system_header
#include "array/hdf5.hpp"
#include "utility/timer.hpp"

#include <algorithm>
#include <iostream>
#include <memory>


#include <boost/noncopyable.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/thread.hpp>

#include "boost/utility/profiler.hpp"
#include <numeric>

#include <fstream>
#include <algorithm>
#include <stdlib.h> //posix_memalign


#include <Eigen/Dense>


#include <ri-energy-terms.hpp>
#include <ri-openshell-work.hpp>


#include "ri-integrals.hpp"
#include "ri-async.hpp"


#include <ctime>



#include <three-center-work.hpp>
#include <ERIContractionEnergyAccum.hpp>

 //lbr
//#undef HAVE_CUBLAS
//#define HAVE_CUBLAS 1

#include <metric.hpp>

//#if HAVE_CUBLAS
#include <device.hpp>
//#endif
#include <math.hpp>
namespace cchem {

    namespace rimp2 {


	typedef boost::numeric::ublas::matrix<
	    double, boost::numeric::ublas::column_major> Matrix;
	typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
	typedef ::rysq::TwoCenterEri TWO_ERI;
	typedef ::rysq::ThreeCenterEri THREE_ERI;

	namespace detail {

	    using cchem::Thread;
	    
	}//detail
    }//rimp2



/*
   @brief ri-mp2/ri-zapt energy evaluation
   @author LBR
   @detail ri-mp2/ri-zapt energy evaluation
   @param wf wavefunction
   @param rt runtime
*/
double rimp2::energy(Wavefunction wf, Runtime &rt) {


    bool debug = 0;

    BOOST_AUTO(const &basis, wf.basis());
    BOOST_AUTO(const &shells, basis.shells());

    BOOST_AUTO(const &auxbasis, wf.auxbasis());
    BOOST_AUTO(const &auxshells, auxbasis.shells());

    // The wf.sort must be on. This deals with the transfer equations in
    //    the rysq quadrature code. You need to have shell pairs in which center
    //    A has a higher moment than center B.  (A B| C)
    //    This is becuase momentum is (see transfer equations)
    //    transfered to center B from A (A --> B). The way in which the code is
    //    currently setup (ThreeCenterWork::getThreeCenterEriBatch)
    //    does not reorder the pairs, so the sort is needed. Otherwise,
    //    the rysq quadrature code will need to be rewritten to transfer momentum
    //    from center B to center A (not worth it)
    //
    wf.sort();
    wf.auxsort();


    double E = 0;
    double e_term[7] = {0};

    const double pseudoTol = rt.get<double>("/rimp2/inversion/cutoff", 1e-10);

    double cutoff = rt.get<double>("/mp2/integrals/cutoff", 1e-10);
    integrals::Screening screening(basis, cutoff);

    const std::string couMeth( rt.get<std::string>("/rimp2/coumeth", "eigen") );
    // std::cout << couMeth << " " << couMeth.size() << std::endl;

    std::cout << std::endl;
    if( !strncmp(couMeth.c_str(),"eigen",5 ) ) std::cout << "Using Eigen Decompostion "; 
    if( !strncmp(couMeth.c_str(),"cholesky",8 ) ) std::cout << "Using Cholesky Decomposition ";
    std::cout << "of Coulomb Metric"<< std::endl;


    //active orbtials coefficients (doubly occupied and singly occupied)
    Matrix Ca = trans(wf.C(wf.active()));

    //size of orbitals basis set
    size_t N = basis.size();


    //number of active orbitals (doubly and singly occupied)
    size_t no = wf.active().size();
    size_t nv = wf.virtuals().size();

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

    //parallel environment
    Parallel pe;

    // suppress output
    if (pe.rank() != 0) rt.set_cout(0);

	BOOST_AUTO(cout, rt.cout());

    //number of GPUs, will reset below
	//int ndevice=0;

    Runtime::Memory &memory = rt.memory();
    double used = 0;

    size_t max_threads = omp_get_max_threads();

    size_t nwords = (memory.available()/sizeof(double));



    int spherical = wf.do_spherical();

//    spherical = 0; //override

    //size of auxiliary basis set
    size_t nl;
    if(spherical){
    	nl = auxbasis.spherical_size();
    }else{
    	nl = auxbasis.size();
    }

    {
	cout << std::endl << "active:     " << no << std::endl;
	cout << "virtual:    " << nv << std::endl;
	cout << "atomic:     " << N << std::endl;
	cout << "auxiliary:  " << nl << std::endl;
	cout << "eri cutoff: " << screening.value() << std::endl;
	cout << "pseudo inversion tolerance : "	<< pseudoTol << std::endl << std::endl;
    }

    std::vector<MapMatrixXd> sph_c;
    if(spherical){
    	cout << std::endl<< "spherical auxiliary basis set will be used to construct (L|M)" << std::endl<< std::endl;
    for(int lang = 0; lang<7; lang++){
    	int dim = (lang+1)*(lang+2)/2;
    	MapMatrixXd temp( wf.spherical_coeff_ptr(lang),dim,dim);
    	sph_c.push_back(temp);
//    	if(lang>1)std::cout << sph_c[2] << std::endl;
    }
//    std::cout << sph_c[2] << std::endl;

    }



    // runtime may be different, ensure consistency
    pe.broadcast(&max_threads, 1, 0);
    pe.broadcast(&nwords, 1, 0);

    {

    	BOOST_AUTO(const &mem_threec, cchem::ri::AuxiliaryThreeCenterInt < THREE_ERI> ::memory(auxbasis, basis, Ca, Cv));
    	BOOST_AUTO(const &mem_twoc, cchem::ri::AuxiliaryTwoCenterInt< TWO_ERI >::memory(auxbasis));

    	size_t required = 0;
    	foreach (size_t n, mem_twoc) required += n;
    	required *= max_threads;
    	required +=nl*nl; //for L^-1 - only one copy.
    	if(spherical) {
    		required += auxbasis.max().size()*auxbasis.max().size()*max_threads;
    	}
    	if (required > nwords)
    		throw cchem::exception("not enough memory to run RI-MP2");

    	required = 0;
    	foreach (size_t n, mem_threec) required += n;
    	required *= max_threads;
    	required +=nl*nl; //for L^-1 - carried over from two-center Eri
    	if(spherical){
    		required += auxbasis.max_spherical_shell_size()*basis.size()*basis.size()*max_threads;
    	}
    	if (required > nwords)
    		throw cchem::exception("not enough memory to run RI-MP2");

    	nwords -= required;

    }


    struct {

	utility::timer::value_type global, coulombMetric, threeCTrans, 
	    coeffTrans, ERIAccum, readTime, writeTime, inversion, 
	    t2Trans, t2Write, t3Trans,  t3Write, t2Readt3Write,
	    twoCenter,energyAccum, eriFormation,eriFormationii,
	    eriFormationij,eriReads,eri,threeCenter,t2Read,cholesky,
	    gpuSync33,gpuSyncEri,gpuSyncTrans;
	
    } profile = {};


    utility::timer globalTimer, stepTimer, fineTimer,fineTimer2;

    globalTimer.reset();








    



/*
  Set up memory for cpu (and gpu)
      CPU: just statically assign it (e.g. ngb)
      GPU/CPU: data buffers are based on the size of device memory
*/


    //this is static (probably should be input variable?)
    size_t ngb1 = 5; //20;
    //amount of words used form 3 center ERIs and 1/3+2/3 transformation
    const size_t eriTransWords = ngb1*1e9/sizeof(double);

    int stride;
    int strideRow = 1;

    size_t freeWords;
    size_t sizePerStride;
    int maxStride;
    int numReads;
 


#if !HAVE_CUBLAS

    const int ndevice = 0;
    size_t ngb = 20;
    freeWords = ngb*1e9/sizeof(double);


    //reduce for coulomb metric
    freeWords -= nl*nl;
    
    sizePerStride = 4*nl*nvs + omp_get_max_threads()*nvs*nvs;
    
    maxStride = freeWords / sizePerStride;
    numReads = no/maxStride + (no%maxStride > 0);
    stride = no/numReads + (no%maxStride > 0);
    stride = maxStride -1;
    if(maxStride > no/2)stride = no/2+1;

#endif


    ::cchem::rimp2::device * pGPU = NULL;

#if HAVE_CUBLAS
    
    //blocks per cuda kernel
    const int cudaBlocks = 8; //do not mess with yet
    //threads per block
    const int cudaThreadsPB = 64; //do not mess with yet
    
    const int nstreams = ::cchem::rimp2::device::getNumberOfStreams();
    const int ndevice = ::cchem::rimp2::device::getNumberOfDevices();

    ::cchem::rimp2::device GPU(ndevice,nstreams, cudaBlocks, cudaThreadsPB, (pe.rank()==0) );
    pGPU = &GPU;

    //available words per GPU for stride
    freeWords = pGPU->maxWords();

    
    //storage of off-diagonal block of B_ia^Q integrals (1 per stream)
    freeWords -= nl*nvs*nstreams;
    //storage for energy contributions (512 is hardwired value 64*8)
    freeWords -= 512*nstreams;
    //occupied orbital energies
    freeWords -= no;
    //singly occupied and virtual orbital energies
    freeWords -= nvs;
    //for Lm1
    freeWords -= nl*nl;
    
    //something of zapt
    if(ns){freeWords -= nd*nv*nstreams;};
    
    
    //cost to store each diagonal block (i.e. stride) on device
    sizePerStride = nl*nvs + nvs*nvs*nstreams;
    
    
    maxStride = freeWords / sizePerStride;
    numReads = no/maxStride + (no%maxStride > 0);
    stride = no/numReads + (no%maxStride > 0);
    stride = maxStride -1;
    if(maxStride > no/2)stride = no/2+1;
    
    size_t strideTemp = stride;
    //in case of heterogenous hardware, let other nodes adjusted to
    //   the smaller stride
    for(size_t ip = 0; ip < pe.size(); ip++){
	
	if(ip == pe.rank()){
	    pe.broadcast(&stride, 1, ip);
	}else{
	    pe.broadcast(&strideTemp, 1, ip);
	}
	if(strideTemp < stride)stride = strideTemp;
    }//ip
    
#endif
    
    if(debug){
	cout << std::endl << "stride size: " << stride
	     << " occupied: " << no << std::endl;;
    }


















/*
       Coulomb Metric Work

 1) build two-center ERIs
 2) decompose ( w/ Cholesky ) two-center ERI matrix L
 3a) invert L ( i.e. L^-1 )  (cpu code)
 3b) invert L and form (P|Q)-1 = Lm1*Lm1^T  (if you have cublas)
*/
    stepTimer.reset();
    double *data_twoc = new double[nl*nl];
    double *data_pqm1 = new double[nl*nl];

    {  //scope: two-center eri

    	fineTimer.reset();
    	MapMatrixXd twoc_eri(data_twoc, nl, nl);

    	::cchem::rimp2_gradient::detail::Metric
		 metric(boost::cref(auxbasis), pe, spherical, nl, data_twoc, data_pqm1);

    	//build (A|B) matrix
    	metric.BuildTwoCenterEri(sph_c);

    	//form L-1 metric
    	metric.choleskyInversion(data_pqm1, data_twoc, (pe.rank() == 0) );


	// //form (P|Q)^(-1/2)
	// double pseudoTol=1e-10;
	// metric.pseudoInversion(pseudoTol, data_pqm1, data_twoc, (pe.rank()==0) );
	// //metric.pseudoInversion(pseudoTol, pLm1, pTwoC, (pe.rank()==0) );
    	

	// memcpy(data_twoc, data_pqm1, nl*nl*sizeof(double) );
	// //form (PQ)^-1 from (P|Q)^-0.5
	// //  nb: A = (P|Q)^0.5 * (Q|P)^0.5
	// //      A^-1 = (Q|P)^-0.5 * (P|Q)^-0.5  (nb swapped transpose)
	// cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
	// 	    nl,nl,nl,
	// 	    1.0, data_twoc, nl,
	// 	    data_twoc, nl,
	// 	    0.0, data_pqm1, nl);


	profile.inversion += fineTimer;


    	// Eigen::MatrixXd I(nl,nl);
    	// I = pqm1*A;
    	// std::cout << I.block(0,0,5,5) << std::endl;

    	// double Anorm =A.lpNorm<Eigen::Infinity>();
    	// double Am1Norm = Am1.lpNorm<Eigen::Infinity>();

    	// std::cout << "log10 of coulomb metric condition number: "
    	// 	  << std::log10(Anorm*Am1Norm) << std::endl << std::endl;

    } //scope: two-center eri

    profile.coulombMetric += stepTimer;

    cout << "time for two-center eri            "
    		<< profile.twoCenter << std::endl;
    cout << "time for Cholesky Decomposition    "
    		<< profile.cholesky <<std::endl;
    cout << "time to form   L^-1                "
    		<< profile.inversion <<std::endl<< std::endl;
    cout << "total time to form Coulomb Metric  "
    		<< profile.coulombMetric <<std::endl<< std::endl;



//set up array to store 2/3 index transformed integrals
{

	size_t dims[2] = {  no*nvs, nl }; 
	size_t chunk[2] = {  nvs, nl };
	rt.arrays().allocate<double,2>("rimp2.vp(b,qs,l)", dims, pe, chunk);

}

cout << "creating arrays: " << std::endl;
Array<double> *Vp = rt.arrays().find< Array<double> >("rimp2.vp(b,qs,l)");
rt.cout() << "   " << *Vp << std::endl<< std::endl;




/*
 this next part determines working storage for 3-center ERIs and transformation
    and allocates memory on host (and GPU if applicable)
*/

	//blockRanges keeps track of how the integrals are chunked over auxiliary basis shells
	std::vector< size_t > blockRanges;
	//start at zero
	blockRanges.push_back( 0 );

	//max domains refer to the number of integrals (AO and MO) in a domain block
	size_t maxDomainSizeAO = 0;
	size_t domainSizeAO = 0;
	size_t maxDomainSizeMO = 0;
	size_t domainSizeMO = 0;

	//keep track of how much memory is left
	size_t memLeft = eriTransWords;

	for (size_t l = 0; l < auxbasis.shells().size(); l++){
		const Basis::Shell &L = auxbasis.shells().at(l);

		size_t tripletSizeAO = L.size()*N*N;

		size_t tripletSizeMO = L.size()*no*nvs;
		if(spherical) tripletSizeMO =  L.sphsize()*no*nvs;

		//how much memory will is need in integrals of l auxiliary shell
		//the '2' refers to the buffers (cpu: one for work, one for IO)
		size_t memChunk = 2 * (tripletSizeAO + tripletSizeMO);
		//		//(cpu: one for work, one for IO)
		//		//(gpu: one for integrals, one for transformed integrals, one for IO)
		//		size_t memChunk = nBuffers * (tripletSizeAO + tripletSizeMO);

		if( memChunk < memLeft){

			memLeft -= memChunk;

			domainSizeAO += L.size();
			domainSizeMO += L.size();
			if(spherical)domainSizeMO += L.sphsize() - L.size();

			maxDomainSizeAO = std::max(maxDomainSizeAO, domainSizeAO);
			maxDomainSizeMO = std::max(maxDomainSizeMO, domainSizeMO);

		}else{

			maxDomainSizeAO = std::max(maxDomainSizeAO, domainSizeAO);
			maxDomainSizeMO = std::max(maxDomainSizeMO, domainSizeMO);
			domainSizeAO = 0;
			domainSizeMO = 0;
			blockRanges.push_back( l );
			memLeft = eriTransWords - memChunk;

			domainSizeAO += L.size();
			domainSizeMO += L.size();
			if(spherical)domainSizeMO += L.sphsize() - L.size();

		}//if( tripletSize < eriTransWords)

	}
	//always push the total number of shell onto blockRanges
	blockRanges.push_back( auxbasis.shells().size() );


	if(debug){
	    cout << "maxDomainSizeAO " << maxDomainSizeAO << std::endl;
	    cout << "maxDomainSizeMO " << maxDomainSizeMO << std::endl;
	    for (int iblock = 0; iblock < blockRanges.size(); iblock++)
		cout << "block " << iblock 
		     << " first aux shell id " << blockRanges[iblock] << std::endl;
	    cout << "NOTE: the last block is only a boundary "<<
		"(not to be computed)" << std::endl;
	}








	//pointers to work/IO buffers
	double *ptr_buff1;
	double *ptr_buff2;
	double *ptr_buff3;
	double *ptr_buff4;

#if !HAVE_CUBLAS

	ptr_buff1 = new double[ maxDomainSizeAO*N*N ];
	ptr_buff2 = new double[ maxDomainSizeAO*N*N ];
	ptr_buff3 = new double[ maxDomainSizeMO*no*nvs ];
	ptr_buff4 = new double[ maxDomainSizeMO*no*nvs ];

#elif HAVE_CUBLAS

	//page-locked memory
	GPUerrchk(cudaMallocHost( &ptr_buff1, (size_t)(maxDomainSizeAO*N*N*sizeof(double)) ));
	GPUerrchk(cudaMallocHost( &ptr_buff2, (size_t)(maxDomainSizeAO*N*N*sizeof(double)) ));
	GPUerrchk(cudaMallocHost( &ptr_buff3, (size_t)(maxDomainSizeMO*no*nvs*sizeof(double)) ));
	GPUerrchk(cudaMallocHost( &ptr_buff4, (size_t)(maxDomainSizeMO*no*nvs*sizeof(double)) ));

#endif




// 1) evaluate three-center ERIs
// 2) perform 1/3 and 2/3 AO to MO transformations
stepTimer.reset();
utility::Progress progress;

cout << std::endl
<< "Working on 3-center ERI and 1/3 + 2/3 transformation" << std::endl;
if (pe.rank() == 0) progress.reset(blockRanges.size()-1 );

{   //scope: three-center Eri and transform


    cchem::ri::async async;

    std::vector<size_t> haveBlock;
    std::vector<size_t>::iterator blockIter;


#if HAVE_CUBLAS
    //GPU transform object
    // cchem::threeCenterInt::threeCenterTransform 
    // 	MOBatch(boost::cref(Ca),
    // 		boost::cref(Cv),
    // 		boost::cref(auxbasis), 
    // 		spherical, blockRanges, nvs, no, N,
    // 		ndevice,nstreams,dv1_handle,dv1_streams,
    // 		async,Vp);
    cchem::threeCenterInt::threeCenterTransform 
	MOBatch(boost::cref(Ca),
		boost::cref(Cv),
		boost::cref(auxbasis), 
		spherical, blockRanges, nvs, no, N,
		ndevice,nstreams,pGPU->cublasHandles(),pGPU->cudaStreams(),
		async,Vp);

#endif



    pe.task().reset();
    detail::Thread::Task<Parallel::Task&> task(pe.task());
    task.reset();

    //this checks to if async.wait() should be called
    bool activeWrite = 0;

    //if(pe.rank() != 0)goto bolt;
    //{
#pragma omp parallel
    //    if (pe.node().rank() == 0 && pe.rank() == 0) {
   if (pe.node().rank() == 0) {
	   
	//this memory is for spherical to cartesian transformation
	double *threec_cart_to_sph = NULL;
	if(spherical) threec_cart_to_sph = 
			  new double[auxbasis.max_spherical_shell_size()*N*N];
	
	//eri object
	cchem::ri::AuxiliaryThreeCenterInt<THREE_ERI> 
	    auxiliary_eri(boost::cref(auxbasis),
			  boost::cref(basis),
			  boost::cref(Ca),
			  boost::cref(Cv) );
	
	//this gets batch of three ceneter ERIs.
	cchem::threeCenterInt::threeCenterEriWork 
	    AOBatch(boost::cref(basis),
		    boost::cref(auxbasis),
		    spherical,threec_cart_to_sph,N);
	
	
#if !HAVE_CUBLAS
	//CPU transform object
	cchem::threeCenterInt::threeCenterTransform 
	    MOBatch(boost::cref(Ca),
		    boost::cref(Cv),
		    boost::cref(auxbasis), 
		    spherical, blockRanges, nvs, no, N,async,Vp);
#endif //!HAVE_CUBLAS
	

//this has to occur for the counter to work correctly (the way I am using it)
#pragma omp master
	{
	    ++task;
	    if(debug)
		std::cout << pe.rank() << " " 
			  << pe.node().rank() 
			  << " working on task " << task 
			  << std::endl <<std::flush;
	}

#pragma omp barrier

	//loop over domains
	for(size_t domain = 0; domain < blockRanges.size()-1; domain++){
	    size_t shellIndex = blockRanges[domain];
	    const Basis::Shell &blockL = auxbasis.shells().at( shellIndex );
	    size_t offsetAO = (blockL.start())*N*N;
	    if(spherical) offsetAO = (blockL.sphstart())*N*N;
	    
	    if(task != domain)continue;

#pragma omp barrier

#pragma omp master
	    if(debug){
		std::cout << pe.rank() << " " 
			  << pe.node().rank() 
			  << " working on task " << task 
			  << std::endl <<std::flush;
	    }

	    //////////////////////////////////////
	    //get a batch of (AO AO|aux) intergrals here
	    //////////////////////////////////////	    
#pragma omp single nowait
	    fineTimer.reset();
	    
	    AOBatch.getThreeCenterEriBatch(domain, blockRanges,
					   auxiliary_eri,sph_c,ptr_buff1,offsetAO);
	    
#pragma omp single //nowait
	    {profile.threeCenter +=fineTimer;
		blockIter = haveBlock.begin();
		haveBlock.insert(blockIter,domain);
		(++task);
	    }
	    
	    
	    
	    //////////////////////////////////////
	    //swap buffers and wait for write to finish 
	    // wait (if domain > 0) and a write is active
	    //////////////////////////////////////
#pragma omp single
	    {
		
		if(haveBlock.size() > 1 && activeWrite)async.wait();
		
		std::swap(ptr_buff1,ptr_buff2);
		std::swap(ptr_buff3,ptr_buff4);
		
		if(haveBlock.size() > 1){
		    
#if HAVE_CUBLAS
		    pGPU->synchronize();
		    //for (int idevice = 0; idevice < ndevice; idevice++){
		    //GPUerrchk( cudaSetDevice( idevice ) );
		    //GPUerrchk( cudaDeviceSynchronize() );
		    //}//idevice
#endif //HAVE_CUBLAS
		    
		    size_t domainTemp = haveBlock.back();
		    MOBatch.writeBlock(domainTemp,ptr_buff3);
		    haveBlock.pop_back();
		    
		    activeWrite = 1;
		    
		}//(haveBlock.size() > 1)
		
	    }//omp single
	    
 
	    
	    
	    
	    
	    
	    
	    //////////////////////////////////////
	    //AO-> MO transform here
	    //////////////////////////////////////
	    size_t offsetMO = (blockL.start())*no*nvs;
	    if(spherical)offsetMO = (blockL.sphstart())*no*nvs;
	    
#pragma omp master
	    fineTimer2.reset();
	    
	    if(ndevice>0 && omp_get_thread_num() != 0)goto skip;
	    
	    MOBatch.transformBlock(domain, offsetAO,offsetMO,
				   ptr_buff2, ptr_buff4 );
	    
	skip: //for gpu: the transformation is setup for one cpu thread only
	    
#pragma omp master
	    {profile.t2Trans +=fineTimer2;
		progress.jump(domain);}
	    //progress.jump(domain+1);}
	    
	}//domain




	    //////////////////////////////////////
	    //Write the last bit of transformed integrals
	    //////////////////////////////////////
	
	//i believe this barrier is needed in case the asyc.put is called AFTER the async.wait
#pragma omp barrier

	
#pragma omp single
	{
	    
	    if(haveBlock.size() > 0){

#if HAVE_CUBLAS
		pGPU->synchronize();
		//for (int idevice = 0; idevice < ndevice; idevice++){
		//GPUerrchk( cudaSetDevice( idevice ) );
		//GPUerrchk( cudaDeviceSynchronize() );
		//}//idevice
#endif //HAVE_CUBLAS
		
		

		size_t domainTemp = haveBlock[0];
		MOBatch.writeBlock(domainTemp,ptr_buff4);
		
		fineTimer2.reset();
		async.wait();
		profile.t2Write +=fineTimer2;
		
		
	    }//(haveBlock.size() > 0)
	    
	}//omp single
	delete [] threec_cart_to_sph;
	
	
    }//pe.node().rank() == 0

   //}
//bolt:;
    
} //scope: three-center Eri and transform


 pe.barrier();

 //make sure to indicate things are finished
 if(pe.rank() == 0)progress.jump(blockRanges.size()-1);

#if !HAVE_CUBLAS

	delete [] ptr_buff1;
	delete [] ptr_buff2;
	delete [] ptr_buff3;
	delete [] ptr_buff4;

#elif HAVE_CUBLAS

GPUerrchk( cudaFreeHost(ptr_buff1) );
GPUerrchk( cudaFreeHost(ptr_buff2) );
GPUerrchk( cudaFreeHost(ptr_buff3) );
GPUerrchk( cudaFreeHost(ptr_buff4) );

#endif

 profile.threeCTrans +=stepTimer;


cout<< "time for three_center eri     " << profile.threeCenter  << std::endl;
cout<< "time for 2/3 transform        " << profile.t2Trans  << std::endl;
cout<< "time for 2/3 write            " << profile.t2Write <<std::endl;
#if HAVE_CUBLAS
cout<< "time for gpu synchronization  " << profile.gpuSyncTrans << std::endl;
#endif
cout << std::endl<< "total 3C-2E and 1/3+2/3    " << profile.threeCTrans <<std::endl<< std::endl;





double *ptr_bj;
double *ptr_bj2;
double *ptr_bi;
double *ptr_bi2;

 //lbr
//#undef HAVE_CUBLAS
//#define HAVE_CUBLAS 0


#if !HAVE_CUBLAS

ptr_bj = new double[(nv+ns)*nl*stride];
ptr_bj2 = new double[(nv+ns)*nl*stride];
ptr_bi = new double[(nv+ns)*nl*stride];
ptr_bi2 = new double[(nv+ns)*nl*stride];

#elif HAVE_CUBLAS

GPUerrchk( cudaMallocHost( &ptr_bj, (size_t)((nv+ns)*nl*stride*sizeof(double)) ));
GPUerrchk( cudaMallocHost( &ptr_bj2, (size_t)((nv+ns)*nl*stride*sizeof(double)) ));
GPUerrchk( cudaMallocHost( &ptr_bi, (size_t)((nv+ns)*nl*stride*sizeof(double)) ));
GPUerrchk( cudaMallocHost( &ptr_bi2, (size_t)((nv+ns)*nl*stride*sizeof(double)) ));

#endif


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//start big comment out here

// #if HAVE_CUBLAS


// //vector of device pointers to (ia|Q) integrals
// std::vector< std::vector<double*> > dv1_IA(ndevice, std::vector<double*>(nstreams));

// //vector of device pointers to cholesky decomposed coulomb metric
// std::vector<double*> dv1_METRIC;

// //vector of device pointers to scratch data
// std::vector< std::vector<double*> > dv1_TEMP(ndevice, std::vector<double*>(nstreams));




// for (int idevice = 0; idevice < ndevice; idevice++) {

//     GPUerrchk( cudaSetDevice( deviceid[idevice] ) );

//     for(int istream = 0; istream < nstreams; istream++){

// 	//device storage for 3/3 transformed ERIs (each stream creates unique integrals)
// 	double *ptr_TEMP = 0;
// 	int matrixSizeTEMP = nvs*nl;
// 	GPUerrchk( cudaMalloc((void **)&ptr_TEMP, matrixSizeTEMP * sizeof(double)) );
// 	dv1_TEMP[idevice][istream] = ptr_TEMP;
	
// 	//device storage for (ia|Q) integrals and then (ia|Q)L-1 tranposed integrals
// 	double *ptr_IA = 0;
// 	int matrixSizeIA = nl*nvs;
// 	GPUerrchk( cudaMalloc((void **)&ptr_IA, matrixSizeIA * sizeof(double)) );
// 	dv1_IA[idevice][istream] = ptr_IA;

//     }//istream
   
//     //device storage for cholesky decomposed coulomb metric 
//     double *ptr_METRIC = 0;
//     int matrixSizeMETRIC = nl*nl;
//     GPUerrchk( cudaMalloc((void **)&ptr_METRIC, matrixSizeMETRIC * sizeof(double)) );
//     dv1_METRIC.push_back(ptr_METRIC);  
    
//  }//idevice

// #endif //HAVE_CUBLAS




// // 3/3 transformation
// // 1) perform 3/3 transformation to create b^ia_n matricies

//  stepTimer.reset();
//  {//scope: 3/3 transformation

//      //L^-1 matrix : Lm1 attaches to data_twoc pointer from way back when
//      //     MapMatrixXd Lm1(data_twoc,nl , nl);
//      //     Lm1.transposeInPlace(); 


// #if HAVE_CUBLAS
//      for(int idevice = 0; idevice < ndevice; idevice++){
// 	 GPUerrchk( cudaSetDevice( idevice ) );
// 	 //copy coulomb metric to device
// 	 GPUerrchk( cudaMemcpy(dv1_METRIC[idevice],data_twoc, nl*nl*sizeof(double),
// 			       cudaMemcpyHostToDevice) );    
//      }//idevice

//  #endif //HAVE_CUBLAS



//      size_t end= std::min(no, (size_t)(stride));
//      size_t start[] = { 0, 0 };
//      size_t finish[] = { end*(nv+ns), nl };
//      cchem::ri::async async;

//      async.get( 0, end*(nv+ns) , nl , *Vp, ptr_bj2);

// #if HAVE_CUBLAS
//      omp_set_num_threads(1);
// #endif //HAVE_CUBLAS

// #pragma omp parallel
//      {

// 	 for(int iocc = 0;iocc<no;iocc+=stride){
// 	     size_t end= std::min(no, (size_t)(iocc+stride));
// 	     size_t start[] = { iocc*(nv+ns), 0 };
// 	     size_t finish[] = { end*(nv+ns), nl };
	     
// 	     std::vector<size_t> startGet(2), finishGet(2), startPut(2), finishPut(2);

// #pragma omp single
// 	     {
// 		 fineTimer2.reset();
// 		 async.wait();
// 		 std::swap(ptr_bj, ptr_bj2);


// 		 size_t future_end   = std::min(no,(size_t)(iocc +stride +stride));
// 		 size_t future_start = std::min((size_t)(iocc +stride),future_end);

// 		 if(stride < no){
// 		     if(iocc == 0 ){

// 			 startGet[0] = future_start*(nv+ns);
// 			 startGet[1] = 0;
// 			 finishGet[0] = future_end*(nv+ns);
// 			 finishGet[1] = nl;

// 			 async.get( startGet, finishGet, *Vp, ptr_bj2);

// 			 profile.t2Read += fineTimer2;
// 		     }else{ //(iocc == 0)

// 			 std::swap(ptr_bi, ptr_bi2);
// 			 if(future_start < no){

// 			     startGet[0] = future_start*(nv+ns);
// 			     startGet[1] = 0;
// 			     finishGet[0] = future_end*(nv+ns);
// 			     finishGet[1] = nl;

// 			     startPut[0] = 0;
// 			     startPut[1] = (iocc-stride)*(nv+ns);
// 			     finishPut[0] = nl;
// 			     finishPut[1] = (iocc)*(nv+ns);

// 			     async.get_put( startGet, finishGet,
//                                             startPut, finishPut,
//                                             *Vp, ptr_bj2,
//                                             *VpT, ptr_bi2);

// 			     profile.t2Readt3Write += fineTimer2;
// 			 }else{ //(future_start < no)
			     
// 			     startPut[0] = 0;
// 			     startPut[1] = (iocc-stride)*(nv+ns);
// 			     finishPut[0] = nl;
// 			     finishPut[1] = (iocc)*(nv+ns);

// 			     async.put( startPut, finishPut, *VpT, ptr_bi2);


// 			     profile.t3Write += fineTimer2;			     
// 			 } //(future_start < no)
			 
// 		     } //(iocc == 0)
// 		 }// (stride < no)

// 	     }//omp single
	     
// 	     int remaining = std::min(stride,(int)(end-iocc) );

// 	     if(omp_get_thread_num() == 0)fineTimer.reset();
  

// #if !HAVE_CUBLAS	     

// #pragma omp for schedule(dynamic,1)
// 	     for (int s=0; s<end-iocc; s++){

// 		 //form B_ia^Q = (ia|Q)L-1 3/3 transformed integrals		 
// 		 cblas_dtrmm(CblasColMajor,
// 		 	     CblasRight,
// 			     CblasUpper,
// 		 	     CblasNoTrans,
// 		 	     CblasNonUnit,
// 		 	     (int)(nv+ns),(int)nl,
// 		 	     1.0,data_twoc,(int)nl,
// 		 	     &ptr_bj[s*(nv+ns)], (int)((nv+ns)*remaining) );

// 		 //this simply transposes B_ia^Q integrals to a more friendly form (for later reads/writes)
// 		 // (ia)x(nl) --> (nl)x(ia)
// 		 for (int ic = 0; ic < (nv+ns); ic++){
// 		     cblas_dcopy( (int)nl, &ptr_bj[s*(nv+ns) +ic], (int)((nv+ns)*remaining), &ptr_bi[s*(nv+ns)*nl + ic*nl] ,(int)1 );
// 		 }//ic
		 
// 	     }//s

//  #elif HAVE_CUBLAS

// #pragma omp single
// 	     {
// 		 for(int idevice = 0; idevice < ndevice; idevice++){
// 		     gpu_counter[idevice] = 0; //reset device counter
// 		 }//idevice
		 
// 	     } //omp single

//    for (int s=0; s<end-iocc; s++){

//        int hard_device = s%ndevice;

//        GPUerrchk( cudaSetDevice( hard_device ) );
//        int new_stream = gpu_counter[hard_device]%nstreams;		 
//        gpu_counter[hard_device] += 1;
//        cublasSetStream(dv1_handle[hard_device], dv1_streams[hard_device][new_stream]);
		 
//        double alpha = 1.0;
//        double beta = 0.0;

//        //copy forward i block of (ia|Q) 2/3 transformed integrals 
//        //    where Q=1,nl and a=1,nvs
//        GPUerrchk( cudaMemcpy2DAsync(dv1_IA[hard_device][new_stream],nvs*sizeof(double),
// 				    &ptr_bj[s*nvs],remaining*nvs*sizeof(double),
// 				    nvs*sizeof(double),nl,
// 				    cudaMemcpyHostToDevice,
// 				    dv1_streams[hard_device][new_stream]) );

//        //form B_ia^Q = (ia|Q)L-1 3/3 transformed integrals
//        cublasDtrmm(dv1_handle[hard_device],
//        		   CUBLAS_SIDE_RIGHT,
//        		   CUBLAS_FILL_MODE_UPPER,
//        		   CUBLAS_OP_N,
//        		   CUBLAS_DIAG_NON_UNIT,
//        		   nvs, nl,
//        		   &alpha,dv1_METRIC[hard_device],nl,
//        		   dv1_IA[hard_device][new_stream], nvs,
//        		   dv1_TEMP[hard_device][new_stream], nvs);

//        //this simply transposes B_ia^Q integrals to a more friendly form (for later reads/writes)
//        // (ia)x(nl) --> (nl)x(ia)
//        cublasStatus_t status = 
//        	   cublasDgeam(dv1_handle[hard_device],
//        		       CUBLAS_OP_T, CUBLAS_OP_T,
//        		       nl, nvs,
//        		       &alpha,
//        		       dv1_TEMP[hard_device][new_stream], nvs,
//        		       &beta, 
//        		       dv1_TEMP[hard_device][new_stream], nvs,
//        		       dv1_IA[hard_device][new_stream], nl );

//        //copy back B_ia^Q 
//        GPUerrchk( cudaMemcpyAsync(&ptr_bi[s*nl*(nv+ns)],
// 				  dv1_IA[hard_device][new_stream], 
// 				  nl*nvs*sizeof(double),
// 				  cudaMemcpyDeviceToHost,
// 				  dv1_streams[hard_device][new_stream]) );    
       
//    }//s


//    fineTimer2.reset();
//    for (int idevice = 0; idevice < ndevice; idevice++){

//        GPUerrchk( cudaSetDevice( idevice ) );
//        GPUerrchk( cudaDeviceSynchronize() );

//    }//idevice
//    profile.gpuSync33 += fineTimer2;

// #endif 

// 	     if(omp_get_thread_num() == 0)profile.t3Trans += fineTimer;
	     
// 	     if( (iocc+stride) >= no){
// #pragma omp single
// 		 { //last write

// 		     fineTimer2.reset();
// 		     async.wait();
// 		     profile.t2Readt3Write += fineTimer2;		     

// 		     fineTimer2.reset();

// 		     size_t start2[] = { 0, iocc*(nv+ns) };
// 		     size_t finish2[] = { nl, end*(nv+ns) };
// 		     VpT->put(ptr_bi, start2, finish2);

// 		     profile.t3Write += fineTimer2;

// 		 } //omp single

// 	     } //(iocc+stride) >= no)
	     
// 	 } //iocc
	 
//      } //omp parallel
     
//  } //scope: 3/3 transformation

//  profile.coeffTrans += stepTimer;

// //done with coulomb metric
// delete data_twoc;






// #if HAVE_CUBLAS
// for(int idevice = 0; idevice < ndevice; idevice++){

//     GPUerrchk( cudaSetDevice( idevice ) );

//     GPUerrchk( cudaFree(dv1_METRIC[idevice]) );

//     for(int istream = 0; istream < nstreams; istream++){

// 	GPUerrchk( cudaFree( dv1_TEMP[idevice][istream]) );
// 	GPUerrchk( cudaFree(dv1_IA[idevice][istream]) );

	
//     }//istream
    
//  }//idevice
// #endif //HAVE_CUBLAS

// std::cout<< std::endl
//          << "time for 2/3 read              " << profile.t2Read  << std::endl;
// std::cout<< "time for 3/3 transform         " << profile.t3Trans << std::endl;
// std::cout<< "time for 2/3 read + 3/3 write  " << profile.t2Readt3Write << std::endl;
// std::cout<< "time for 3/3 write             " << profile.t3Write <<std::endl;
// #if HAVE_CUBLAS
// std::cout<< "time gpu snchronization        " << profile.gpuSync33 <<std::endl<< std::endl;
// #endif
// std::cout<< "total 2/3 to  3/3 time         " << profile.coeffTrans <<std::endl<< std::endl;


//end big comment out here
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// pe.barrier() ;
//     std::cout << pe.rank() << " made it! " << std::endl << std::flush;
// pe.barrier();

	// 1) form four-center ERIs
	// 2) compute ri-mp2 or ri-zapt energies

	{//scope: form four-center ERIs and compute energy

		BOOST_AUTO(const &ea, wf.e(wf.double_occ()));
		BOOST_AUTO(const &ev, wf.e(wf.single_virtual()));
 
		double *ptr_eri = new double[(nv+ns)*(nv+ns)];

		//e_m: (ts|ts) exchange integrals (t,s: single occcupied)
		//  used to modify energy denominators open-shell terms
		//pfock: used for fock like term
		std::vector<double> e_m;
		double *pfock;
//				std::vector<double *> pfock;
		detail::ri_openshell_work *os_work_ptr;
		if(ns>0){
			for (int i = 0; i < ns; ++i) e_m.push_back(0.0);
			os_work_ptr = new detail::ri_openshell_work(nd,ns,nv,nl,memory);
			os_work_ptr->ri_exchange_eri(Vp,ptr_bi,ptr_eri,e_m, data_twoc);
//			os_work_ptr->ri_exchange_eri2(VpT,ptr_bi,ptr_eri,e_m);
			pfock = os_work_ptr->get_pfock();
		}//(ns>0)

		delete [] ptr_eri;








//utility::Progress progress;
 if (pe.rank() == 0) progress.reset(no);


 cchem::ri::async async;


 stepTimer.reset();

#if HAVE_CUBLAS
 omp_set_num_threads(1);
#endif //HAVE_CUBLAS


   pe.task().reset();
   detail::Thread::Task<Parallel::Task&> task(pe.task());
   task.reset();

   //intdications how many Jblocks are stored in buffer (if 0, nothing to do)

   std::vector<int>haveTask;
   std::vector<int>::iterator taskIter;


   omp_lock_t *pLocks = new omp_lock_t[ omp_get_max_threads() ];

   // for(int i = 0; i < omp_get_max_threads(); i++){
   //     omp_init_lock(&pLocks[i]);
   //     //       omp_set_lock(&pLocks[i]);
   // }
 
   //   omp_unset_lock(&pLocks[0]);

   // omp_set_num_threads(1);
#pragma omp parallel
 //if (pe.node().rank() == 0 && pe.rank() == 0) { //node (local) rank 
 //   if (pe.node().rank() == 0 && pe.rank() == 1) { //node (local) rank 0
 if (pe.node().rank() == 0) { //node (local) rank 0




#if !HAVE_CUBLAS

   detail::Energy_Term energy_term(nd,ns,nv,nl, stride, strideRow);

   detail::RIMP2Engine CPUEngine(nd, ns, nv, nl,
   		   	   stride, strideRow, data_twoc,
   			   &(ea.data()[0]),&(ev.data()[nd]),
				 e_m, energy_term,pLocks);

#endif

#if HAVE_CUBLAS

   if(debug){
       std::cout << " ndevice       " <<  ndevice << " " << std::endl
		 << " nstreams      " << nstreams << " "<< std::endl
		 << " cudaBlocks    " << cudaBlocks << " " << std::endl
		 << " cudaThreadsPB " << cudaThreadsPB << std::endl << std::flush;
   }

   detail::RIMP2Engine GPUEngine(nv,ns,nl,nd,stride,no,
   				 &(ea.data()[0]),&(ev.data()[0]),
   				 e_m, data_pqm1,
				 pGPU);
#endif



   // //intdications how many Jblocks are stored in buffer (if 0, nothing to do)

   // std::vector<int>haveTask;
   // std::vector<int>::iterator taskIter;

// #pragma omp master 
//    ++task;





   //for IO
   std::vector<size_t> startGet(2), finishGet(2), startPut(2), finishPut(2);

#pragma omp single
   {


     // fineTimer.reset();
     // size_t jblock_start = 0;
     // size_t jblock_end =   std::min( (size_t)(stride2),no);
     

     // startGet[0] = jblock_start*nvs;
     // startGet[1] = 0;
     // finishGet[0] = jblock_end*nvs;
     // finishGet[1] = nl;

     // async.get( startGet, finishGet, *Vp, ptr_bj2);

     // profile.eriReads += fineTimer;




       ++task;
       size_t jblock_start = task*stride;
       size_t jblock_end = std::min( (size_t)(task+1)*stride, no );
       if(jblock_start < no){
     

     startGet[0] = jblock_start*nvs;
     startGet[1] = 0;
     finishGet[0] = jblock_end*nvs;
     finishGet[1] = nl;

     async.get( startGet, finishGet, *Vp, ptr_bj2);

     profile.eriReads += fineTimer;

     haveTask.push_back(task);
       }//(jblock_start < no)



   }//omp single



   // for (int jblock = 0; jblock < no; jblock+=stride2){
   //   size_t jblock_end = std::min( (size_t)(jblock+stride2) ,no);


   while( haveTask.size() ){
   int jblock = haveTask.back()*stride;
   size_t jblock_end = std::min( (size_t)( jblock+stride) ,no);

#pragma omp single
     {

       fineTimer2.reset();
       // size_t j_future_start = std::min( (size_t)(jblock+stride2) ,no);
       // size_t j_future_end = std::min( (size_t)(jblock+stride2+stride2) ,no);


       // async.wait();
       // std::swap(ptr_bj,ptr_bj2);


       // startGet[0] = j_future_start*nvs;
       // startGet[1] = 0;
       // finishGet[0] = j_future_end*nvs;
       // finishGet[1] = nl;

       // if(jblock_end != no){

       // 	   async.get( startGet, finishGet, *Vp, ptr_bj2);

       // }//(jblock_end != no)





       async.wait();
       std::swap(ptr_bj,ptr_bj2);

       ++task;

       //       std::cout << pe.rank() << " " << task << std::endl << std::flush;

        size_t j_future_start = std::min( (size_t)(task*stride), no );
        size_t j_future_end  = std::min( (size_t)((task+1)*stride), no );
        if( j_future_start < no){



       startGet[0] = j_future_start*nvs;
       startGet[1] = 0;
       finishGet[0] = j_future_end*nvs;
       finishGet[1] = nl;

       async.get( startGet, finishGet, *Vp, ptr_bj2);

        taskIter = haveTask.begin();
        haveTask.insert(taskIter,task);
        }//( j_future_start < no)

       profile.eriReads += fineTimer2;
     }//omp single

     int j_remain = std::min(stride,(int)(jblock_end-jblock) );

#if !HAVE_CUBLAS

     if(omp_get_thread_num() == 0)fineTimer2.reset();

#pragma omp single
     CPUEngine.resetLocks();

     CPUEngine(ptr_bj,jblock,jblock_end,j_remain);

     if(omp_get_thread_num() == 0)profile.eriFormationii += fineTimer2;

#elif HAVE_CUBLAS

    GPUEngine(jblock_end, jblock,ptr_bj,nvs,j_remain);

#endif //HAVE_CUBLAS



#pragma omp single
     if(jblock+stride < no){
       fineTimer.reset();
       size_t iblock_start = std::min( (size_t)(jblock+stride),no);
       size_t iblock_end =   std::min( (size_t)(jblock+stride+stride),no);

       async.wait();

       startGet[0] = iblock_start*nvs;
       startGet[1] = 0;
       finishGet[0] = iblock_end*nvs;
       finishGet[1] = nl;
       async.get( startGet, finishGet, *Vp, ptr_bi2);

       profile.eriReads += fineTimer;
     }//(jblock+stride < no)


     for (int iblock = jblock+stride; iblock < no; iblock+=stride){
       size_t iblock_end = std::min( (size_t)(iblock+stride) ,no);
       int i_remain = std::min(stride,(int)(iblock_end-iblock) );

#pragma omp single
       {
	 fineTimer2.reset();
	 size_t i_future_start = std::min( (size_t)(iblock+stride) ,no);
	 size_t i_future_end = std::min( (size_t)(iblock+stride+stride) ,no);

	 async.wait();
	 std::swap(ptr_bi,ptr_bi2);


	 startGet[0] = i_future_start*nvs;
	 startGet[1] = 0;
	 finishGet[0] = i_future_end*nvs;
	 finishGet[1] = nl;

	 if(iblock_end != no){
	   async.get( startGet, finishGet, *Vp, ptr_bi2);
	 }


	 profile.eriReads += fineTimer2;
       }//omp single

       int num_chunks = i_remain/strideRow + ( (i_remain %strideRow) != 0);

#if !HAVE_CUBLAS

       if(omp_get_thread_num() == 0)fineTimer2.reset();

       CPUEngine(ptr_bj,ptr_bi,num_chunks,i_remain,j_remain,iblock,jblock,jblock_end);

       if(omp_get_thread_num() == 0)profile.eriFormationij += fineTimer2;

#elif HAVE_CUBLAS

       GPUEngine(num_chunks,i_remain,j_remain,strideRow,ptr_bi,ptr_bj,nvs,iblock,jblock);


#endif //HAVE_CUBLAS

     }//iblock

#pragma omp master
     progress.jump(jblock+j_remain-1);
     //progress.jump(jblock+j_remain);


#if HAVE_CUBLAS
	  fineTimer.reset();
	  for (int idevice = 0; idevice < ndevice; idevice++){
	      GPUerrchk( cudaSetDevice( idevice ) );
	      GPUerrchk( cudaDeviceSynchronize() );
	  }//idevice
	  profile.gpuSyncEri += fineTimer;
#endif //HAVE_CUBLAS


	  //remove block that we are done with
#pragma omp single
	  haveTask.pop_back();




   }//jblock


   //#pragma omp barrier
#if !HAVE_CUBLAS
       CPUEngine.reduceThreadEnergies( &(e_term[0]), pfock, pe );
#endif


     if(ns > 0){
#pragma omp single
       delete os_work_ptr; //must delete since os_work_ptr is not a smart pointer
     }//(ns > 0)

#if HAVE_CUBLAS
     //gathers energy contributions and deallocates memory on gpu
     GPUEngine.Deallocate(&e_term[0], pe);

#endif //HAVE_CUBLAS


 }//(pe.node().rank() == 0)

 if(pe.rank() == 0)progress.jump(no);

  profile.eri += stepTimer;

}//scope: form four-center ERIs and compute energy



cout<< "time for:"<<std::endl;
cout<< "   total eri                  "<<profile.eri<<std::endl;
cout<< "   eri_reads                  "<<profile.eriReads<<std::endl;
cout<< "   eri i=j pair contribution  "<<profile.eriFormationii<< std::endl;
cout<< "   eri i>j pair contribution  "<<profile.eriFormationij<< std::endl;
#if !HAVE_CUBLAS
cout<< "   CPU ERI dgemms             "<<profile.eriFormation<< std::endl;
cout<< "   CPU energy accumulation    "<<profile.energyAccum<< std::endl;
#elif HAVE_CUBLAS
cout<< "   gpu synchonization         "<<profile.gpuSyncEri<< std::endl;
#endif


//	cout <<std::setprecision(10);
//std::cout << std::endl << pe.rank() <<  ": corr energy: " << e_term[0] << std::endl;
pe.reduce("+",&e_term[0],(size_t)7);

#if !HAVE_CUBLAS

	delete [] ptr_bi;
	delete [] ptr_bi2;
	delete [] ptr_bj;
	delete [] ptr_bj2;

#elif HAVE_CUBLAS

GPUerrchk( cudaFreeHost(ptr_bi) );
GPUerrchk( cudaFreeHost(ptr_bi2) );
GPUerrchk( cudaFreeHost(ptr_bj) );
GPUerrchk( cudaFreeHost(ptr_bj2) );

#endif

E = std::accumulate(e_term,e_term+7,0.0);

used = memory.used();
memory.clear();

    {
	cout << "    memory: " << used/(1<<20) << " MB" << std::endl;
	cout << std::endl;
	BOOST_PROFILE_DUMP(cout);
    }


    //print out energy decomposition for ZAPT (i.e. ns > 0)
    if(ns > 0){
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

profile.global += globalTimer;
cout << std::endl << "Total time spent inside RI-MP2: " << profile.global << std::endl;

        // std::cout << std::endl << "erasing array ri-mp2.vp(b,qs,l)" << std::endl;
        // rt.arrays().erase< Array<double> >("rimp2.vp(b,qs,l)");
        // rt.arrays().erase< Array<double> >("rimp2.vp2(b,qs,l)");


//herr_t hstatus;
//= H5Fclose ("rimp2.vp(b,qs,l)");
//H5close();

    // std::cout << "done erasing arrays" << std::endl;

    return E;

}//rimp2::energy(Wavefunction wf, Runtime &rt)
} // namespace cchem

