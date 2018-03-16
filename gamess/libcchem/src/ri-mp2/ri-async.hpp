/*
 * ri-async.hpp
 *
 *  Created on: Jun 16, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_MP2_RI_ASYNC_HPP_
#define LIBCCHEM_SRC_MP2_RI_ASYNC_HPP_


#include "god.hpp"


//#include <stdlib.h>
//#include "runtime.hpp"

#include "thread.hpp"
#include "core/wavefunction.hpp"
#include "parallel.hpp"

#if HAVE_CUBLAS
#include <device.hpp>
#endif


namespace cchem{
namespace ri{






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
	//    	void get(size_t k, size_t m, size_t n, size_t b,
	//    		 const Array<double> &V, double *data) {
	void get(size_t i, size_t j, size_t l,
			const Array<double> &V, double *data) {
		thread_.reset(new boost::thread(launch, i, j, l, boost::cref(V), data));
	}

	void get(std::vector<size_t> start_get, std::vector<size_t> finish_get,
			const Array<double> &V, double *data) {
		thread_.reset(new boost::thread(launch_get_2, start_get, finish_get, boost::cref(V), data));
	}

	void get_put(size_t i, size_t j, size_t l,size_t iput, size_t jput,
			Array<double> &V, double *data, double *data_put) {
		thread_.reset(new boost::thread(launch_get_put, i, j, l, iput, jput, boost::ref(V), data, data_put));
	}

	void get_put(std::vector<size_t> start_get, std::vector<size_t> finish_get,
			std::vector<size_t> start_put, std::vector<size_t> finish_put,
			Array<double> &V, double *data_get,
			Array<double> &W, double *data_put) {
		thread_.reset(new boost::thread(launch_get_put_2,
				start_get, finish_get,
				start_put, finish_put,
				boost::ref(V), data_get,
				boost::ref(W), data_put));
	}

	void put(size_t l,size_t iput, size_t jput,
			Array<double> &V, double *data_put) {
		thread_.reset(new boost::thread(launch_put, l, iput, jput, boost::ref(V), data_put));
	}

	void put(std::vector<size_t> start_put, std::vector<size_t> finish_put,
			Array<double> &V, double *data_put) {
		thread_.reset(new boost::thread(launch_put_2, start_put, finish_put, boost::ref(V), data_put));
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
	//		void launch(size_t k, size_t m, size_t n, size_t b,
	//				const Array<double> &V, double *data) {
	void launch(size_t i, size_t j, size_t l,
			const Array<double> &V, double *data) {
		//    		BOOST_PROFILE_LINE;
	    //		utility::timer timer;
		size_t start[] = { i, 0 };
		size_t finish[] = { j, l };
		V.get(data, start, finish);
		//std::cout << "I/O: " << (N*N*B*8)/(timer*1e6) << std::endl;
	}


	static
	void launch_get_put(size_t i, size_t j, size_t l, size_t iput, size_t jput,
			Array<double> &V, double *data, double *data_put) {
		//    		BOOST_PROFILE_LINE;
	    //		utility::timer timer;
		size_t start[] = { i, 0 };
		size_t finish[] = { j, l };
		V.get(data, start, finish);

		size_t start2[] = { iput, 0 };
		size_t finish2[] = { jput, l };
		V.put(data_put, start2, finish2);

	} //launch_get_put

	static
	void launch_put(size_t l, size_t iput, size_t jput,
			Array<double> &V,double *data_put) {

		size_t start[] = { iput, 0 };
		size_t finish[] = { jput, l };
		V.put(data_put, start, finish);

	} //launch_get_put





	static
	void launch_get_2( std::vector<size_t> start_get, std::vector<size_t> finish_get,
			const Array<double> &V, double *data) {
		V.get(data, start_get, finish_get);
	}//launch_get_2

	static
	void launch_get_put_2(std::vector<size_t> start_get, std::vector<size_t> finish_get,
			std::vector<size_t> start_put, std::vector<size_t> finish_put,
			Array<double> &V, double *data,
			Array<double> &W, double *data_put) {
		V.get(data, start_get, finish_get);
		W.put(data_put, start_put, finish_put);
	} //launch_get_put_2

	static
	void launch_put_2(std::vector<size_t> start_put, std::vector<size_t> finish_put,
			Array<double> &V,double *data_put) {
		V.put(data_put, start_put, finish_put);
	} //launch_put_2


}; //async

}//namespace ri
}//namespace cchem
















//#include <ri-lagrangian.hpp> //for friend class




namespace cchem{
namespace rimp2_gradient{
namespace detail {

    //using cchem::Thread;


/*
  @brief Class that encapsulates asychronous IO with chosen functor
  @author LBR
  @detail Class that encapsulates asychronous IO with chosen functor.
  The contructors are overloaded depending on the type of IO requested.
  @param NMB number of MB alloted for double buffered IO
  @param ArrayOne array pointer to data
  @param NA1RowB number of blocks the row are divided into
  @param NA1RowBstart which row to start accessing data pointed to by ArrayOne
  @param NA1ColB bumber of blocks the columns are divided into
  @param NA1ColBStart which column to start accessing data pointed to by ArrayOne
  @param ReadRow flag to specify if rows(1) or columns(0) are to be read
  @param debug debug flag (0,1)
  @param pe parallel environment
  @param printArray printArray access tables/information (default=0)
  @param MagicSwap swap arrays bounds when setting up access tables
  @param read_restrict restrict access tables from multiple arrays to be identical
  @param ArrayTwo array pointer to data
  @param NA2RowB number of blocks the row are divided into
  @param NA2RowBstart which row to start accessing data pointed to by ArrayTwo
  @param NA2ColB bumber of blocks the columns are divided into
  @param NA2ColBStart which column to start accessing data pointed to by ArrayTwo
  @param loopRestrict to force array reads to be read2>=read1
*/
    class DFAsync {

	typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
	
	typedef std::vector < std::pair <std::vector <size_t>, std::vector <size_t> > > ReadPos;
	typedef std::vector < std::pair < size_t, size_t> > Blocks;
	
	struct ArrayData{

	    //the array we're messing with
	    Array<double> *Array_;
	    //number of row blocks , number of rows per block
	    size_t NAB_, NABRows_;
	    //number of column blocks, number of columns per block
	    size_t NAColB_, NABCols_;
	    //read start and stop position  [e.g., { i, j } --> { k, l }]
	    ReadPos	AReadPos_;
	    //vector of start stop block positions (chunk of blocks)  [ 0 : 8 ]
	    Blocks ABlocks_;
	    
	    //largest chunk of blocks to be accessed (row or column), number of reads for for array
	    size_t MaxNAB_, NAReads_;
	    
	    //currently not used, but for unique reads (i.e. read(i) & read(j>=i) )
	    bool read_restrict_;
	    
	    bool loopRestrict_;
	    std::string arrayType_;
	    
	    ArrayData(Array<double> *Array, size_t NAB, size_t NAColB) :
		Array_(Array),
		NAB_(NAB),NABRows_(0),
		NAColB_(NAColB),NABCols_(0),
		AReadPos_(),ABlocks_(),
		MaxNAB_(0),NAReads_(0),
		read_restrict_(0),
		loopRestrict_(0),
		arrayType_(Array_->str())
	    {};
	    
	};//ArrayData
	
	
	size_t NMB_;
	Array<double> *ArrayOne_;
	Array<double> *ArrayTwo_;
	Array<double> *ArrayThr_;
	Array<double> *ArrayFou_;
	std::vector<ArrayData> Array_;
	
	//parallel environment
	Parallel &pe_;
	
    public:
	
	~DFAsync(){}
	
	DFAsync(size_t NMB,
		Array<double> *ArrayOne,
		size_t NA1RowB, size_t NA1RowBStart, size_t NA1ColB, size_t NA1ColBStart,
		bool ReadRow, bool debug, Parallel &pe, bool printArray = 0) :
	    NMB_(NMB),
	    ArrayOne_(ArrayOne),
	    ArrayTwo_(NULL),
	    ArrayThr_(NULL),
	    ArrayFou_(NULL),
	    pe_(pe)
	{
	    
	    Array_.push_back(ArrayData(ArrayOne_,NA1RowB-NA1RowBStart, NA1ColB-NA1ColBStart));
	    
	    size_t NWorkingBytes = NMB_*1000*1000;
	    bool async = 1;
	    
	    this->set_up_array_info(ReadRow, async, NWorkingBytes,
				    Array_[0], NA1RowBStart, NA1ColBStart);
	    
	    if(printArray){
		this->PrintArray(0);
	    }//printArray
	    
	}//RIAsync

	DFAsync(size_t NMB,
		Array<double> *ArrayOne,
		size_t NA1RowB, size_t NA1RowBStart, size_t NA1ColB, size_t NA1ColBStart,
		Array<double> *ArrayTwo,
		size_t NA2RowB, size_t NA2RowBStart, size_t NA2ColB, size_t NA2ColBStart,
		bool read_restrict, bool loopRestrict, bool ReadRow, bool MagicSwap,
		Parallel &pe, bool arrayDebug = 0) :
	    NMB_(NMB),
	    ArrayOne_(ArrayOne),
	    ArrayTwo_(ArrayTwo),
	    ArrayThr_(NULL),
	    ArrayFou_(NULL),
	    pe_(pe)
	{
	    
	    if(MagicSwap){
		std::swap(ArrayOne_,ArrayTwo_);
		std::swap(NA1RowB,NA2RowB);
		std::swap(NA1ColB,NA2ColB);
		std::swap(NA1RowBStart,NA2RowBStart);
		std::swap(NA1ColBStart,NA2ColBStart);
	    }
	    
	    Array_.push_back(ArrayData(ArrayOne_,NA1RowB-NA1RowBStart, NA1ColB-NA1ColBStart));
	    Array_.push_back(ArrayData(ArrayTwo_,NA2RowB-NA2RowBStart, NA2ColB-NA2ColBStart));
	    
	    size_t NWorkingBytes = NMB_*1000*1000;
	    bool async = 1;
	    
	    this->set_up_pair_array_info(read_restrict, ReadRow, async, NWorkingBytes,
					 Array_[0], NA1RowBStart, NA1ColBStart,
					 Array_[1], NA2RowBStart, NA2ColBStart);
	    
	    Array_[0].loopRestrict_ = loopRestrict;
	    
	    if(MagicSwap){
		std::iter_swap(Array_.begin() + 0,Array_.begin() + 1);
		std::swap(ArrayOne_,ArrayTwo_);
		std::swap(NA1RowB,NA2RowB);
		std::swap(NA1ColB,NA2ColB);
		std::swap(NA1RowBStart,NA2RowBStart);
		std::swap(NA1ColBStart,NA2ColBStart);
	    }
	    
	    if(arrayDebug){
		this->PrintArray(0);
		this->PrintArray(1);
	    }//arrayDebug

	}//DFAsync
	
	size_t getMaxNAB(int i){return Array_[i].MaxNAB_;}
	
	size_t getNABRows(int i){return Array_[i].NABRows_;}
	
	size_t getNABCols(int i){return Array_[i].NABCols_;}

	void set_up_array_info(bool &ReadRow, bool &async, size_t &NWorkingBytes,
			       ArrayData &A1,  size_t NA1RowBStart, size_t NA1ColBStart);
	
	void set_up_pair_array_info(bool &restrict, bool &ReadRow , bool &async, size_t & NWorkingBytes,
				    ArrayData &A1,  size_t NA1RowBStart,size_t NA1ColBStart,
				    ArrayData &A2,  size_t NA2RowBStart,size_t NA2ColBStart);
	
	void set_up_block_info(ArrayData &A,  size_t NARowBStart, size_t NAColBStart, bool ReadRow);
	
	void PrintArray(int i);
	
	template< typename OP >
	void DoOpAsync_R_W(OP Functor);

	template< typename OP >
	void DoOpAsync_R_W_Test(OP &Functor);

	template< typename OP >
	void DoOpAsync_R_W_NODE_PARALLEL(OP &Functor);

	template <typename FUNCTOR >
	void DoOpAsync_RRSerial(FUNCTOR &Functor);

	template <typename FUNCTOR >
	void DoOpAsync_RR(FUNCTOR &Functor);

	template< typename OP >
	void DoOpAsync_R(OP &Functor);

	template< typename OP >
	void DoOpAsync_R_NODE_PARALLEL(OP &Functor);

};//class DFAsync

}//detail
}//rimp2_gradient
}//cchem

template< typename OP >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_R_W(OP Functor) {

    ArrayData &A1 = Array_[0];
    ArrayData &A2 = Array_[1];

    //#if !HAVE_CUBLAS
    //storage for '2/3' transformed integrals
    double *ptr_ints = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];


    //storage for '2/3' transformed integrals
    double *ptr_ints2 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
	
    //storage for '3/3' transformed integrals
    double *ptr_Coeff = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
	
    //storage for '3/3' transformed integrals
    double *ptr_Coeff2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

    // #elif HAVE_CUBLAS

    // 		double *ptr_ints;
    // 		GPUerrchk(cudaMallocHost( &ptr_ints, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_ints2;
    // 		GPUerrchk(cudaMallocHost( &ptr_ints2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_Coeff;
    // 		GPUerrchk(cudaMallocHost( &ptr_Coeff, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_Coeff2;
    // 		GPUerrchk(cudaMallocHost( &ptr_Coeff2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // #endif

    utility::timer timer;

    utility::timer::value_type waitTime, mixedTime, readTime, writeTime, firstReadTime,
	finalWriteTime;;

    ::cchem::ri::async async;
    async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

    timer.reset();
    async.wait(); firstReadTime += timer;

#pragma omp parallel
    {
	    
	for(size_t read = 0; read < A1.NAReads_; read++ ){

#pragma omp single
	    {
		timer.reset();
		async.wait(); mixedTime += timer;

		std::swap(ptr_ints,ptr_ints2);
		    
		if(A1.NAReads_ > 1){
			
		    if(read == 0){
			async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
				  *(A1.Array_), ptr_ints2);
		    }else{
			std::swap(ptr_Coeff,ptr_Coeff2);
			if(read < A1.NAReads_-1){
				
			    async.get_put(A1.AReadPos_[read+1].first, 
					  A1.AReadPos_[read+1].second,
					  A2.AReadPos_[read-1].first, 
					  A2.AReadPos_[read-1].second,
					  *(A1.Array_), ptr_ints2,
					  *(A2.Array_), ptr_Coeff2);
				
				
			    // async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
			    // 	  *(A1.Array_), ptr_ints2);

			    // timer.reset();
			    // async.wait(); readTime += timer;

			    // async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
			    // 	  *(A2.Array_), ptr_Coeff2);

			    // timer.reset();
			    // async.wait(); writeTime += timer;
				
				
			}else{

			    async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
				      *(A2.Array_), ptr_Coeff2);

			    // timer.reset();
			    // async.wait(); writeTime += timer;

			}//(read < A1.NAReads_-1)
		    }//(read == 0)
		}//(A1.NAReads_ > 0)
	    }//omp single

	    Functor(ptr_ints,A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABRows_,
		    A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABCols_,
		    ptr_Coeff,A2.ABlocks_[read].first,A2.ABlocks_[read].second,A2.NABRows_,
		    A2.ABlocks_[read].first,A2.ABlocks_[read].second,A2.NABCols_);
		
	}//read

    }//omp parallel

    //NOT asynchronous
    timer.reset();
    A2.Array_->put(ptr_Coeff, A2.AReadPos_[A2.NAReads_-1].first,
		   A2.AReadPos_[A2.NAReads_-1].second);

    async.wait(); 

    finalWriteTime += timer;

    //	async.wait(); waitTime += timer;
    //async.wait(); writeTime += timer;
	
    //#if !HAVE_CUBLAS

    delete [] ptr_Coeff2;
    delete [] ptr_Coeff;
    delete [] ptr_ints2;
    delete [] ptr_ints;

    // #elif HAVE_CUBLAS

    // 		GPUerrchk( cudaFreeHost(ptr_Coeff2) );
    // 		GPUerrchk( cudaFreeHost(ptr_Coeff) );
    // 		GPUerrchk( cudaFreeHost(ptr_ints2) );
    // 		GPUerrchk( cudaFreeHost(ptr_ints) );

    // #endif

    waitTime = readTime + writeTime + mixedTime + firstReadTime;

    // if(pe_.rank() == 0){
    // 	std::cout << std::endl << "I/O total wait time : " << waitTime << std::endl;
    // 	std::cout << "read I/O time : " << readTime << std::endl;
    // 	std::cout << "write I/O time : " << writeTime << std::endl;
    // 	std::cout << "mixed I/O time : " << mixedTime << std::endl;
    // 	std::cout << "first read I/O time : " << firstReadTime << std::endl;
    // 	std::cout << "final write I/O time : " << finalWriteTime << std::endl;
    // }//(pe_.rank() == 0)

}//DoOpAsync_R_W



template< typename OP >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_R_W_Test(OP &Functor) {

    ArrayData &A1 = Array_[0];
    ArrayData &A2 = Array_[1];

    //#if !HAVE_CUBLAS
    //storage for '2/3' transformed integrals
    double *ptr_ints = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];


    //storage for '2/3' transformed integrals
    double *ptr_ints2 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
	
    //storage for '3/3' transformed integrals
    double *ptr_Coeff = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
	
    //storage for '3/3' transformed integrals
    double *ptr_Coeff2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

    // #elif HAVE_CUBLAS

    // 		double *ptr_ints;
    // 		GPUerrchk(cudaMallocHost( &ptr_ints, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_ints2;
    // 		GPUerrchk(cudaMallocHost( &ptr_ints2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_Coeff;
    // 		GPUerrchk(cudaMallocHost( &ptr_Coeff, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_Coeff2;
    // 		GPUerrchk(cudaMallocHost( &ptr_Coeff2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // #endif

    utility::timer timer;

    utility::timer::value_type waitTime, mixedTime, readTime, writeTime, firstReadTime,
	finalWriteTime;;

    ::cchem::ri::async async;
    async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

    timer.reset();
    async.wait(); firstReadTime += timer;

#pragma omp parallel
    {
	    
	for(size_t read = 0; read < A1.NAReads_; read++ ){

#pragma omp single
	    {
		timer.reset();
		async.wait(); mixedTime += timer;

		std::swap(ptr_ints,ptr_ints2);
		    
		if(A1.NAReads_ > 1){
			
		    if(read == 0){
			async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
				  *(A1.Array_), ptr_ints2);
		    }else{
			std::swap(ptr_Coeff,ptr_Coeff2);
			if(read < A1.NAReads_-1){
				
			    async.get_put(A1.AReadPos_[read+1].first, 
					  A1.AReadPos_[read+1].second,
					  A2.AReadPos_[read-1].first, 
					  A2.AReadPos_[read-1].second,
					  *(A1.Array_), ptr_ints2,
					  *(A2.Array_), ptr_Coeff2);
				
				
			    // async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
			    // 	  *(A1.Array_), ptr_ints2);

			    // timer.reset();
			    // async.wait(); readTime += timer;

			    // async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
			    // 	  *(A2.Array_), ptr_Coeff2);

			    // timer.reset();
			    // async.wait(); writeTime += timer;
				
				
			}else{

			    async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
				      *(A2.Array_), ptr_Coeff2);

			    // timer.reset();
			    // async.wait(); writeTime += timer;

			}//(read < A1.NAReads_-1)
		    }//(read == 0)
		}//(A1.NAReads_ > 0)
	    }//omp single

	    Functor(ptr_ints,A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABRows_,
		    A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABCols_,
		    ptr_Coeff,A2.ABlocks_[read].first,A2.ABlocks_[read].second,A2.NABRows_,
		    A2.ABlocks_[read].first,A2.ABlocks_[read].second,A2.NABCols_);
		
	}//read

    }//omp parallel

    //NOT asynchronous
    timer.reset();

    // std::cout << (A2.AReadPos_[A2.NAReads_-1].first[0]) 
    // 	  << " " << (A2.AReadPos_[A2.NAReads_-1].first[1]) 
    // 	  << " " << (A2.AReadPos_[A2.NAReads_-1].second[0])
    // 	  << " " << (A2.AReadPos_[A2.NAReads_-1].second[1])
    // 	  << std::flush << std::endl;
    A2.Array_->put(ptr_Coeff, A2.AReadPos_[A2.NAReads_-1].first,
		   A2.AReadPos_[A2.NAReads_-1].second);



    //async.wait(); 

    finalWriteTime += timer;

	
    //#if !HAVE_CUBLAS

    delete [] ptr_Coeff2;
    delete [] ptr_Coeff;
    delete [] ptr_ints2;
    delete [] ptr_ints;

    // #elif HAVE_CUBLAS

    // 		GPUerrchk( cudaFreeHost(ptr_Coeff2) );
    // 		GPUerrchk( cudaFreeHost(ptr_Coeff) );
    // 		GPUerrchk( cudaFreeHost(ptr_ints2) );
    // 		GPUerrchk( cudaFreeHost(ptr_ints) );

    // #endif

    waitTime = readTime + writeTime + mixedTime + firstReadTime;

    // if(pe_.rank() == 0){
    // 	std::cout << std::endl << "I/O total wait time : " << waitTime << std::endl;
    // 	std::cout << "read I/O time : " << readTime << std::endl;
    // 	std::cout << "write I/O time : " << writeTime << std::endl;
    // 	std::cout << "mixed I/O time : " << mixedTime << std::endl;
    // 	std::cout << "first read I/O time : " << firstReadTime << std::endl;
    // 	std::cout << "final write I/O time : " << finalWriteTime << std::endl;
    // }//(pe_.rank() == 0)

}//DoOpAsync_R_W_Test





template< typename OP >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_R_W_NODE_PARALLEL(OP &Functor) {

    ArrayData &A1 = Array_[0];
    ArrayData &A2 = Array_[1];

    //#if !HAVE_CUBLAS
    //storage for '2/3' transformed integrals
    double *ptr_ints = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];


    //storage for '2/3' transformed integrals
    double *ptr_ints2 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
	
    //storage for '3/3' transformed integrals
    double *ptr_Coeff = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
	
    //storage for '3/3' transformed integrals
    double *ptr_Coeff2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

    // #elif HAVE_CUBLAS

    // 		double *ptr_ints;
    // 		GPUerrchk(cudaMallocHost( &ptr_ints, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_ints2;
    // 		GPUerrchk(cudaMallocHost( &ptr_ints2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_Coeff;
    // 		GPUerrchk(cudaMallocHost( &ptr_Coeff, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // 		double *ptr_Coeff2;
    // 		GPUerrchk(cudaMallocHost( &ptr_Coeff2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));

    // #endif

    utility::timer timer;

    utility::timer::value_type waitTime, mixedTime, readTime, writeTime, firstReadTime,
	finalWriteTime;;

    ::cchem::ri::async async;

    //	async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);
    size_t readInitial = pe_.rank();
    //	int R = A2.NAReads_%pe_.size();
    size_t lastRead;
    if(pe_.rank() < A1.NAReads_){
	async.get(A1.AReadPos_[pe_.rank()].first,
		  A1.AReadPos_[pe_.rank()].second,
		  *(A1.Array_), ptr_ints2);
    }else{//(pe_.rank() < A1.NAReads_)
	goto nothingToDo;
    }//(pe_.rank() < A1.NAReads_)

    timer.reset();
    async.wait(); firstReadTime += timer;

#pragma omp parallel
    {
	    
	//	    for(size_t read = 0; read < A1.NAReads_; read++ ){
	for(size_t read = readInitial; read < A1.NAReads_; read += pe_.size() ){

#pragma omp single
	    {
		timer.reset();
		async.wait(); mixedTime += timer;

		std::swap(ptr_ints,ptr_ints2);
		lastRead = read;
		if(A1.NAReads_ > 1){
			
		    //if(read == 0){
		    if(read == pe_.rank() && read+pe_.size() < A1.NAReads_){
			// async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
			// 	      *(A1.Array_), ptr_ints2);
			async.get(A1.AReadPos_[read+pe_.size()].first, 
				  A1.AReadPos_[read+pe_.size()].second,
				  *(A1.Array_), ptr_ints2);

		    }else if(read == pe_.rank() && read+pe_.size() >= A1.NAReads_){

			//do nothing

		    }else{

			std::swap(ptr_Coeff,ptr_Coeff2);
			//if(read < A1.NAReads_-1){
			if(read+pe_.size() < A1.NAReads_){
				
			    // async.get_put(A1.AReadPos_[read+1].first, 
			    // 	      A1.AReadPos_[read+1].second,
			    // 	      A2.AReadPos_[read-1].first, 
			    // 	      A2.AReadPos_[read-1].second,
			    // 	      *(A1.Array_), ptr_ints2,
			    // 	      *(A2.Array_), ptr_Coeff2);
			    async.get_put(A1.AReadPos_[read+pe_.size()].first, 
					  A1.AReadPos_[read+pe_.size()].second,
					  A2.AReadPos_[read-pe_.size()].first, 
					  A2.AReadPos_[read-pe_.size()].second,
					  *(A1.Array_), ptr_ints2,
					  *(A2.Array_), ptr_Coeff2);
				
				
			}else{

			    // async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
			    // 	  *(A2.Array_), ptr_Coeff2);
			    async.put(A2.AReadPos_[read-pe_.size()].first, 
				      A2.AReadPos_[read-pe_.size()].second,
				      *(A2.Array_), ptr_Coeff2);

			    // timer.reset();
			    // async.wait(); writeTime += timer;

			}//(read < A1.NAReads_-1)
		    }//(read == 0)
		}//(A1.NAReads_ > 0)
	    }//omp single

	    Functor(ptr_ints,A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABRows_,
		    A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABCols_,
		    ptr_Coeff,A2.ABlocks_[read].first,A2.ABlocks_[read].second,A2.NABRows_,
		    A2.ABlocks_[read].first,A2.ABlocks_[read].second,A2.NABCols_);
		
	}//read

    }//omp parallel

    //NOT asynchronous
    timer.reset();
    // A2.Array_->put(ptr_Coeff, A2.AReadPos_[A2.NAReads_-1].first,
    // 	       A2.AReadPos_[A2.NAReads_-1].second);

    // A2.Array_->put(ptr_Coeff, A2.AReadPos_[A2.NAReads_-(pe_.size()+pe_.rank() -R) ].first,
    // 	       A2.AReadPos_[A2.NAReads_- (pe_.size()+pe_.rank()-R) ].second);

    A2.Array_->put(ptr_Coeff, A2.AReadPos_[ lastRead ].first,
		   A2.AReadPos_[ lastRead ].second);

    async.wait();


    finalWriteTime += timer;



 nothingToDo:;

    //	pe_.barrier();
    //	async.wait(); waitTime += timer;
    //async.wait(); writeTime += timer;
	
    //#if !HAVE_CUBLAS

    delete [] ptr_Coeff2;
    delete [] ptr_Coeff;
    delete [] ptr_ints2;
    delete [] ptr_ints;

    // #elif HAVE_CUBLAS

    // 		GPUerrchk( cudaFreeHost(ptr_Coeff2) );
    // 		GPUerrchk( cudaFreeHost(ptr_Coeff) );
    // 		GPUerrchk( cudaFreeHost(ptr_ints2) );
    // 		GPUerrchk( cudaFreeHost(ptr_ints) );

    // #endif

    waitTime = readTime + writeTime + mixedTime + firstReadTime;

    // if(pe_.rank() == 0){
    // 	std::cout << std::endl << "I/O total wait time : " << waitTime << std::endl;
    // 	std::cout << "read I/O time : " << readTime << std::endl;
    // 	std::cout << "write I/O time : " << writeTime << std::endl;
    // 	std::cout << "mixed I/O time : " << mixedTime << std::endl;
    // 	std::cout << "first read I/O time : " << firstReadTime << std::endl;
    // 	std::cout << "final write I/O time : " << finalWriteTime << std::endl;
    // }//(pe_.rank() == 0)

}//DoOpAsync_R_W_NODE_PARALLEL






template <typename FUNCTOR >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_RRSerial(FUNCTOR &Functor){

    ArrayData &A1 = Array_[0];
    ArrayData &A2 = Array_[1];


#if !HAVE_CUBLAS
    double *ptr_a1 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    double *ptr_a1_2 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    double *ptr_a2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
    double *ptr_a2_2 = NULL;
#elif HAVE_CUBLAS
    // std::cout <<A1.MaxNAB_*A1.NABRows_*A1.NABCols_ << " " 
    // 	  <<A2.MaxNAB_*A2.NABRows_*A2.NABCols_ << std::flush << std::endl;

    //size_t bytes = A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double);

    double *ptr_a1_2 = NULL;
    double *ptr_a1 = NULL;
    GPUerrchk(cudaMallocHost( &ptr_a1, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double)) );
    //cudaHostAlloc( &ptr_a1, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double),0);
    // std::cout <<cudaMallocHost( &ptr_a1, bytes )
    // 	      <<std::flush << std::endl;


    GPUerrchk(cudaMallocHost( &ptr_a1_2, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double)) );
    //cudaHostAlloc( &ptr_a1_2, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double),0);
    // std::cout 
    // 	<< cudaMallocHost( &ptr_a1_2, bytes )
    // 	<<std::flush << std::endl;;

    double *ptr_a2 = NULL;
    GPUerrchk(cudaMallocHost( &ptr_a2, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double)) );
    //cudaHostAlloc( &ptr_a2, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double),0);

    // double *pUber;
    // GPUerrchk(cudaMallocHost( (void**)&pUber, ((A2.MaxNAB_*A2.NABRows_*A2.NABCols_)
    // 					+ (A1.MaxNAB_*A1.NABRows_*A1.NABCols_)
    // 					+ (A1.MaxNAB_*A1.NABRows_*A1.NABCols_))*sizeof(double)) );
 
    // double *ptr_a1 = pUber;
    // double *ptr_a1_2 = &pUber[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    // //double *ptr_a1_2 = &pUber[A1.MaxNAB_];

    // double *ptr_a2 = &pUber[2*A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

    double *ptr_a2_2;
#endif

		
    ::cchem::ri::async async_1;
    ::cchem::ri::async async_2;

    if(A2.NAReads_ > 1){
#if !HAVE_CUBLAS
	ptr_a2_2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
#elif HAVE_CUBLAS
	GPUerrchk(cudaMallocHost( &ptr_a2_2, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double)) );
	//cudaHostAlloc( &ptr_a2_2, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double),0);
#endif
    }else{

	async_2.get(A2.AReadPos_[0].first,
		    A2.AReadPos_[0].second,
		    *(A2.Array_), ptr_a2);
	async_2.wait();
    }


    async_1.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_a1_2);
    //async_1.wait(); 
    //async_1.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_a1);
    //async_1.wait(); 

    // std::cout << A1.AReadPos_[0].first[0] << " " 
    // 	      << A1.AReadPos_[0].first[1] << std::endl
    // 	      << A1.AReadPos_[0].second[0] << " "
    // 	      << A1.AReadPos_[0].second[1] << std::flush << std::endl;
    // std::cout << A1.NAReads_ << std::flush << std::endl;
    //    ptr_a1 = ptr_a1_2;
    //double *temp;
    //temp = ptr_a1;
    //ptr_a1 = ptr_a1_2;
    //ptr_a1_2 = temp;

    //std::swap(ptr_a1,ptr_a1_2);
    //std::swap(ptr_a1,ptr_a1_2);
    {
	for(size_t reada1 = 0; reada1 < A1.NAReads_ ; reada1++ ){

	    async_1.wait();
	    std::swap(ptr_a1,ptr_a1_2);
	    if(reada1 < A1.NAReads_-1 && A1.NAReads_ > 1)
		async_1.get(A1.AReadPos_[reada1+1].first,
			    A1.AReadPos_[reada1+1].second,
			    *(A1.Array_), ptr_a1_2);
			
	    const size_t a1_start = A1.ABlocks_[reada1].first;
	    const size_t a1_stop = A1.ABlocks_[reada1].second;
	    const size_t a1_length = a1_stop - a1_start;
			
	    if(A2.NAReads_ > 1){
		async_2.get(A2.AReadPos_[A2.read_restrict_*reada1].first,
			    A2.AReadPos_[A2.read_restrict_*reada1].second,
			    *(A2.Array_), ptr_a2_2);
	    }//if(A2.NAReads_ > 1)


	    for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){

		if(A2.NAReads_ > 1){
		    async_2.wait();
		    std::swap(ptr_a2,ptr_a2_2);
		    if(reada2 < A2.NAReads_-1 && A2.NAReads_ > 1)
			async_2.get(A2.AReadPos_[reada2+1].first,
				    A2.AReadPos_[reada2+1].second,
				    *(A2.Array_), ptr_a2_2);
		}//if(A2.NAReads_ > 1)
			    
		Functor(ptr_a1,A1.ABlocks_[reada1].first,A1.ABlocks_[reada1].second,A1.NABRows_,
			A1.ABlocks_[reada1].first,A1.ABlocks_[reada1].second,A1.NABCols_,
			ptr_a2,A2.ABlocks_[reada2].first,A2.ABlocks_[reada2].second,A2.NABRows_,
			A2.ABlocks_[reada2].first,A2.ABlocks_[reada2].second,A2.NABCols_);

		// Functor(ptr_a1_2,A1.ABlocks_[reada1].first,A1.ABlocks_[reada1].second,A1.NABRows_,
		// 	    A1.ABlocks_[reada1].first,A1.ABlocks_[reada1].second,A1.NABCols_,
		// 	    ptr_a2,A2.ABlocks_[reada2].first,A2.ABlocks_[reada2].second,A2.NABRows_,
		// 	    A2.ABlocks_[reada2].first,A2.ABlocks_[reada2].second,A2.NABCols_);
			    
	    }//reada2

			
	}//reada1
		    
    }//omp parallel

#if !HAVE_CUBLAS		
    delete [] ptr_a2;
    delete [] ptr_a1;
    delete [] ptr_a2_2;
    delete [] ptr_a1_2;
#elif HAVE_CUBLAS
    //GPUerrchk( cudaFreeHost(pUber) );
    GPUerrchk( cudaFreeHost(ptr_a2) );
    GPUerrchk( cudaFreeHost(ptr_a1) );
    if(A2.NAReads_ > 1)GPUerrchk( cudaFreeHost(ptr_a2_2) );
    GPUerrchk( cudaFreeHost(ptr_a1_2) );
#endif
		
}//DoOpAsync_RRSerial





template <typename FUNCTOR >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_RR(FUNCTOR &Functor){

    ArrayData &A1 = Array_[0];
    ArrayData &A2 = Array_[1];

    utility::timer timer;
    utility::timer::value_type readTime;
    utility::timer::value_type initTime;




#if !HAVE_CUBLAS
    //std::cout << A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) << std::endl;
    double *ptr_a1 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    double *ptr_a1_2 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    double *ptr_a2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
    double *ptr_a2_2 = NULL;
#elif HAVE_CUBLAS

    double *ptr_a1 = NULL;
    GPUerrchk(cudaMallocHost( &ptr_a1, (size_t)(A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double))) );
    if(ptr_a1 == NULL) std::cout << "problem" << std::endl;

    double *ptr_a1_2 = NULL;
    GPUerrchk(cudaMallocHost( &ptr_a1_2, (size_t)(A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double)) ));
    if(ptr_a1_2 == NULL) std::cout << "problem" << std::endl;

    double *ptr_a2 = NULL;
    GPUerrchk(cudaMallocHost( &ptr_a2, (size_t)(A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double)) ));
    if(ptr_a2 == NULL) std::cout << "problem" << std::endl;

    double *ptr_a2_2;
#endif
    ::cchem::ri::async async_1;
    ::cchem::ri::async async_2;



    //    bool yes = 1;
    //    GA_Summarize(yes);
    if(A2.NAReads_ > 1){
#if !HAVE_CUBLAS
	ptr_a2_2 = new double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
#elif HAVE_CUBLAS
 	GPUerrchk(cudaMallocHost( &ptr_a2_2, (size_t)(A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double)) ));
#endif
    }else{
	async_2.get(A2.AReadPos_[0].first,A2.AReadPos_[0].second,*(A2.Array_), ptr_a2);
	timer.reset();
	async_2.wait(); initTime += timer;
    }//(A2.NAReads_ > 1)

    bool doThis = 1;    		
    size_t readInitial = pe_.rank();
    if(pe_.rank() < A1.NAReads_){
	async_1.get(A1.AReadPos_[pe_.rank()].first,
		    A1.AReadPos_[pe_.rank()].second,
		    *(A1.Array_), ptr_a1_2);
    }else{//(pe_.rank() < A1.NAReads_)
	goto nothingToDo;
    }//(pe_.rank() < A1.NAReads_)
    

#pragma omp parallel
    {
	for(size_t reada1 = readInitial; reada1 < A1.NAReads_ ; reada1 += pe_.size() ){
    

#pragma omp single
	    {

		timer.reset();
		async_1.wait(); initTime += timer;
		std::swap(ptr_a1,ptr_a1_2);
		if(reada1+pe_.size() < A1.NAReads_ && A1.NAReads_ > 1)
		    async_1.get(A1.AReadPos_[reada1+pe_.size()].first,
				A1.AReadPos_[reada1+pe_.size()].second,
				*(A1.Array_), ptr_a1_2);
	    }//omp single
	    

#pragma omp single
	    {

		//if(A2.NAReads_ > 1 && A1.Array_ == A2.Array_ && A2.read_restrict_){
		if(A2.NAReads_ > 1 && A1.Array_ == A2.Array_ && A2.loopRestrict_){
		   

		    size_t bytes = (A1.AReadPos_[reada1].first[0]-A1.AReadPos_[reada1].second[0])
			*(A1.AReadPos_[reada1].first[1]-A1.AReadPos_[reada1].second[1])*sizeof(double);
		    // size_t n = (A1.AReadPos_[reada1].first[0]-A1.AReadPos_[reada1].second[0])
		    // 	*(A1.AReadPos_[reada1].first[1]-A1.AReadPos_[reada1].second[1]);
		    // if(pe_.rank()==0) std::cout << "short cut " << bytes << std::flush << std::endl;
		    // if(pe_.rank()==0)std::cout << A1.AReadPos_[reada1].first[0] << " " 
		    // 	      << A1.AReadPos_[reada1].first[1] << " " 
		    // 	      << A1.AReadPos_[reada1].second[1] << " " 
		    // 	      << A1.AReadPos_[reada1].second[0] << std::endl;

		    //cblas_dcopy(1, ptr_a1, 1, ptr_a2_2, 1);

		    
		    // async_2.get(A2.AReadPos_[A2.read_restrict_*reada1].first,
		    // 		A2.AReadPos_[A2.read_restrict_*reada1].second,
		    // 		*(A1.Array_), ptr_a2_2);

		    // async_2.wait();
		    //memset(ptr_a2_2, 0, bytes);
		    doThis = 0;
		    memcpy(ptr_a2_2, ptr_a1, bytes);



		}else if(A2.NAReads_ > 1){
		    //if(A2.NAReads_ > 1){
		    async_2.get(A2.AReadPos_[A2.loopRestrict_*reada1].first,
				A2.AReadPos_[A2.loopRestrict_*reada1].second,
				*(A2.Array_), ptr_a2_2);
		}//if(A2.NAReads_ > 1)

	    }//omp single
	    
	    //for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){
	    for(size_t reada2 = A2.loopRestrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){


		bool restriction = 0;
		if(A2.read_restrict_ && reada2 == reada1) restriction = 1;
		
#pragma omp single
		{
		    if(A2.NAReads_ > 1){

			timer.reset();
			if(doThis)
			    async_2.wait();
			doThis = 1;
			readTime += timer;
			std::swap(ptr_a2,ptr_a2_2);
			if(reada2 < A2.NAReads_-1 && A2.NAReads_ > 1)
			    async_2.get(A2.AReadPos_[reada2+1].first,
					A2.AReadPos_[reada2+1].second,
					*(A2.Array_), ptr_a2_2);

		    }//if(A2.NAReads_ > 1)
		}//omp single

		Functor(ptr_a1,A1.ABlocks_[reada1].first,A1.ABlocks_[reada1].second,A1.NABRows_,
			A1.ABlocks_[reada1].first,A1.ABlocks_[reada1].second,A1.NABCols_,
			ptr_a2,A2.ABlocks_[reada2].first,A2.ABlocks_[reada2].second,A2.NABRows_,
			A2.ABlocks_[reada2].first,A2.ABlocks_[reada2].second,A2.NABCols_);

#pragma omp barrier
		
	    }//reada2
	    
	}//reada1
	
    }//omp parallel
    

 nothingToDo:;
    
    // if(pe_.rank() == 0)
    // 	std::cout << std::endl
    // 		  << "initial I/O wait time "
    // 		  << initTime
    // 		  << std::endl;

    // if(pe_.rank() == 0)    
    // 	std::cout << "other I/O wait time "
    // 		  << readTime
    // 		  << std::endl;
    
    
#if !HAVE_CUBLAS
    delete[] ptr_a2;
    delete[] ptr_a1;
    delete[] ptr_a2_2;
    delete[] ptr_a1_2;
#elif HAVE_CUBLAS
    GPUerrchk( cudaFreeHost(ptr_a2) );
    GPUerrchk( cudaFreeHost(ptr_a1) );
    if(A2.NAReads_ > 1)GPUerrchk( cudaFreeHost(ptr_a2_2) );
    GPUerrchk( cudaFreeHost(ptr_a1_2) );
#endif
}//DoOpAsync_RR








template< typename OP >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_R(OP &Functor) {

    ArrayData &A1 = Array_[0];
    //	std::cout << A1.NABRows_ << " " << A1.NABCols_ << std::endl;

#if !HAVE_CUBLAS
    //storage for '2/3' transformed integrals
    void * aints = NULL;
    posix_memalign(&aints, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
    double *ptr_ints = new(aints) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
	    
    //storage for '2/3' transformed integrals
    void * aints2 = NULL;
    posix_memalign(&aints2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
    double *ptr_ints2 = new(aints2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

#elif HAVE_CUBLAS
    double *ptr_ints;
    double *ptr_ints2;
    GPUerrchk(cudaMallocHost( &ptr_ints, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));
    GPUerrchk(cudaMallocHost( &ptr_ints2, (size_t)A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) ));
#endif


    ::cchem::ri::async async;
    async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

#pragma omp parallel
    {
	for(size_t read = 0; read < A1.NAReads_; read++ ){

#pragma omp single
	    {
		async.wait();
		std::swap(ptr_ints,ptr_ints2);

		if(A1.NAReads_ > 1 && read < A1.NAReads_-1){
		    async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
			      *(A1.Array_), ptr_ints2);
		}//(A1.NAReads_ > 0)
	    }//omp single

#pragma omp barrier

	    Functor(ptr_ints,A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABRows_,
		    A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABCols_);

	}//read

    }//omp parallel

#if !HAVE_CUBLAS
    delete [] ptr_ints2;
    delete [] ptr_ints;
#elif HAVE_CUBLAS
    GPUerrchk( cudaFreeHost(ptr_ints2) );
    GPUerrchk( cudaFreeHost(ptr_ints) );
#endif

}//DoOpAsync_R






template< typename OP >
void ::cchem::rimp2_gradient::detail::DFAsync::DoOpAsync_R_NODE_PARALLEL(OP &Functor) {

    utility::timer timer;
    utility::timer timer2;
    utility::timer::value_type readTime;
    timer2.reset();

    ArrayData &A1 = Array_[0];

    std::vector< size_t > myBlock;
    std::string gaString ("Array::GA");
    //if(pe_.rank()==0)std::cout<< "Array Type: " << A1.arrayType_ << std::flush << std::endl;
    int* lo = new int[2];
    int* hi = new int[2];

    if(A1.arrayType_.find(gaString) != std::string::npos){
    	//if(pe_.rank()==0)std::cout << "This is a GA array" << std::flush << std::endl;
    	A1.Array_->getArrayDistribution(pe_.rank(),lo, hi);
	// if(pe_.rank() == 0)std::cout << pe_.rank() << " " 
	// 	  << lo[0] << " " 
	// 	  << hi[0] << " " 
	// 	  << lo[1] << " " 
	// 	  << hi[1] << " " 
	// 	  << std::flush << std::endl;
	// if(pe_.rank() == 0)std::cout << A1.NAB_ << " " <<A1.NAColB_ << std::flush << std::endl;

	if( A1.NAB_ > 1){
	    
	}else{

	    for(size_t i = 0; i < A1.NAReads_; i++){
		if(A1.AReadPos_[i].first[1] >= lo[0] && 
		   A1.AReadPos_[i].first[1] <= hi[0])
		    myBlock.push_back(i);

		// if(A1.AReadPos_[i].first[1] >= lo[0] && 
		//    A1.AReadPos_[i].first[1] < hi[0])
		//     std::cout << A1.AReadPos_[i].first[1] << " "
		// 	      << lo[0] << " "
		// 	      << lo[1] << std::flush << std::endl;

		 // if(pe_.rank() == 0)std::cout << A1.AReadPos_[i].first[1] << " "
		 // 				 << lo[0] << " "
		 // 				 << lo[1] << std::flush << std::endl;

	    }
			       
	    // if(pe_.rank() == 3)
	    // 	for (int i = 0; i < myBlock.size(); i++)std::cout << myBlock[i] << " " << std::flush;
	    // if(pe_.rank() == 3)std::cout <<std::flush << std::endl;
	}
	
    }else{
    	//std::cout << "This is an HDF5 array" << std::flush << std::endl;
	if(pe_.size() > 1)std::cout << "This is an HDF5 array:" << std::endl 
				    << "this  part of the code is not setup to work with parallelHDF5" 
				    << __FILE__ << ":" << __LINE__ << std::endl;
	for(size_t i = 0; i < A1.NAReads_; i++)myBlock.push_back(i);


    }
		  
    delete [] lo;
    delete [] hi;
    
    // for(int i = 0; i < pe_.size(); i++){
    // 	if(pe_.rank() == i){
    // 	    for (int i = 0; i < myBlock.size(); i++)std::cout << myBlock[i] << " " << std::flush;
    // 	    std::cout <<std::flush << std::endl;
    // 	}//(pe_.rank() == i)
    // 	pe_.barrier();
    // }//i

    //#if !HAVE_CUBLAS
    //storage for '2/3' transformed integrals
    // void * aints = NULL;
    // posix_memalign(&aints, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
    // double *ptr_ints = new(aints) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    double *ptr_ints = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    
    //storage for '2/3' transformed integrals
    // void * aints2 = NULL;
    // posix_memalign(&aints2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
    // double *ptr_ints2 = new(aints2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    double *ptr_ints2 = new double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
    
    //     #elif HAVE_CUBLAS
    //     double *ptr_ints;
    //     double *ptr_ints2;
    //     GPUerrchk(cudaMallocHost( (void **)&ptr_ints,
    // 			      (size_t)(A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double)) ));
    //     GPUerrchk(cudaMallocHost( (void **)&ptr_ints2,
    // 			      (size_t)(A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double)) ));
    // #endif
    
    
    ::cchem::ri::async async;
    //async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

    

    // if(pe_.rank() < A1.NAReads_){
    // 	async.get(A1.AReadPos_[pe_.rank()].first,
    // 		  A1.AReadPos_[pe_.rank()].second,
    // 		  *(A1.Array_), ptr_ints2);
    // }else{//(pe_.rank() < A1.NAReads_)
    // 	goto nothingToDo;
    // }//(pe_.rank() < A1.NAReads_)

    if(myBlock.size()>0){
    	async.get(A1.AReadPos_[myBlock[0]].first,
    		  A1.AReadPos_[myBlock[0]].second,
    		  *(A1.Array_), ptr_ints2);
    }else{//(pe_.rank() < A1.NAReads_)
    	goto nothingToDo;
    }//(pe_.rank() < A1.NAReads_)

#pragma omp parallel
    {
	//for(size_t read = 0; read < A1.NAReads_; read++ ){
	//for(size_t read = pe_.rank(); read < A1.NAReads_; read += pe_.size() ){
	for(size_t read = 0; read < myBlock.size(); read++ ){
	    
#pragma omp single
	    {
		timer.reset();
		async.wait();
		readTime += timer;
		std::swap(ptr_ints,ptr_ints2);
		
		//if(A1.NAReads_ > 1 && read < A1.NAReads_-1){
		//async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
		//*(A1.Array_), ptr_ints2);

		// if(A1.NAReads_ > 1 && read+pe_.size() < A1.NAReads_){
		//     async.get(A1.AReadPos_[read+pe_.size()].first,
		// 	      A1.AReadPos_[read+pe_.size()].second,
		// 	      *(A1.Array_), ptr_ints2);
		// }//(A1.NAReads_ > 0)

		if(myBlock.size() > 1 && read < myBlock.size()-1 )
		    async.get(A1.AReadPos_[myBlock[read+1]].first,
			      A1.AReadPos_[myBlock[read+1]].second,
			      *(A1.Array_), ptr_ints2);
			

	    }//omp single
	    
	    // Functor(ptr_ints,A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABRows_,
	    //  	    A1.ABlocks_[read].first,A1.ABlocks_[read].second,A1.NABCols_);
	    Functor(ptr_ints,A1.ABlocks_[myBlock[read]].first,
		    A1.ABlocks_[myBlock[read]].second, A1.NABRows_,
		    A1.ABlocks_[myBlock[read]].first,
		    A1.ABlocks_[myBlock[read]].second, A1.NABCols_);

	}//read


    }//omp parallel

 nothingToDo:;






    pe_.barrier();

    //#if !HAVE_CUBLAS
    delete [] ptr_ints2;
    delete [] ptr_ints;
    // #elif HAVE_CUBLAS
    //     GPUerrchk( cudaFreeHost(ptr_ints2) );
    //     GPUerrchk( cudaFreeHost(ptr_ints) );
    // #endif

    // std::cout << pe_.rank() 
    // 	      << ": Wait I/O time: " 
    // 	      << readTime << std::flush << std::endl;

    // pe_.barrier();
    
    // std::cout << pe_.rank() 
    // 	      << ": Total time: " 
    // 	      << timer2 << std::flush << std::endl;


}//DoOpAsync_R_NODE_PARALLEL













































// namespace cchem{
// namespace rimp2_gradient{
// namespace detail{

// class RIAsync  {

// 	//    using cchem::Thread;


// 	//public:
// 	//	void get(size_t i, size_t j, size_t l,
// 	//			const Array<double> &V, double *data) {
// 	//
// 	//		thread_.reset(new boost::thread(launch, i, j, l, boost::cref(V), data));
// 	//	}
// 	//
// 	//
// 	//	void wait() {
// 	//    thread_->join();
// 	//	}
// 	//
// 	//private:
// 	//	std::auto_ptr<boost::thread> thread_;
// 	//
// 	//	static
// 	////		void launch(size_t k, size_t m, size_t n, size_t b,
// 	////				const Array<double> &V, double *data) {
// 	//	void launch(size_t i, size_t j, size_t l,
// 	//			const Array<double> &V, double *data) {
// 	////    		BOOST_PROFILE_LINE;
// 	//		utility::timer timer;
// 	//		size_t start[] = { i, 0 };
// 	//		size_t finish[] = { j, l };
// 	//		V.get(data, start, finish);
// 	//		//std::cout << "I/O: " << (N*N*B*8)/(timer*1e6) << std::endl;
// 	//	}




// 	typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;

// 	typedef std::vector < std::pair <std::vector <size_t>, std::vector <size_t> > > ReadPos;
// 	typedef std::vector < std::pair < size_t, size_t> > Blocks;

// 	struct ArrayData{

// 		Array<double> *Array_;
// 		//number of blocks , number of rows per block , number of columns per block
// 		size_t NAB_, NABRows_;
// 		size_t NAColB_, NABCols_;
// 		ReadPos	AReadPos_;
// 		Blocks ABlocks_;

// 		ReadPos	AReadColPos_;
// 		Blocks AColBlocks_;


// 		//larget block possible
// 		size_t MaxNAB_, NAReads_;

// 		size_t MaxNAColB_, NAColReads_;
// 		//		size_t restrict_;
// 		// used to restict arrays reads to depend on other arrays (user defined really)
// 		bool read_restrict_;

// 		ArrayData(Array<double> *Array, size_t NAB) :
// 			Array_(Array),
// 			NAB_(NAB),NABRows_(0),
// 			NAColB_(0),NABCols_(0),
// 			AReadPos_(),ABlocks_(),
// 			AReadColPos_(),AColBlocks_(),
// 			MaxNAB_(0),NAReads_(0),
// 			MaxNAColB_(0),NAColReads_(0),
// 			read_restrict_(0)
// 		{};

// 	};//ArrayData


// 	size_t NMB_;
// 	Array<double> *ArrayOne_;
// 	Array<double> *ArrayTwo_;
// 	Array<double> *ArrayThr_;
// 	Array<double> *ArrayFou_;
// 	std::vector<ArrayData> Array_;



// public:



// 	RIAsync(size_t NMB,
// 			Array<double> *ArrayOne, size_t &NA1B,
// 			bool read_rows, bool async) :
// 				NMB_(NMB),
// 				ArrayOne_(ArrayOne),
// 				ArrayTwo_(NULL),
// 				ArrayThr_(NULL),
// 				ArrayFou_(NULL)
// {

// 		Array_.push_back(ArrayData(ArrayOne_,NA1B));

// 		//		size_t NWorkingBytes = NMB_*1024*1024;
// 		size_t NWorkingBytes = NMB_*1000*1000;

// 		this->set_up_single_array_info(read_rows,async, NWorkingBytes,
// 				Array_[0], 0);

// 		Array_[1].read_restrict_ = 0; //useless here

// 		//		std::cout << NAB_ << std::endl;
// //		this->PrintArray(0);

// }//RIAsync



// 	RIAsync(size_t NMB,
// 			Array<double> *ArrayOne, size_t &NA1B, bool A1_access_rows,
// 			Array<double> *ArrayTwo, size_t &NA2B, bool A2_access_rows,
// 			bool async) :
// 				NMB_(NMB),
// 				ArrayOne_(ArrayOne),
// 				ArrayTwo_(ArrayTwo),
// 				ArrayThr_(NULL),
// 				ArrayFou_(NULL)
// 	{

// 		Array_.push_back(ArrayData(ArrayOne_,NA1B));
// 		Array_.push_back(ArrayData(ArrayTwo_,NA2B));

// 		//		size_t NWorkingBytes = NMB_*1024*1024;
// 		size_t NWorkingBytes = NMB_*1000*1000;

// 		this->set_up_two_array_info(A1_access_rows,async, NWorkingBytes,
// 				Array_[0], 0,Array_[1], 0);

// 		Array_[1].read_restrict_ = 1;
// 		//
// 		//		std::cout << NAB_ << std::endl;
// //		this->PrintArray(0);
// //		this->PrintArray(1);
// 	}//RIAsync







// 	RIAsync(size_t NMB,
// 			Array<double> *ArrayOne, size_t &NA1B,
// 			Array<double> *ArrayTwo, size_t &NA2B,
// 			bool read_restrict, bool async) :
// 				NMB_(NMB),
// 				ArrayOne_(ArrayOne),
// 				ArrayTwo_(ArrayTwo),
// 				ArrayThr_(NULL),
// 				ArrayFou_(NULL)
// 	{
// 		Array_.push_back(ArrayData(ArrayOne_,NA1B));
// 		Array_.push_back(ArrayData(ArrayTwo_,NA2B));

// 		//		bool async = 0;
// 		//		size_t NWorkingBytes = NMB_*1024*1024;
// 		size_t NWorkingBytes = NMB_*1000*1000;

// 		this->set_up_array_info(read_restrict,async, NWorkingBytes,
// 				Array_[0], 0,Array_[1], 0);
// 		Array_[1].read_restrict_ = read_restrict;

// 		//		this->PrintArray(0);
// 		//		this->PrintArray(1);
// 	}//RIAsync


// 	RIAsync(size_t NMB,
// 			Array<double> *ArrayOne, size_t NA1BS, size_t &NA1B,
// 			Array<double> *ArrayTwo, size_t NA2BS, size_t &NA2B,
// 			bool read_restrict, bool async, int tag) :
// 				NMB_(NMB),
// 				ArrayOne_(ArrayOne),
// 				ArrayTwo_(ArrayTwo),
// 				ArrayThr_(NULL),
// 				ArrayFou_(NULL)
// 	{
// 		std::swap(ArrayOne_,ArrayTwo_);
// 		std::swap(NA1B,NA2B);
// 		std::swap(NA1BS,NA2BS);

// 		Array_.push_back(ArrayData(ArrayOne_, (NA1B-NA1BS) ) );
// 		Array_.push_back(ArrayData(ArrayTwo_, (NA2B-NA2BS) ) );

// 		//		bool async = 0;
// 		//		size_t NWorkingBytes = NMB_*1024*1024;
// 		size_t NWorkingBytes = NMB_*1000*1000;

// 		this->set_up_array_info(read_restrict,async, NWorkingBytes,
// 				Array_[0], NA1BS,Array_[1], NA2BS);
// 		Array_[1].read_restrict_ = read_restrict;

// 		std::iter_swap(Array_.begin() + 0,Array_.begin() + 1);
// 		std::swap(ArrayOne_,ArrayTwo_);
// 		std::swap(NA1B,NA2B);
// 		std::swap(NA1BS,NA2BS);

// 		//		this->PrintArray(0);
// 		//		this->PrintArray(1);
// 	}//RIAsync

// 	RIAsync(size_t NMB,
// 			Array<double> *ArrayOne, size_t &NA1B,
// 			Array<double> *ArrayTwo, size_t &NA2B, bool read_restrict_array2,
// 			Array<double> *ArrayThr, size_t &NA3B) :
// 				NMB_(NMB),
// 				ArrayOne_(ArrayOne),
// 				ArrayTwo_(ArrayTwo),
// 				ArrayThr_(ArrayThr),
// 				ArrayFou_(NULL)
// 	{
// 		Array_.push_back(ArrayData(ArrayOne_,NA1B));
// 		Array_.push_back(ArrayData(ArrayTwo_,NA2B));
// 		Array_.push_back(ArrayData(ArrayThr_,NA3B));

// 		bool async = 0;

// 		//		size_t NWorkingBytes = NMB_*1024*1024;
// 		size_t NWorkingBytes = NMB_*1000*1000;
// 		this->set_up_array_info(read_restrict_array2,async, NWorkingBytes,
// 				Array_[0], 0,Array_[1], 0);
// 		Array_[1].read_restrict_ = read_restrict_array2;
// 		//old-not work with sym		this->set_up_array_info(restrict,async,
// 		//				Array_[2],Array_[3]);

// 		this->set_up_unrestricted_array_info(async, NWorkingBytes,
// 				Array_[0],Array_[2]);
// 		Array_[2].read_restrict_ = 0;
// 		//		this->PrintArray(0);
// 		//		this->PrintArray(1);
// 		//		this->PrintArray(2);
// 		//		this->PrintArray(3);

// 	}//RIAsync


// 	//	RIAsync(size_t NMB,
// 	//			Array<double> *ArrayOne, size_t &NA1B,
// 	//			Array<double> *ArrayTwo, size_t &NA2B,
// 	//			Array<double> *ArrayThr, size_t &NA3B,
// 	//			Array<double> *ArrayFou, size_t &NA4B) :
// 	//				NMB_(NMB),
// 	//				ArrayOne_(ArrayOne),
// 	//				ArrayTwo_(ArrayTwo),
// 	//				ArrayThr_(ArrayThr),
// 	//				ArrayFou_(ArrayFou)
// 	//	{
// 	//		Array_.push_back(ArrayData(ArrayOne_,NA1B));
// 	//		Array_.push_back(ArrayData(ArrayTwo_,NA2B));
// 	//		Array_.push_back(ArrayData(ArrayThr_,NA1B));
// 	//		Array_.push_back(ArrayData(ArrayFou_,NA4B));
// 	//
// 	//		bool restrict = 0;
// 	//		bool async = 0;
// 	//
// 	//		this->set_up_array_info(restrict,async,
// 	//				Array_[0],Array_[1]);
// 	//
// 	//		this->set_up_array_info(restrict,async,
// 	//				Array_[2],Array_[3]);
// 	//
// 	//		//		this->PrintArray(0);
// 	//		//		this->PrintArray(1);
// 	//		//		this->PrintArray(2);
// 	//		//		this->PrintArray(3);
// 	//	}//RIAsync

// 	void set_up_array_info(bool &restrict, bool &async, size_t & NWorkingBytes,
// 			ArrayData &A1,  size_t NAB1S, ArrayData &A2,  size_t NAB2S);


// 	void set_up_single_array_info(bool &row, bool &async, size_t & NWorkingBytes,
// 			ArrayData &A1,  size_t NAB1S);

// 	void set_up_two_array_info(bool &row, bool &async, size_t & NWorkingBytes,
// 			ArrayData &A1,  size_t NAB1S, ArrayData &A2,  size_t NAB2S);

// 	//	void set_up_array_info(bool &restrict, bool &async, size_t & NWorkingBytes,
// 	//			ArrayData &A1,  size_t NAB1S, ArrayData &A2,  size_t NAB2S){
// 	////		ArrayData &A1,  ArrayData &A2){
// 	//
// 	//		//we need to see how many blocks of A1 and A2 can
// 	//		// be loaded into the buffers at the same time
// 	////		size_t NBytes = NMB_*1024*1024;
// 	//		size_t NBytes = NWorkingBytes;
// 	//
// 	//		A1.NABRows_ = A1.Array_->shape()[0]/(A1.NAB_+NAB1S);
// 	//		A1.NABCols_ = A1.Array_->shape()[1];
// 	//		A2.NABRows_ = A2.Array_->shape()[0]/(A2.NAB_+NAB2S);
// 	//		A2.NABCols_ = A2.Array_->shape()[1];
// 	//
// 	//		size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);
// 	//		size_t A2BlockSize = A2.NABRows_*A2.NABCols_*sizeof(double);
// 	//
// 	//		if(restrict){
// 	//
// 	//			//we'll need at least one block from A1 and A2
// 	//			A1.MaxNAB_=1;
// 	//
// 	//			NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
// 	//			//see how many CijQ/BARE_IA blocks we can fit into the buffer
// 	//			while(NBytes > (1+async)*(A1BlockSize + A2BlockSize) ){
// 	//				A1.MaxNAB_++;
// 	//				NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
// 	//				if(A1.MaxNAB_ == A1.NAB_)break;
// 	//			}
// 	//			A2.MaxNAB_ = A1.MaxNAB_;
// 	//
// 	//			//figure out how many reads will need to occur
// 	//			A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 	//			A2.NAReads_ = A1.NAReads_;
// 	//
// 	//			this->set_up_block_info(A1,  NAB1S);
// 	//
// 	//			//restrict: A1 to be the same as A2
// 	//			A2.AReadPos_ = A1.AReadPos_;
// 	//			A2.ABlocks_ = A1.ABlocks_;
// 	//
// 	//			//			A1.restrict_ = 1;
// 	//			//			A2.restrict_ = 1;
// 	//
// 	//		}else{//(restrict)
// 	//
// 	//			//we'll need at least one block from A1 and A2
// 	//			A1.MaxNAB_=1;
// 	//			A2.MaxNAB_=1;
// 	//
// 	//			NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
// 	//			//see how many A1.MaxNAB_ blocks we can fit into the buffer
// 	//			while(NBytes > (1+async)*A1BlockSize ){
// 	//				A1.MaxNAB_++;
// 	//				NBytes -= (1+async)*A1BlockSize;
// 	//				if(A1.MaxNAB_ == A1.NAB_)break;
// 	//			}
// 	//
// 	//			//see how many A2.MaxNAB_ blocks we can fit into the buffer
// 	//			while(NBytes > (1+async)*A2BlockSize ){
// 	//				A2.MaxNAB_++;
// 	//				NBytes -= (1+async)*A2BlockSize;
// 	//				if(A2.MaxNAB_ == A2.NAB_)break;
// 	//			}
// 	//
// 	//			//figure out how many reads will need to occur
// 	//			A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 	//			A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);
// 	//
// 	//			//try to balance out read size
// 	//			A1.MaxNAB_ = A1.NAB_/A1.NAReads_ + (A1.NAB_%A1.NAReads_>0);
// 	//			A2.MaxNAB_ = A2.NAB_/A2.NAReads_ + (A2.NAB_%A2.NAReads_>0);
// 	//
// 	//			this->set_up_block_info(A1, NAB1S);
// 	//			this->set_up_block_info(A2, NAB2S);
// 	//
// 	//			//			A1.restrict_ = 0;
// 	//			//			A2.restrict_ = 0;
// 	////			std::cout << NAB1S << " "<< NAB2S << std::endl;
// 	//		}//(restrict)
// 	//
// 	//	}//set_up_array_info

// 	void set_up_unrestricted_array_info(bool &async, size_t &NWorkingBytes,
// 			ArrayData &A1, ArrayData &A2);

// 	//	void set_up_unrestricted_array_info(bool &async, size_t &NWorkingBytes,
// 	//			ArrayData &A1, ArrayData &A2) {
// 	//
// 	//		//we need to see how many blocks of A1 and A2 can
// 	//		// be loaded into the buffers at the same time
// 	////		size_t NBytes = NMB_*1024*1024;
// 	//		size_t NBytes = NWorkingBytes;
// 	//
// 	//		A2.NABRows_ = A2.Array_->shape()[0]/A2.NAB_;
// 	//		A2.NABCols_ = A2.Array_->shape()[1];
// 	//
// 	//		size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);
// 	//		size_t A2BlockSize = A2.NABRows_*A2.NABCols_*sizeof(double);
// 	//
// 	//		//		if(restrict){
// 	//
// 	//		//we'll need at least one block from A2 (A1 is already known)
// 	//		A2.MaxNAB_=1;
// 	//
// 	//		NBytes -= A1BlockSize*A1.MaxNAB_;
// 	//		NBytes -= (A2BlockSize);
// 	//		//see how many A1,A2,A3 blocks we can fit into the buffer
// 	//		while(NBytes > (A2BlockSize) ){
// 	//			A2.MaxNAB_++;
// 	//			NBytes -= (A2BlockSize);
// 	//			if(A2.MaxNAB_ == A2.NAB_)break;
// 	//		}
// 	//
// 	//		//figure out how many reads will need to occur
// 	//		A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);
// 	//
// 	//		this->set_up_block_info(A2, 0);
// 	//
// 	//		//this just means that the buffers are restrict to NMB
// 	//		//			A2.restrict_ = 0;
// 	//
// 	//
// 	//		//		}else{//(restrict)
// 	//		//
// 	//		//			//we'll need at least one block from A1 and A2
// 	//		//			A1.MaxNAB_=1;
// 	//		//			A2.MaxNAB_=1;
// 	//		//
// 	//		//			NBytes -= (A1BlockSize + A2BlockSize);
// 	//		//			//see how many A1.MaxNAB_ blocks we can fit into the buffer
// 	//		//			while(NBytes > A1BlockSize ){
// 	//		//				A1.MaxNAB_++;
// 	//		//				NBytes -= A1BlockSize;
// 	//		//				if(A1.MaxNAB_ == A1.NAB_)break;
// 	//		//			}
// 	//		//
// 	//		//			//see how many A2.MaxNAB_ blocks we can fit into the buffer
// 	//		//			while(NBytes > A2BlockSize ){
// 	//		//				A2.MaxNAB_++;
// 	//		//				NBytes -= A2BlockSize;
// 	//		//				if(A2.MaxNAB_ == A2.NAB_)break;
// 	//		//			}
// 	//		//
// 	//		//			//figure out how many reads will need to occur
// 	//		//			A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 	//		//			A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);
// 	//		//
// 	//		//			this->set_up_block_info(A1);
// 	//		//			this->set_up_block_info(A2);
// 	//		//
// 	//		//			A1.restrict_ = 0;
// 	//		//			A2.restrict_ = 0;
// 	//		//
// 	//		//		}//(restrict)
// 	//
// 	//	}//set_up_3array_info


// 	void set_up_block_info(ArrayData &A1,  size_t NAB1S);

// 	void set_up_block_info_cols(ArrayData &A1,  size_t NAB1S);

// 	//	void set_up_block_info(ArrayData &A1,  size_t NAB1S){
// 	//
// 	//		if(A1.NAReads_ == 1){
// 	//			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 	//			temp.first.push_back(NAB1S*A1.NABRows_); temp.first.push_back(0);
// 	//			temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABRows_); temp.second.push_back(A1.NABCols_);
// 	//			A1.AReadPos_.push_back(temp);
// 	//			A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S,NAB1S+A1.NAB_) );
// 	//		}else{
// 	//			for(int i =0; i<A1.NAReads_-1; i++){
// 	//				std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 	//				temp.first.push_back( (NAB1S+i*A1.MaxNAB_)*A1.NABRows_); temp.first.push_back(0);
// 	//				temp.second.push_back( (NAB1S+(i+1)*A1.MaxNAB_)*A1.NABRows_); temp.second.push_back(A1.NABCols_);
// 	//				A1.AReadPos_.push_back(temp);
// 	//				A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+i*A1.MaxNAB_ , NAB1S+(i+1)*A1.MaxNAB_) );
// 	//			}//i
// 	//			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 	//			temp.first.push_back( (NAB1S+(A1.NAReads_-1)*A1.MaxNAB_) *A1.NABRows_); temp.first.push_back(0);
// 	//			temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABRows_ ); temp.second.push_back(A1.NABCols_);
// 	//			A1.AReadPos_.push_back(temp);
// 	//			A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+(A1.NAReads_-1)*A1.MaxNAB_ , NAB1S+A1.NAB_) );
// 	//		}//(A1.NAReads_ == 1)
// 	//
// 	//	}//set_up_block_info

// 	void PrintArray(int i);

// 	//	void PrintArray(int i){
// 	//
// 	//		ArrayData &A = Array_[i];
// 	//		std::cout << std::endl << "printing information for Array_[" << i << "]:" << std::endl ;
// 	//		std::cout << "  full array "<< i<<" looks like " << A.Array_->shape()[0] << " x " << A.Array_->shape()[1] << std::endl;
// 	//		std::cout << "  single array block size: " << A.NABRows_  << " x " <<  A.NABCols_  << std::endl;
// 	//		std::cout << "  need to perform " << A.NAReads_ << " NAReads_" << std::endl;
// 	//		std::cout << "  MaxNAB_ " << A.MaxNAB_ << std::endl;
// 	//		std::cout << "  array reads look like: " << std::endl;
// 	//		for(int i = 0; i < A.NAReads_; i++){
// 	//			std::cout << "    " <<"["<< i <<"]"<< "   { " << A.AReadPos_[i].first[0] << " " << A.AReadPos_[i].first[1] << " }"
// 	//			<<" --> { "<< A.AReadPos_[i].second[0] << " " << A.AReadPos_[i].second[1] << " }" << std::endl;
// 	//		}
// 	//		std::cout << "  array blocks range like: " << std::endl;
// 	//		for(int i = 0; i < A.NAReads_; i++){
// 	//			std::cout << "    " <<"["<< i  <<"]"<< "   " << A.ABlocks_[i].first << " : " << A.ABlocks_[i].second << std::endl;
// 	//		}
// 	//		std::cout << std::endl;
// 	//	}//Print

// 	//	template< typename OP >
// 	//	void DoOperation(OP &Functor);

// 	//	template< typename OP >
// 	//	void DoOperation(OP Functor) {
// 	//
// 	//		ArrayData &A1 = Array_[0];
// 	//		ArrayData &A2 = Array_[1];
// 	//
// 	//		//		::cchem::ri::async async;
// 	//
// 	//		//storage for '2/3' transformed integrals
// 	//		void * aints = NULL;
// 	//		posix_memalign(&aints, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 	//		double *ptr_ints = new(aints) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
// 	//
// 	//		//storage for '3/3' transformed integrals
// 	//		void * aCoeff = NULL;
// 	//		posix_memalign(&aCoeff, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 	//		double *ptr_Coeff = new(aCoeff) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
// 	//
// 	//#pragma omp parallel
// 	//		{
// 	//			for(size_t read = 0; read < A1.NAReads_; read++ ){
// 	//
// 	//#pragma omp single
// 	//				A1.Array_->get(ptr_ints, A1.AReadPos_[read].first, A1.AReadPos_[read].second);
// 	//				size_t b_start = A1.ABlocks_[read].first;
// 	//				size_t b_stop = A1.ABlocks_[read].second;
// 	//				size_t b_length = b_stop - b_start;
// 	//
// 	//				MapMatrixXd Integrals(ptr_ints,b_length*A1.NABRows_,A1.NABCols_);   // '2/3' transformed integrals (ij|Q)
// 	//				MapMatrixXd Coeff(ptr_Coeff,b_length*A2.NABRows_,A2.NABCols_);   // '3/3' transformed integrals (ij|Q)
// 	//
// 	//				Functor(Coeff, b_start, b_stop, A1.NABRows_, A1.NABCols_,
// 	//						Integrals, b_start, b_stop, A1.NABRows_, A1.NABCols_);
// 	//
// 	//#pragma omp single
// 	//				A2.Array_->put(ptr_Coeff, A2.AReadPos_[read].first, A2.AReadPos_[read].second);
// 	//
// 	//			}//read
// 	//
// 	//		}//omp parallel
// 	//
// 	//		delete ptr_Coeff;
// 	//		delete ptr_ints;
// 	//
// 	//	}//DoOperation






// 	template< typename OP >
// 	void DoOperationAsync(OP Functor) {
// 		//luke
// 		ArrayData &A1 = Array_[0];
// 		ArrayData &A2 = Array_[1];

// 		//storage for '2/3' transformed integrals
// 		void * aints = NULL;
// 		posix_memalign(&aints, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_ints = new(aints) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		//storage for '2/3' transformed integrals
// 		void * aints2 = NULL;
// 		posix_memalign(&aints2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_ints2 = new(aints2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		//storage for '3/3' transformed integrals
// 		void * aCoeff = NULL;
// 		posix_memalign(&aCoeff, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		double *ptr_Coeff = new(aCoeff) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		//storage for '3/3' transformed integrals
// 		void * aCoeff2 = NULL;
// 		posix_memalign(&aCoeff2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		double *ptr_Coeff2 = new(aCoeff2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		//		std::cout << "this far" << std::flush << std::endl;
// 		::cchem::ri::async async;
// 		async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

// #pragma omp parallel
// 		{

// 			//#pragma omp single
// 			//			async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

// 			//			std::cout << "start reads for A1" << std::endl;
// 			for(size_t read = 0; read < A1.NAReads_; read++ ){
// 				//				std::cout << "done with read " << read << " of " << A1.NAReads_ << std::endl;
// 				//				#pragma omp single
// 				//								A1.Array_->get(ptr_ints, A1.AReadPos_[read].first, A1.AReadPos_[read].second);

// #pragma omp single
// 				{
// 					async.wait();
// 					std::swap(ptr_ints,ptr_ints2);

// 					if(A1.NAReads_ > 1){

// 						if(read == 0){
// 							async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
// 									*(A1.Array_), ptr_ints2);
// 						}else{
// 							std::swap(ptr_Coeff,ptr_Coeff2);
// 							if(read < A1.NAReads_-1){
// 								async.get_put(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
// 										A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
// 										*(A1.Array_), ptr_ints2,
// 										*(A2.Array_), ptr_Coeff2);
// 							}else{
// 								async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
// 										*(A2.Array_), ptr_Coeff2);
// 							}//(read < A1.NAReads_-1)
// 						}//(read == 0)
// 					}//(A1.NAReads_ > 0)
// 				}//omp master

// #pragma omp barrier

// 				size_t b_start = A1.ABlocks_[read].first;
// 				size_t b_stop = A1.ABlocks_[read].second;
// 				size_t b_length = b_stop - b_start;

// 				MapMatrixXd Integrals(ptr_ints,b_length*A1.NABRows_,A1.NABCols_);   // '2/3' transformed integrals (ij|Q)
// 				MapMatrixXd Coeff(ptr_Coeff,b_length*A2.NABRows_,A2.NABCols_);   // '3/3' transformed integrals (ij|Q)

// 				Functor(Coeff, b_start, b_stop, A1.NABRows_, A1.NABCols_,
// 						Integrals, b_start, b_stop, A1.NABRows_, A1.NABCols_);

// 				//				Functor(ptr_coeff, b_start, b_stop, A1.NABRows_, A1.NABCols_,
// 				//						ptr_ints, b_start, b_stop, A1.NABRows_, A1.NABCols_);




// 				//				#pragma omp single
// 				//								A2.Array_->put(ptr_Coeff, A2.AReadPos_[read].first, A2.AReadPos_[read].second);

// 				//#pragma omp barrier

// 				//#pragma omp single
// 				//				if(read == A1.NAReads_-1){
// 				//					async.wait();
// 				//					A2.Array_->put(ptr_Coeff, A2.AReadPos_[read].first, A2.AReadPos_[read].second);
// 				//				}//(read == A1.NAReads_-1)

// 			}//read

// 		}//omp parallel

// 		async.wait();
// 		A2.Array_->put(ptr_Coeff, A2.AReadPos_[A2.NAReads_-1].first,
// 				A2.AReadPos_[A2.NAReads_-1].second);

// 		delete ptr_Coeff2;
// 		delete ptr_Coeff;
// 		delete ptr_ints2;
// 		delete ptr_ints;
// 		//std::cout << "leaving " << std::flush << std::endl;

// 	}//DoOperationAsync




// 	//	template <typename FUNCTOR_ONE, typename FUNCTOR_TWO >
// 	//	void DoOperation2(FUNCTOR_ONE FunctorOne, FUNCTOR_TWO FunctorTwo){
// 	//
// 	//		ArrayData &A1 = Array_[0];
// 	//		ArrayData &A2 = Array_[1];
// 	//		ArrayData &A3 = Array_[2];
// 	//
// 	//		void * a1 = NULL;
// 	//		posix_memalign(&a1, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 	//		double *ptr_a1 = new(a1) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
// 	//
// 	//		size_t largest_buffer =
// 	//				std::max(A2.MaxNAB_*A2.NABRows_*A2.NABCols_,A3.MaxNAB_*A3.NABRows_*A3.NABCols_);
// 	//		void * a2 = NULL;
// 	//		posix_memalign(&a2, 16, largest_buffer*sizeof(double) );
// 	//		double *ptr_a2 = new(a2) double[largest_buffer];
// 	//
// 	//#pragma omp parallel
// 	//		{
// 	//			size_t largest_eri = std::max(A1.NABRows_*A2.NABRows_,A1.NABRows_*A3.NABRows_);
// 	//			void * aeri = NULL;
// 	//			posix_memalign(&aeri, 16, largest_eri*sizeof(double) );
// 	//			double *ptr_eri = new(aeri) double[largest_eri];
// 	//
// 	//			MapMatrixXd eri_vooo(ptr_eri,A1.NABRows_,A2.NABRows_);
// 	//			MapMatrixXd eri_vvvo(ptr_eri,A1.NABRows_,A3.NABRows_);
// 	//
// 	//			for(size_t reada1 = 0; reada1 < A1.NAReads_ ; reada1++ ){
// 	//#pragma omp single
// 	//				A1.Array_->get(ptr_a1, A1.AReadPos_[reada1].first, A1.AReadPos_[reada1].second);
// 	//				size_t a1_start = A1.ABlocks_[reada1].first;
// 	//				size_t a1_stop = A1.ABlocks_[reada1].second;
// 	//				size_t a1_length = a1_stop - a1_start;
// 	//				MapMatrixXd subA1(ptr_a1,a1_length*A1.NABRows_,A1.NABCols_);
// 	//
// 	//				for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){
// 	//					bool loop_restriction = 0;
// 	//					if(A2.read_restrict_ && reada2 == reada1)loop_restriction = 1;
// 	//					//				for(size_t reada2 = 0; reada2 < A2.NAReads_ ; reada2++ ){
// 	//					//					bool loop_restriction = 0;
// 	//
// 	//#pragma omp single
// 	//					A2.Array_->get(ptr_a2, A2.AReadPos_[reada2].first, A2.AReadPos_[reada2].second);
// 	//					size_t a2_start = A2.ABlocks_[reada2].first;
// 	//					size_t a2_stop = A2.ABlocks_[reada2].second;
// 	//					size_t a2_length = a2_stop - a2_start;
// 	//					MapMatrixXd subA2(ptr_a2,a2_length*A2.NABRows_,A2.NABCols_);
// 	//
// 	//					FunctorOne(eri_vooo,
// 	//							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 	//							subA2,a2_start,a2_stop, A2.NABRows_,A2.NABCols_,
// 	//							loop_restriction);
// 	//				}//reada2
// 	//
// 	//				//				for(size_t reada3 = A3.restrict_*reada1; reada3 < A3.NAReads_ ; reada3++ ){
// 	//				//					bool loop_restriction = 0;
// 	//				//					if(A3.restrict_ && reada3 == reada1)loop_restriction = 1;
// 	//				for(size_t reada3 = 0; reada3 < A3.NAReads_ ; reada3++ ){
// 	//					bool loop_restriction = 0;
// 	//
// 	//#pragma omp single
// 	//					A3.Array_->get(ptr_a2, A3.AReadPos_[reada3].first, A3.AReadPos_[reada3].second);
// 	//					size_t a3_start = A3.ABlocks_[reada3].first;
// 	//					size_t a3_stop = A3.ABlocks_[reada3].second;
// 	//					size_t a3_length = a3_stop - a3_start;
// 	//					MapMatrixXd subA3(ptr_a2,a3_length*A3.NABRows_,A3.NABCols_);
// 	//
// 	//					//#pragma omp barrier
// 	//					FunctorTwo(eri_vvvo,
// 	//							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 	//							subA3,a3_start,a3_stop, A3.NABRows_,A3.NABCols_,
// 	//							loop_restriction);
// 	//				}//reada3
// 	//
// 	//			}//reada1
// 	//
// 	//
// 	//
// 	//			delete ptr_eri;
// 	//
// 	//		}//omp parallel
// 	//
// 	//		delete ptr_a2;
// 	//		delete ptr_a1;
// 	//
// 	//	}//DoOperatrion2










// 	template <typename FUNCTOR_ONE, typename FUNCTOR_TWO >
// 	void DoOperation2Async(FUNCTOR_ONE FunctorOne,FUNCTOR_TWO FunctorTwo){

// 		ArrayData &A1 = Array_[0];
// 		ArrayData &A2 = Array_[1];
// 		ArrayData &A3 = Array_[2];

// 		void * a1 = NULL;
// 		posix_memalign(&a1, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_a1 = new(a1) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		void * a1_2 = NULL;
// 		posix_memalign(&a1_2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_a1_2 = new(a1_2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		//		void * a2 = NULL;
// 		//		posix_memalign(&a2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		//		double *ptr_a2 = new(a2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		//		void * a3 = NULL;
// 		//		posix_memalign(&a3, 16, A3.MaxNAB_*A3.NABRows_*A3.NABCols_*sizeof(double) );
// 		//		double *ptr_a3 = new(a3) double[A3.MaxNAB_*A3.NABRows_*A3.NABCols_];

// 		size_t largest_buffer =
// 				std::max(A2.MaxNAB_*A2.NABRows_*A2.NABCols_,A3.MaxNAB_*A3.NABRows_*A3.NABCols_);

// 		void * a23 = NULL;
// 		posix_memalign(&a23, 16, largest_buffer*sizeof(double) );
// 		double *ptr_a23 = new(a23) double[largest_buffer];

// 		void * a23_2 = NULL;
// 		posix_memalign(&a23_2, 16, largest_buffer*sizeof(double) );
// 		double *ptr_a23_2 = new(a23_2) double[largest_buffer];


// 		::cchem::ri::async async_1;
// 		::cchem::ri::async async_2;
// 		::cchem::ri::async async_3;

// 		async_1.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_a1_2);

// #pragma omp parallel
// 		{
// 			size_t largest_eri = std::max(A1.NABRows_*A2.NABRows_,A1.NABRows_*A3.NABRows_);
// 			void * aeri = NULL;
// 			posix_memalign(&aeri, 16, largest_eri*sizeof(double) );
// 			double *ptr_eri = new(aeri) double[largest_eri];

// 			MapMatrixXd eri_vooo(ptr_eri,A1.NABRows_,A2.NABRows_);
// 			MapMatrixXd eri_vvvo(ptr_eri,A1.NABRows_,A3.NABRows_);

// 			for(size_t reada1 = 0; reada1 < A1.NAReads_ ; reada1++ ){
// 				//#pragma omp single
// 				//				A1.Array_->get(ptr_a1, A1.AReadPos_[reada1].first, A1.AReadPos_[reada1].second);

// #pragma omp single
// 				{
// 					async_1.wait();
// 					std::swap(ptr_a1,ptr_a1_2);
// 					if(reada1 < A1.NAReads_-1 && A1.NAReads_ > 1)
// 						async_1.get(A1.AReadPos_[reada1+1].first,
// 								A1.AReadPos_[reada1+1].second,
// 								*(A1.Array_), ptr_a1_2);
// 				}//omp single

// 				size_t a1_start = A1.ABlocks_[reada1].first;
// 				size_t a1_stop = A1.ABlocks_[reada1].second;
// 				size_t a1_length = a1_stop - a1_start;
// 				MapMatrixXd subA1(ptr_a1,a1_length*A1.NABRows_,A1.NABCols_);

// #pragma omp single
// 				async_2.get(A2.AReadPos_[A2.read_restrict_*reada1].first,
// 						A2.AReadPos_[A2.read_restrict_*reada1].second,
// 						*(A2.Array_), ptr_a23_2);

// 				for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){
// 					bool loop_restriction = 0;
// 					if(A2.read_restrict_ && reada2 == reada1)loop_restriction = 1;
// 					//				for(size_t reada2 = 0; reada2 < A2.NAReads_ ; reada2++ ){
// 					//					bool loop_restriction = 0;

// 					//#pragma omp single
// 					//					A2.Array_->get(ptr_a23, A2.AReadPos_[reada2].first, A2.AReadPos_[reada2].second);
// #pragma omp single
// 					{
// 						async_2.wait();
// 						std::swap(ptr_a23,ptr_a23_2);
// 						if(reada2 < A2.NAReads_-1 && A2.NAReads_ > 1)
// 							async_2.get(A2.AReadPos_[reada2+1].first,
// 									A2.AReadPos_[reada2+1].second,
// 									*(A2.Array_), ptr_a23_2);
// 					}//omp single

// 					size_t a2_start = A2.ABlocks_[reada2].first;
// 					size_t a2_stop = A2.ABlocks_[reada2].second;
// 					size_t a2_length = a2_stop - a2_start;
// 					MapMatrixXd subA2(ptr_a23,a2_length*A2.NABRows_,A2.NABCols_);

// 					FunctorOne(eri_vooo,
// 							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 							subA2,a2_start,a2_stop, A2.NABRows_,A2.NABCols_,
// 							loop_restriction);
// 				}//reada2


// #pragma omp single
// 				async_3.get(A3.AReadPos_[0].first,
// 						A3.AReadPos_[0].second,
// 						*(A3.Array_), ptr_a23_2);

// 				//not implemented for functor				for(size_t reada3 = A3.restrict_*reada1; reada3 < A3.NAReads_ ; reada3++ ){
// 				//not implemented for functor				bool loop_restriction = 0;
// 				//not implemented for functor					if(A3.restrict_ && reada3 == reada1)loop_restriction = 1;
// 				for(size_t reada3 = 0; reada3 < A3.NAReads_ ; reada3++ ){
// 					bool loop_restriction = 0;

// 					//#pragma omp single
// 					//					A3.Array_->get(ptr_a23, A3.AReadPos_[reada3].first, A3.AReadPos_[reada3].second);
// #pragma omp single
// 					{
// 						async_3.wait();
// 						std::swap(ptr_a23,ptr_a23_2);
// 						if(reada3 < A3.NAReads_-1 && A3.NAReads_ > 1)
// 							async_2.get(A3.AReadPos_[reada3+1].first,
// 									A3.AReadPos_[reada3+1].second,
// 									*(A3.Array_), ptr_a23_2);
// 					}//omp single

// 					size_t a3_start = A3.ABlocks_[reada3].first;
// 					size_t a3_stop = A3.ABlocks_[reada3].second;
// 					size_t a3_length = a3_stop - a3_start;
// 					MapMatrixXd subA3(ptr_a23,a3_length*A3.NABRows_,A3.NABCols_);

// 					//#pragma omp barrier
// 					FunctorTwo(eri_vvvo,
// 							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 							subA3,a3_start,a3_stop, A3.NABRows_,A3.NABCols_,
// 							loop_restriction);
// 				}//reada3

// 			}//reada1



// 			delete ptr_eri;

// 		}//omp parallel

// 		delete [] ptr_a23;
// 		delete [] ptr_a23_2;
// 		delete [] ptr_a1;
// 		delete [] ptr_a1_2;

// 	}//DoOperatrion2Async







// 	//	template <typename FUNCTOR >
// 	//	void DoOperation2(FUNCTOR Functor){
// 	//
// 	//		ArrayData &A1 = Array_[0];
// 	//		ArrayData &A2 = Array_[1];
// 	//
// 	//		void * a1 = NULL;
// 	//		posix_memalign(&a1, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 	//		double *ptr_a1 = new(a1) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];
// 	//
// 	//		void * a2 = NULL;
// 	//		posix_memalign(&a2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 	//		double *ptr_a2 = new(a2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];
// 	//
// 	//#pragma omp parallel
// 	//		{
// 	//			void * aeri = NULL;
// 	//			posix_memalign(&aeri, 16, A1.NABRows_*A2.NABRows_*sizeof(double) );
// 	//			double *ptr_eri = new(aeri) double[A1.NABRows_*A2.NABRows_];
// 	//
// 	//			MapMatrixXd eri(ptr_eri,A1.NABRows_,A2.NABRows_);
// 	//
// 	//			for(size_t reada1 = 0; reada1 < A1.NAReads_ ; reada1++ ){
// 	//#pragma omp single
// 	//				A1.Array_->get(ptr_a1, A1.AReadPos_[reada1].first, A1.AReadPos_[reada1].second);
// 	//				size_t a1_start = A1.ABlocks_[reada1].first;
// 	//				size_t a1_stop = A1.ABlocks_[reada1].second;
// 	//				size_t a1_length = a1_stop - a1_start;
// 	//				MapMatrixXd subA1(ptr_a1,a1_length*A1.NABRows_,A1.NABCols_);
// 	//
// 	//				for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){
// 	//					bool restriction = 0;
// 	//					if(A2.read_restrict_ && reada2 == reada1) restriction = 1;
// 	//					//					for(size_t reada2 = 0; reada2 < A2.NAReads_ ; reada2++ ){
// 	//					//						bool restriction = 0;
// 	//
// 	//
// 	//#pragma omp single
// 	//					A2.Array_->get(ptr_a2, A2.AReadPos_[reada2].first, A2.AReadPos_[reada2].second);
// 	//					size_t a2_start = A2.ABlocks_[reada2].first;
// 	//					size_t a2_stop = A2.ABlocks_[reada2].second;
// 	//					size_t a2_length = a2_stop - a2_start;
// 	//					MapMatrixXd subA2(ptr_a2,a2_length*A2.NABRows_,A2.NABCols_);
// 	//
// 	//					Functor(eri,
// 	//							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 	//							subA2,a2_start,a2_stop, A2.NABRows_,A2.NABCols_,
// 	//							restriction);
// 	//
// 	//				}//reada2
// 	//
// 	//			}//reada1
// 	//
// 	//
// 	//
// 	//			delete ptr_eri;
// 	//
// 	//		}//omp parallel
// 	//
// 	//		delete ptr_a2;
// 	//		delete ptr_a1;
// 	//
// 	//	}//DoOperatrion2






// 	template <typename FUNCTOR >
// 	void DoOperation2Async(FUNCTOR Functor){

// 		ArrayData &A1 = Array_[0];
// 		ArrayData &A2 = Array_[1];


// 		void * a1 = NULL;
// 		posix_memalign(&a1, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_a1 = new(a1) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		void * a1_2 = NULL;
// 		posix_memalign(&a1_2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_a1_2 = new(a1_2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		void * a2 = NULL;
// 		posix_memalign(&a2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		double *ptr_a2 = new(a2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		double *ptr_a2_2 = NULL;

// 		::cchem::ri::async async_1;
// 		::cchem::ri::async async_2;

// 		if(A2.NAReads_ > 1){

// 			void * a2_2 = NULL;
// 			posix_memalign(&a2_2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 			ptr_a2_2 = new(a2_2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		}else{

// 			async_2.get(A2.AReadPos_[0].first,
// 					A2.AReadPos_[0].second,
// 					*(A2.Array_), ptr_a2);
// 			async_2.wait();
// 		}


// 		async_1.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_a1_2);

// #pragma omp parallel
// 		{
// 			void * aeri = NULL;
// 			posix_memalign(&aeri, 16, A1.NABRows_*A2.NABRows_*sizeof(double) );
// 			double *ptr_eri = new(aeri) double[A1.NABRows_*A2.NABRows_];

// 			MapMatrixXd eri(ptr_eri,A1.NABRows_,A2.NABRows_);

// 			for(size_t reada1 = 0; reada1 < A1.NAReads_ ; reada1++ ){
// 				//#pragma omp single
// 				//				A1.Array_->get(ptr_a1, A1.AReadPos_[reada1].first, A1.AReadPos_[reada1].second);

// #pragma omp single
// 				{
// 					async_1.wait();
// 					std::swap(ptr_a1,ptr_a1_2);
// 					if(reada1 < A1.NAReads_-1 && A1.NAReads_ > 1)
// 						async_1.get(A1.AReadPos_[reada1+1].first,
// 								A1.AReadPos_[reada1+1].second,
// 								*(A1.Array_), ptr_a1_2);
// 				}//omp single

// 				size_t a1_start = A1.ABlocks_[reada1].first;
// 				size_t a1_stop = A1.ABlocks_[reada1].second;
// 				size_t a1_length = a1_stop - a1_start;
// 				MapMatrixXd subA1(ptr_a1,a1_length*A1.NABRows_,A1.NABCols_);

// #pragma omp single
// 				{
// 					if(A2.NAReads_ > 1){
// 						async_2.get(A2.AReadPos_[A2.read_restrict_*reada1].first,
// 								A2.AReadPos_[A2.read_restrict_*reada1].second,
// 								*(A2.Array_), ptr_a2_2);
// 					}//if(A2.NAReads_ > 1)
// 				}//omp single

// 				for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){
// 					bool restriction = 0;
// 					if(A2.read_restrict_ && reada2 == reada1) restriction = 1;
// 					//					for(size_t reada2 = 0; reada2 < A2.NAReads_ ; reada2++ ){
// 					//						bool restriction = 0;


// 					//#pragma omp single
// 					//					A2.Array_->get(ptr_a2, A2.AReadPos_[reada2].first, A2.AReadPos_[reada2].second);

// #pragma omp single
// 					{
// 						if(A2.NAReads_ > 1){
// 							async_2.wait();
// 							std::swap(ptr_a2,ptr_a2_2);
// 							if(reada2 < A2.NAReads_-1 && A2.NAReads_ > 1)
// 								async_2.get(A2.AReadPos_[reada2+1].first,
// 										A2.AReadPos_[reada2+1].second,
// 										*(A2.Array_), ptr_a2_2);
// 						}//if(A2.NAReads_ > 1)
// 					}//omp single


// 					size_t a2_start = A2.ABlocks_[reada2].first;
// 					size_t a2_stop = A2.ABlocks_[reada2].second;
// 					size_t a2_length = a2_stop - a2_start;
// 					MapMatrixXd subA2(ptr_a2,a2_length*A2.NABRows_,A2.NABCols_);

// 					Functor(eri,
// 							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 							subA2,a2_start,a2_stop, A2.NABRows_,A2.NABCols_,
// 							restriction);

// 				}//reada2

// 			}//reada1



// 			delete ptr_eri;

// 		}//omp parallel

// 		delete [] ptr_a2;
// 		delete [] ptr_a1;
// 		delete [] ptr_a2_2;
// 		delete [] ptr_a1_2;

// 	}//DoOperatrion2Async
























// 	template <typename FUNCTOR >
// 	void DoOpAsync_RCRC(FUNCTOR Functor){ //luke

// 		ArrayData &A1 = Array_[0];
// 		ArrayData &A2 = Array_[1];


// 		void * a1 = NULL;
// 		posix_memalign(&a1, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_a1 = new(a1) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		void * a1_2 = NULL;
// 		posix_memalign(&a1_2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_a1_2 = new(a1_2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		void * a2 = NULL;
// 		posix_memalign(&a2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		double *ptr_a2 = new(a2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		double *ptr_a2_2 = NULL;

// 		::cchem::ri::async async_1;
// 		::cchem::ri::async async_2;

// 		if(A2.NAReads_ > 1){

// 			void * a2_2 = NULL;
// 			posix_memalign(&a2_2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 			ptr_a2_2 = new(a2_2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		}else{

// 			async_2.get(A2.AReadPos_[0].first,
// 					A2.AReadPos_[0].second,
// 					*(A2.Array_), ptr_a2);
// 			async_2.wait();
// 		}


// 		async_1.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_a1_2);

// 		//#pragma omp parallel
// 		{
// 			//			void * aeri = NULL;
// 			//			posix_memalign(&aeri, 16, A1.NABRows_*A2.NABRows_*sizeof(double) );
// 			//			double *ptr_eri = new(aeri) double[A1.NABRows_*A2.NABRows_];
// 			//
// 			//			MapMatrixXd eri(ptr_eri,A1.NABRows_,A2.NABRows_);

// 			for(size_t reada1 = 0; reada1 < A1.NAReads_ ; reada1++ ){
// 				//				std::cout << "reada1:" << reada1 << std::endl;
// 				//#pragma omp single
// 				//				A1.Array_->get(ptr_a1, A1.AReadPos_[reada1].first, A1.AReadPos_[reada1].second);

// #pragma omp single
// 				{
// 					async_1.wait();
// 					std::swap(ptr_a1,ptr_a1_2);
// 					if(reada1 < A1.NAReads_-1 && A1.NAReads_ > 1)
// 						async_1.get(A1.AReadPos_[reada1+1].first,
// 								A1.AReadPos_[reada1+1].second,
// 								*(A1.Array_), ptr_a1_2);
// 				}//omp single

// 				size_t a1_start = A1.ABlocks_[reada1].first;
// 				size_t a1_stop = A1.ABlocks_[reada1].second;
// 				size_t a1_length = a1_stop - a1_start;
// 				//				MapMatrixXd subA1(ptr_a1,a1_length*A1.NABRows_,A1.NABCols_);

// #pragma omp single
// 				{
// 					if(A2.NAReads_ > 1){
// 						async_2.get(A2.AReadPos_[A2.read_restrict_*reada1].first,
// 								A2.AReadPos_[A2.read_restrict_*reada1].second,
// 								*(A2.Array_), ptr_a2_2);
// 					}//if(A2.NAReads_ > 1)
// 				}//omp single

// 				for(size_t reada2 = A2.read_restrict_*reada1; reada2 < A2.NAReads_ ; reada2++ ){
// 					//					std::cout << "reada2:" << reada2 << std::endl;
// 					bool restriction = 0;
// 					if(A2.read_restrict_ && reada2 == reada1) restriction = 1;
// 					//					for(size_t reada2 = 0; reada2 < A2.NAReads_ ; reada2++ ){
// 					//						bool restriction = 0;


// 					//#pragma omp single
// 					//					A2.Array_->get(ptr_a2, A2.AReadPos_[reada2].first, A2.AReadPos_[reada2].second);

// #pragma omp single
// 					{
// 						if(A2.NAReads_ > 1){
// 							async_2.wait();
// 							std::swap(ptr_a2,ptr_a2_2);
// 							if(reada2 < A2.NAReads_-1 && A2.NAReads_ > 1)
// 								async_2.get(A2.AReadPos_[reada2+1].first,
// 										A2.AReadPos_[reada2+1].second,
// 										*(A2.Array_), ptr_a2_2);
// 						}//if(A2.NAReads_ > 1)
// 					}//omp single


// 					size_t a2_start = A2.ABlocks_[reada2].first;
// 					size_t a2_stop = A2.ABlocks_[reada2].second;
// 					size_t a2_length = a2_stop - a2_start;
// 					//					MapMatrixXd subA2(ptr_a2,a2_length*A2.NABRows_,A2.NABCols_);

// 					//					Functor(eri,
// 					//							subA1, a1_start,a1_stop, A1.NABRows_,A1.NABCols_,
// 					//							subA2,a2_start,a2_stop, A2.NABRows_,A2.NABCols_,
// 					//							restriction);

// 					Functor(ptr_a1, a1_start,a1_stop,
// 							ptr_a2,a2_start,a2_stop);

// 				}//reada2

// 			}//reada1



// 			//			delete ptr_eri;

// 		}//omp parallel

// 		delete [] ptr_a2;
// 		delete [] ptr_a1;
// 		delete [] ptr_a2_2;
// 		delete [] ptr_a1_2;

// 	}//DoOpAsync_RCRC























// 	template< typename OP >
// 	void DoOpAsync_RC(OP Functor) {

// 		ArrayData &A1 = Array_[0];
// 		std::cout << A1.NABRows_ << " " << A1.NABCols_ << std::endl;
// 		//storage for '2/3' transformed integrals
// 		void * aints = NULL;
// 		posix_memalign(&aints, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_ints = new(aints) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		//storage for '2/3' transformed integrals
// 		void * aints2 = NULL;
// 		posix_memalign(&aints2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_ints2 = new(aints2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		std::cout << A1.AReadPos_[0].first.at(0) << " " << A1.AReadPos_[0].first.at(1) << std::endl;
// 		std::cout << A1.AReadPos_[0].second.at(0) << " " << A1.AReadPos_[0].second.at(1) << std::endl;
// 		::cchem::ri::async async;
// 		async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

// #pragma omp parallel
// 		{

// 			for(size_t read = 0; read < A1.NAReads_; read++ ){


// #pragma omp single
// 				{
// 					async.wait();
// 					std::swap(ptr_ints,ptr_ints2);

// 					if(A1.NAReads_ > 1 && read < A1.NAReads_-1){
// 						async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
// 								*(A1.Array_), ptr_ints2);
// 					}//(A1.NAReads_ > 0)
// 				}//omp single

// #pragma omp barrier

// 				size_t b_start = A1.ABlocks_[read].first;
// 				size_t b_stop = A1.ABlocks_[read].second;
// 				size_t b_length = b_stop - b_start;



// 				Functor(ptr_ints, A1.NABRows_, b_start, b_stop, A1.NABCols_);


// 			}//read

// 		}//omp parallel



// 		delete [] ptr_ints2;
// 		delete [] ptr_ints;
// 		//std::cout << "leaving " << std::flush << std::endl;

// 	}//DoOperationAsync_RC







// 	template< typename OP >
// 	void DoOpAsync_RCWC(OP Functor) {
// 		//read column, write column
// 		ArrayData &A1 = Array_[0];
// 		ArrayData &A2 = Array_[1];

// 		//array one buffers
// 		void * aints = NULL;
// 		posix_memalign(&aints, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_ints = new(aints) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		void * aints2 = NULL;
// 		posix_memalign(&aints2, 16, A1.MaxNAB_*A1.NABRows_*A1.NABCols_*sizeof(double) );
// 		double *ptr_ints2 = new(aints2) double[A1.MaxNAB_*A1.NABRows_*A1.NABCols_];

// 		//array two buffers
// 		void * aCoeff = NULL;
// 		posix_memalign(&aCoeff, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		double *ptr_Coeff = new(aCoeff) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];

// 		void * aCoeff2 = NULL;
// 		posix_memalign(&aCoeff2, 16, A2.MaxNAB_*A2.NABRows_*A2.NABCols_*sizeof(double) );
// 		double *ptr_Coeff2 = new(aCoeff2) double[A2.MaxNAB_*A2.NABRows_*A2.NABCols_];


// 		::cchem::ri::async async;
// 		async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

// #pragma omp parallel
// 		{
// 			//			std::cout << "check " << omp_get_thread_num() << std::endl;
// #pragma omp barrier
// 			//#pragma omp single
// 			//			async.get(A1.AReadPos_[0].first, A1.AReadPos_[0].second, *(A1.Array_), ptr_ints2);

// 			//			std::cout << "start reads for A1" << std::endl;
// 			for(size_t read = 0; read < A1.NAReads_; read++ ){
// 				//				std::cout << "done with read " << read << " of " << A1.NAReads_ << std::endl;
// 				//				#pragma omp single
// 				//								A1.Array_->get(ptr_ints, A1.AReadPos_[read].first, A1.AReadPos_[read].second);

// #pragma omp single
// 				{
// 					async.wait();
// 					std::swap(ptr_ints,ptr_ints2);

// 					if(A1.NAReads_ > 1){

// 						if(read == 0){
// 							async.get(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
// 									*(A1.Array_), ptr_ints2);
// 						}else{
// 							std::swap(ptr_Coeff,ptr_Coeff2);
// 							if(read < A1.NAReads_-1){
// 								async.get_put(A1.AReadPos_[read+1].first, A1.AReadPos_[read+1].second,
// 										A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
// 										*(A1.Array_), ptr_ints2,
// 										*(A2.Array_), ptr_Coeff2);
// 							}else{
// 								async.put(A2.AReadPos_[read-1].first, A2.AReadPos_[read-1].second,
// 										*(A2.Array_), ptr_Coeff2);
// 							}//(read < A1.NAReads_-1)
// 						}//(read == 0)
// 					}//(A1.NAReads_ > 0)
// 				}//omp master

// #pragma omp barrier

// 				size_t b_start = A1.ABlocks_[read].first;
// 				size_t b_stop = A1.ABlocks_[read].second;
// 				size_t b_length = b_stop - b_start;

// 				//				MapMatrixXd Integrals(ptr_ints,b_length*A1.NABRows_,A1.NABCols_);   // '2/3' transformed integrals (ij|Q)
// 				//				MapMatrixXd Coeff(ptr_Coeff,b_length*A2.NABRows_,A2.NABCols_);   // '3/3' transformed integrals (ij|Q)

// 				Functor(ptr_ints, b_start, b_stop,
// 						ptr_Coeff);


// 			}//read

// 		}//omp parallel

// 		async.wait();
// 		A2.Array_->put(ptr_Coeff, A2.AReadPos_[A2.NAReads_-1].first,
// 				A2.AReadPos_[A2.NAReads_-1].second);

// 		delete [] ptr_Coeff2;
// 		delete [] ptr_Coeff;
// 		delete [] ptr_ints2;
// 		delete [] ptr_ints;
// 		//std::cout << "leaving " << std::flush << std::endl;

// 	}//DoOpAsync_RCWC



// };//ri-async




// }//namespace detail
// }//namespace rimp2_gradient
// }//namespace cchem


#endif /* LIBCCHEM_SRC_MP2_RI_ASYNC_HPP_ */
