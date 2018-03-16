

#ifndef SRC_RI_MP2_DEVICE_HPP_
#define SRC_RI_MP2_DEVICE_HPP_

//if defined(HAVE_CUBLAS) && defined(CCHEM_INTEGRALS_ERI_CUDA)

//#define CCHEM_MP2_CUDA
//#warning GPU MP2 disabled due to being slow
// #undef CCHEM_MP2_CUDA
//test#endif
//#include <cuda_profiler_api.h>
//#include "/share/apps/cuda7/include/cuda_profiler_api.h"
#if HAVE_CUBLAS
#include "cublas.hpp"
#include <cuda_runtime.h>
#include <cuda.h>
#endif


namespace cchem {

    namespace rimp2 {
    class device;
    }
} //namespace cchem





#if HAVE_CUBLAS

namespace cchem {

/*
  @brief wrapper function to check for GPU errors
  @author LBR
  @detail wrapper function to check for GPU errors
  @param code cuda error code
  @param file filename where function was called from
  @param line where function was called from
  @param Abort abort switch
*/
    inline void GPUassert(cudaError_t code, char const * file, int line, bool Abort=true)
    {
	if (code != 0) {
	    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),file,line);
	    if (Abort) exit(code);
	}       
    }//GPUassert
    
#define GPUerrchk(ans) { GPUassert((ans), __FILE__, __LINE__); }
    
    namespace rimp2 {




/*
  @brief device class
  @author LBR
  @detail this sets up a device object
  @param nDevices number of devices (on node)
  @param nStreams number of streams per device
  @param cudaBlocks number of cuda blocks per kernel launch
  @param cudaThreadsPB number of threads per cuda block
*/
    	class device : boost::noncopyable {

    	    const int nDevices_;
    	    const int nStreams_;
	    const int cudaBlocks_;
	    const int cudaThreadsPB_;
    	    std::vector< cudaStream_t * > pVecStreams_;
    	    std::vector< cublasHandle_t > vecHandles_;
	    int cudaRuntimeVersion_;
    	    size_t maxDeviceMemory_;
	    bool doPrint_;
    	public:


	    static int getNumberOfDevices(){
		int ndevice;
		GPUerrchk(cudaGetDeviceCount(&ndevice ));
		return ndevice;
	    };


	    static int getNumberOfStreams(){
		char *ev = std::getenv("CCHEM_NUMBER_CUDA_STREAMS");
		int nstreams = (ev ? atoi(ev) : 1);
		return nstreams;
	    };


    	    device (const int nDevices, const int nStreams,
		    const int cudaBlocks, const int cudaThreadsPB,
		    bool doPrint): 
    		nDevices_(nDevices), nStreams_(nStreams),
		cudaBlocks_(cudaBlocks), cudaThreadsPB_(cudaThreadsPB),
    		pVecStreams_(std::vector< cudaStream_t * >(nDevices) ),
    		vecHandles_( std::vector< cublasHandle_t >(nDevices) ),
		cudaRuntimeVersion_(0),
    		maxDeviceMemory_(0),
		doPrint_(doPrint)
    	    { 
    		maxDeviceMemory_ = 0;

		int cudaRuntimeVersion;
		cudaRuntimeGetVersion(&cudaRuntimeVersion);
		cudaRuntimeVersion_ = cudaRuntimeVersion;

		if(doPrint)std::cout << "found " 
				     << nDevices_ << " devices: " 
				     << "each device will setup " 
				     << nStreams_ << " cuda streams" << std::endl;
		
    		for(int idevice = 0; idevice < nDevices_; idevice++){



    		    //switch device
    		    GPUerrchk(cudaSetDevice( idevice ));

		    //reset device for good measure
		    cudaDeviceReset();


    		    cudaDeviceProp prop;
    		    GPUerrchk(cudaGetDeviceProperties(&prop, idevice));
    		    //std::printf("             Device Number: %d\n", idevice);
    		    //std::printf("               Device name: %s\n", prop.name);
		    //std::printf("Unified Virtual Addressing: %i\n",prop.unifiedAddressing);
		    
                    if(doPrint)std::cout << "             Device Number: "
					 << idevice << std::endl;
                    if(doPrint)std::cout << "               Device Name: "
					 << prop.name << std::endl;
                    if(doPrint)std::cout << "      CUDA runtime verison: "
					 << cudaRuntimeVersion_/1000 << "." 
					 << (cudaRuntimeVersion_%100)/10 << std::endl;			


     //     std::printf("  Memory Clock Rate (KHz): %d\n",
     //		 prop.memoryClockRate);
     //     std::printf("  Memory Bus Width (bits): %d\n",
     //		 prop.memoryBusWidth);
     //     std::printf("  Peak Memory Bandwidth (GB/s): %f\n",
     //		 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
     //     std::printf("  deviceOverlap? %d\n",prop.deviceOverlap);
     //     std::printf("  multiprocessor count %d\n",prop.multiProcessorCount);
     //     std::printf("  PCI bus ID  %d\n",prop.pciBusID);
     //     std::printf("  PCI device ID %d\n",prop.pciDeviceID);

    		    //allocate memory for streams
		    pVecStreams_[idevice] = new cudaStream_t[nStreams_*sizeof(cudaStream_t)];

    		    //create nStreams_ streams on current device
    		    for(int istream = 0; istream < nStreams_; istream++){
    			GPUerrchk(cudaStreamCreate(&pVecStreams_[idevice][istream]));
    			//GPUerrchk(cudaStreamCreateWithFlags(&pVecStreams_[idevice][istream],
			//cudaStreamNonBlocking));
    		    }//istream

    		    //create cublas handles on current device
    		    if (cublasCreate(&vecHandles_[idevice]) != CUBLAS_STATUS_SUCCESS)
    			{ fprintf(stdout, "CUBLAS: cublasCreate failed!\n");
    			    cudaDeviceReset();
    			    exit(EXIT_FAILURE); }

    		    //get the amount of free memory on the current device
    		    size_t memory;
    		    size_t total;
    		    GPUerrchk( cudaMemGetInfo(&memory,  &total ));
    		    if(idevice == 0)maxDeviceMemory_ = memory;

                    if(doPrint)std::cout << "       Free Memory (bytes): "
					 << memory << std::endl;


		    
    		    //if more than one device is present, the max working memory is
    		    // set to the device with the smallest amount of memory available 
    		    // (e.g. k20 and k80, max memory will be that of the k20)
    		    maxDeviceMemory_ = std::min(maxDeviceMemory_,memory);
		    
    		    if(doPrint)std::cout << std::endl;
    		}//idevice
		
    	    }//device()
	    
	    
    	    ~device(){
    		for(int idevice = 0; idevice < nDevices_; idevice++){

		    for(int istream = 0; istream < nStreams_; istream++){
			GPUerrchk( cudaStreamDestroy(pVecStreams_[idevice][istream]) );
		    }//istream

		    delete [] pVecStreams_[idevice];

    		    if (cublasDestroy(vecHandles_[idevice]) != CUBLAS_STATUS_SUCCESS)
    			{ fprintf(stdout, "CUBLAS: cublasDestroy failed!\n");
    			    cudaDeviceReset();
    			    exit(EXIT_FAILURE); }
		    
    		}//idevice
    	    }

    	    //accessor to cuda streams
    	    std::vector< cudaStream_t* >& cudaStreams(){return pVecStreams_;}
	    
    	    //accessor to cublas handles
    	    std::vector< cublasHandle_t >& cublasHandles(){return vecHandles_;}

    	    //accessor to working device memory (bytes)
    	    size_t maxBytes(){ return maxDeviceMemory_; }

    	    //accessor to working device memory (bytes)
    	    size_t maxWords(){ return maxDeviceMemory_/sizeof(double); }

    	    void synchronize(){
    		for (int idevice = 0; idevice < nDevices_; idevice++){
    		    GPUerrchk( cudaSetDevice( idevice ) );
    		    GPUerrchk( cudaDeviceSynchronize() );
    		}//idevice
    	    }//synchronize()

	    size_t availableBytes(){

		size_t freeMemory;
    		for(int idevice = 0; idevice < nDevices_; idevice++){

    		    size_t memory;
    		    size_t total;

    		    GPUerrchk( cudaSetDevice( idevice ) );
    		    GPUerrchk( cudaMemGetInfo(&memory,  &total ));
		    if(idevice == 0)freeMemory = memory;
    		    std::cout << "    Free Memory : " << memory << std::endl;
		    if(memory < freeMemory)freeMemory = memory;

		}//idevice
		return freeMemory;
	    }//

	    int getNumDevices(){
		return nDevices_;
	    }

	    int getNumStreams(){
		return nStreams_;
	    }

	    int getCudaBlocks(){
		return cudaBlocks_;
	    }

	    int getCudaThreadsPB(){
		return cudaThreadsPB_;
	    }

	    int getCudaRuntimeVersion(){
		return cudaRuntimeVersion_;
	    }


	    
    	};//class device
	
	
    } //namespace rimp2

    } //namespace cchem

#endif

#endif //SRC_RI_MP2_DEVICE_HPP_
