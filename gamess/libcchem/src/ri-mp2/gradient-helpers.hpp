/*
 * gradient-helpers.hpp
 *
 *  Created on: Nov 4, 2015
 *      Author: luke
 */

#ifndef SRC_RI_MP2_GRADIENT_HELPERS_HPP_
#define SRC_RI_MP2_GRADIENT_HELPERS_HPP_

#include "thread.hpp"
#include <Eigen/Dense>
//#include <boost/numeric/ublas/vector.hpp>
#include "core/wavefunction.hpp"



#if HAVE_CUBLAS
#include <device.hpp>
#endif

#include <math.hpp>
namespace cchem{
namespace rimp2_gradient{
namespace detail {

typedef Eigen::Map<Eigen::MatrixXd,Eigen::AutoAlign> MapMatrixXd;
typedef const Eigen::Map<const Eigen::MatrixXd,Eigen::AutoAlign> ConstMapMatrixXd;

typedef boost::numeric::ublas::matrix<
		double, boost::numeric::ublas::column_major> Matrix;






class gradientTerms
{

	std::vector<Eigen::MatrixXd> gradients_;
	unsigned int centers_;
	int nTerms_;
public:

	enum term_ {DFIntGradient, densityForce, twoCenterNonSeparable, threeCenterNonSeparable, nuclearForce,
		T, V };

	gradientTerms(unsigned int centers):
		centers_(centers),nTerms_(7)
	{

		for(int i = 0; i < nTerms_; i++){
			gradients_.push_back( Eigen::MatrixXd(centers_,3) );
			gradients_[i].setZero();
		}//i

	};//gradientTerms

	~gradientTerms(){};

	void gradientSum(MapMatrixXd &egGlobal){
		egGlobal.setZero();
		for (int i = 0; i < nTerms_; i++){
			egGlobal += gradients_[i];
		}//i
	};//gradientSum


	double *getPointer(int tag){ return gradients_[tag].data(); };

	Eigen::MatrixXd &getMatrix(term_ tag){ return gradients_[tag]; };

	void printGradientTerms(){

		std::cout << "4-center JK derivative contribution " << std::endl;
		std::cout << this->gradients_[0] << std::endl;

		std::cout << std::endl << "gradient of density forces" << std::endl;
		std::cout << this->gradients_[1] << std::endl;

		std::cout << std::endl << "gradient of two-center forces" << std::endl;
		std::cout << this->gradients_[2] << std::endl;

		std::cout << std::endl << "gradient of three-center forces" << std::endl;
		std::cout << this->gradients_[3] << std::endl;

		std::cout << std::endl << "gradient of nuclear derivative" << std::endl;
		std::cout << this->gradients_[4] << std::endl;

		std::cout << std::endl << "gradient of kinetic energy" << std::endl;
		std::cout << this->gradients_[5] << std::endl;

		std::cout << std::endl << "gradient of electron-nuclear (w/Hellman) energy" << std::endl;
		std::cout << this->gradients_[6] << std::endl;

	}//printGradientTerms

};



//for energy weighted density
struct wmat_vovo_functor{
private:
	size_t no_,nv_;
	MapMatrixXd &pmat_;
	MapMatrixXd &wmat_;
public:
	wmat_vovo_functor(size_t no, size_t nv, MapMatrixXd &pmat, MapMatrixXd &wmat):
		no_(no),nv_(nv), pmat_(pmat), wmat_(wmat) {}
	void operator() (MapMatrixXd &eri,
			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
			MapMatrixXd &bba, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
			bool &loop_restriction);
}; //vvvo_functor

struct wmat_vooo_functor{
private:
	size_t no_,nv_;
	MapMatrixXd &pmat_;
	MapMatrixXd &wmat_;
public:
	wmat_vooo_functor(size_t no, size_t nv, MapMatrixXd &pmat, MapMatrixXd &wmat):
		no_(no),nv_(nv), pmat_(pmat), wmat_(wmat) {}
	void operator() (MapMatrixXd &eri,
			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
			MapMatrixXd &bba, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
			bool &loop_restriction);
}; //vvvo_functor

struct wmat_vvoo_functor{
private:
	size_t no_,nv_;
	MapMatrixXd &pmat_;
	MapMatrixXd &wmat_;
public:
	wmat_vvoo_functor(size_t no, size_t nv, MapMatrixXd &pmat, MapMatrixXd &wmat):
		no_(no),nv_(nv), pmat_(pmat), wmat_(wmat) {}
	void operator() (MapMatrixXd &eri,
			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
			MapMatrixXd &bba, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
			bool &loop_restriction);
}; //vvoo_functor

struct wmat_oooo_functor{
private:
	size_t no_,nv_;
	MapMatrixXd &pmat_;
	MapMatrixXd &wmat_;
public:
	wmat_oooo_functor(size_t no, size_t nv, MapMatrixXd &pmat, MapMatrixXd &wmat):
		no_(no),nv_(nv), pmat_(pmat), wmat_(wmat) {}
	void operator() (MapMatrixXd &eri,
			MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
			MapMatrixXd &bba, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
			bool &loop_restriction);
}; //oooo_functor


//BARE_IA -> CIAQ (coefficient) or IA (Cholesky operation)
struct transform_functor{
private:
	MapMatrixXd &matrix_;
public:
	transform_functor(MapMatrixXd &matrix):
		matrix_(matrix) {}
	void operator() (MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
			MapMatrixXd &bba, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols);
	//	void operator() (double *ptr1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &A1_RBSize,
	//								   size_t &A1_BCStart, size_t &A1_BCStop, size_t &A1_CBSize,
	//                   double *ptr2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &A2_RBSize,
	//								   size_t &A2_BCStart, size_t &A2_BCStop, size_t &A2_CBSize)
}; //transform_functor


//BARE_IA -> CIAQ (coefficient) or IA (Cholesky operation)
struct transform_functor_new{
private:
//	MapMatrixXd &matrix_;
	double *matrix_;
public:
//	transform_functor_new(MapMatrixXd &matrix):
		transform_functor_new(double *matrix):
		matrix_(matrix) {}
//	void operator() (MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
//			MapMatrixXd &bba, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols);
		void operator() (double *ptr1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &A1_RBSize,
									   size_t &A1_BCStart, size_t &A1_BCStop, size_t &A1_CBSize,
	                     double *ptr2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &A2_RBSize,
									   size_t &A2_BCStart, size_t &A2_BCStop, size_t &A2_CBSize);
}; //transform_functor_new



    //BARE_IA -> CIAQ (coefficient) or IA (Cholesky operation)
    struct transform_functor_new2 : boost::noncopyable{
    private:
	double *matrix_;
	utility::timer timerTest;
	utility::timer::value_type CPUTime;
    public:
	transform_functor_new2(double *matrix):
	    matrix_(matrix){}

	~transform_functor_new2()
	{
	    //std::cout << std::endl << "VVVVV actual CPU time for transformation (below): " 
	    //<< CPUTime << std::endl;

	}//~
	void operator() (double *ptr1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &A1_RBSize,
			 size_t &A1_BCStart, size_t &A1_BCStop, size_t &A1_CBSize,
			 double *ptr2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &A2_RBSize,
			 size_t &A2_BCStart, size_t &A2_BCStop, size_t &A2_CBSize);
    }; //transform_functor_new2






// #if HAVE_CUBLAS

// //BARE_IA -> CIAQ (coefficient) or IA (Cholesky operation) -- GPU
//     struct transformFunctorGPU{ // : boost::noncopyable {
//     private:
// 	double *metric_;
// 	const size_t gpuFreeB_;
// 	const size_t rows_, nl_;
// 	const int nDevices_, nStreams_;
// 	std::vector<double *> pDeviceVectorMetric_;
// 	std::vector<double*> pDeviceVectorBare_;
// 	std::vector<double*> pDeviceVectorDressed_;
// 	std::vector<cublasHandle_t> &dVHandles_;

// 	//	std::vector< std::vector<double*> > pDeviceVectorBare_;
// 	//	std::vector< std::vector<double*> > pDeviceVectorDressed_;

//     public:

// 	transformFunctorGPU(double *metric, const size_t gpuFreeB, const size_t rows,
// 			    const size_t nl, const int nDevices, const int nStreams,
// 			    std::vector<cublasHandle_t> &dVHandles):
// 	    metric_(metric), gpuFreeB_(gpuFreeB),rows_(rows), nl_(nl), 
// 	    nDevices_(nDevices), nStreams_(nStreams),
// 	    pDeviceVectorMetric_( std::vector<double *> (nDevices) ),
// 	    pDeviceVectorBare_( std::vector<double *> (nDevices) ),
// 	    pDeviceVectorDressed_( std::vector<double *> (nDevices) ),
// 	    dVHandles_ ( dVHandles )



// 	{
// 	    std::cout << " initializing trasnformFunctorGPU" << std::flush << std::endl;
// 	    const int maxRowsPerDevice = (rows_/nDevices_) + ( rows_%nDevices_); 
// 	    //	    std::cout << "maxRowsPerDevice: "<< maxRowsPerDevice << std::endl;
// 	    const int maxMatrixSizePerDevice = maxRowsPerDevice*nl_;

// 	    //	    std::cout << "maxMatrixSizePerDevice: "<< maxMatrixSizePerDevice << std::endl;
// 	    const size_t matrixSizeMetric = nl_*nl_;


// 	    for(int iDevice = 0; iDevice < nDevices_; iDevice++){

// 		//switch device
// 		GPUerrchk(cudaSetDevice( iDevice ));

// 		//for metric
// 	    	GPUerrchk( cudaMalloc((void **)&pDeviceVectorMetric_[iDevice], 
// 	    			      matrixSizeMetric * sizeof(double)) );
// 		//copy metric 
// 		GPUerrchk( cudaMemcpy(pDeviceVectorMetric_[iDevice], metric,
// 				      matrixSizeMetric * sizeof(double), cudaMemcpyHostToDevice) );
		
// 	    	GPUerrchk( cudaMalloc((void **)&pDeviceVectorBare_[iDevice], 
// 	    			      maxMatrixSizePerDevice * sizeof(double)) );
// 	    	GPUerrchk( cudaMalloc((void **)&pDeviceVectorDressed_[iDevice], 
// 	    			      maxMatrixSizePerDevice * sizeof(double)) );

// 	    }//iDevice

	    


// 	}//transformFunctorGPU

// 	~transformFunctorGPU() {
	    
// 	    //release memory allocate on GPU
// 	    for(int iDevice = 0; iDevice < nDevices_; iDevice++){
// 		std::cout << " cleaning up memory on device: " << iDevice << std::flush << std::endl;
// 		//switch device
// 		GPUerrchk(cudaSetDevice( iDevice ));


// 		GPUerrchk( cudaFree(pDeviceVectorMetric_[iDevice]) );
// 		GPUerrchk( cudaFree(pDeviceVectorBare_[iDevice]) );
// 		GPUerrchk( cudaFree(pDeviceVectorDressed_[iDevice]) );
// 	    }//iDevice

// }
	
// 	void operator() (double *ptr1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &A1_RBSize,
// 			 size_t &A1_BCStart, size_t &A1_BCStop, size_t &A1_CBSize,
// 			 double *ptr2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &A2_RBSize,
// 			 size_t &A2_BCStart, size_t &A2_BCStop, size_t &A2_CBSize);
//     }; //transformFunctorGPU


// #endif



//nuclear derivative contribution to gradient
void nuclear_deriv(Molecule molecule, double *ptr_eg_global);

//back transformation
struct backtransform : boost::noncopyable {
private:
	boost::reference_wrapper<const Matrix> C_;
	boost::reference_wrapper<const Matrix> Cv_;
	size_t N_;
	Thread thread_;

public:

	static
	boost::array<size_t,1> memory(const size_t &N){
		boost::array<size_t,1> memory;
		memory[0] = N*N;
		return memory;
	}

	backtransform(size_t N,
			boost::reference_wrapper<const Matrix> C,
			boost::reference_wrapper<const Matrix> Cv)
	:C_(C), Cv_(Cv), N_(N)
	{
		BOOST_AUTO(memory, this->memory(N_) );
		for (int i = 0; i < memory.size(); ++i) {
			thread_.data[i] = thread_.malloc(memory[i]);
		}
	};

	~backtransform(){thread_.free();};

	void operator()(MapMatrixXd &wmat);

};//backtransform



//struct energy_weighted_density : boost::noncopyable {
//
//private:
//	Thread thread_;
//	size_t nd_,ns_,nv_,nvs_,nl_,no_;
//public:
//
//	//	static
//	//	boost::array<size_t,5> memory(const size_t &nvs, const size_t &nl){
//	//		boost::array<size_t,5> memory;
//	//		memory[0] = nvs*nvs*nl;
//	//		memory[1] = nvs*nvs;
//	//		memory[2] = nvs*nvs*nl;
//	//		memory[3] = nvs*nvs;
//	//		memory[4] = nvs*nvs;
//	//		return memory;
//	//	}
//
//	energy_weighted_density(const size_t &nd, const size_t &ns, const size_t &nv, const size_t &nl)
//:nd_(nd),ns_(ns),nv_(nv),nl_(nl),nvs_(nv+ns),no_(nd+ns)
//{
//		//		BOOST_AUTO(memory, this->memory(nvs_, nl_));
//		//		for (int i = 0; i < memory.size(); ++i) {
//		//			thread_.data[i] = thread_.malloc(memory[i]);
//		//		}
//};
//
//	~energy_weighted_density(){thread_.free();};
//
//	double *operator()(int i){return thread_.data[i];}
//
//	template<typename T, typename U, typename V>
//	void build_core_act(T &eri_ij, T &eri_jk, T &eri_ik, U &wij, V &ea , V &ev, int iocc, int jocc, int kocc){
//
//		double value = 0.0;
//		double value2 = 0.0;
//		double denom = ea(iocc) + ea(jocc);
//
//		for (int a=0; a<nv_; a++){
//			for (int b=0; b<nv_; b++){
//
//				// ( 2*(ia|jb) - (ib|ja) ) * (jb|ka)
//				value += ( 2*eri_ij(a,b) - eri_ij(b,a) ) * eri_jk(b,a)
//														/ ( denom - ev(a) - ev(b) );
//
//
//				// ( 2*(ib|ja) - (ia|jb) ) * (ib|ka)
//				if(iocc != jocc){
//					value2 += ( 2*eri_ij(b,a) - eri_ij(a,b) ) * eri_ik(b,a)
//														/ ( denom - ev(a) - ev(b) );
//				}//(iocc != jocc)
//
//
//			}//b
//		}//a
//
//		wij(kocc,iocc) -= value;
//		//		wij(kocc,iocc) += value/2.0;
//
//		if(iocc != jocc){
//			wij(kocc,jocc) -= value2;
//			//			wij(jocc,kocc) += value2;
//		}
//
//
//	}//build_wij
//
//
//
//	template<typename T, typename U, typename V, typename W>
//	void build_act_act(T &eri_ij, T &eri_jk, T &eri_ik, U &wij, V &ea , W &ev, int iocc, int jocc, int kocc){
//
//		double value = 0.0;
//		double value2 = 0.0;
//		double denom = ea(iocc) + ea(jocc);
//
//		double denom_jk = ea(kocc) + ea(iocc);
//		double denom_ik = ea(iocc) + ea(kocc);
//
//		for (int a=0; a<nv_; a++){
//			for (int b=0; b<nv_; b++){
//
//				// ( 2*(ia|jb)-(ib|ja) ) * (ka|jb)
//				value += ( 2*eri_ij(a,b) - eri_ij(b,a) ) * eri_jk(b,a)
//														/ ( denom - ev(a) - ev(b) );
//
//				// ( 2*(ka|jb)-(kb|ja) ) * (ia|jb)
//				//				value += ( 2*eri_jk(b,a) - eri_jk(a,b) ) * eri_ij(a,b)
//				//						/ ( denom_jk - ev(a) - ev(b) );
//
//				if(iocc != jocc){
//					// ( 2*(ja|ib)-(jb|ia) ) * (ka|ib)
//					value2 += ( 2*eri_ij(b,a) - eri_ij(a,b) ) * eri_ik(b,a)
//															/ ( denom - ev(a) - ev(b) );
//
//					// ( 2*(ka|ib)-(kb|ia) ) * (ja|ib)
//					//					value2 += ( 2*eri_ik(b,a) - eri_ik(a,b) ) * eri_ij(b,a)
//					//							/ ( denom_ik - ev(a) - ev(b) );
//
//				}//(iocc != jocc)
//
//			}//b
//		}//a
//
//		wij(kocc,iocc) -= value/2.0;
//		wij(iocc,kocc) -= value/2.0;
//
//		if(iocc != jocc){
//			wij(kocc,jocc) -= value2/2.0;
//			wij(jocc,kocc) -= value2/2.0;
//		}
//
//	}//build_wij
//
//};//energy_weighted_density









//struct gamma_matrix_c: boost::noncopyable{
//
//private:
//	Thread thread_;
//	size_t nd_,ns_,nv_,nl_,nvs_,no_;
//
//public:
//
//	static
//	boost::array<size_t,2> memory(const size_t &nv, const size_t &nl){
//		boost::array<size_t,2> memory;
//		memory[0] = nv*nv; //eri
//		memory[1] = nv*nv; //t_ij^ab
//
//
//		return memory;
//	}
//
//	gamma_matrix_c(const size_t &nd, const size_t &ns, const size_t &nv, const size_t &nl)
//	:nd_(nd),ns_(ns),nv_(nv),nl_(nl),nvs_(nv+ns),no_(nd+ns)
//	{
//		BOOST_AUTO(memory, this->memory(nv_, nl_));
//		for (int i = 0; i < memory.size(); ++i) {
//			thread_.data[i] = thread_.malloc(memory[i]);
//		}
//	};
//
//
//	~gamma_matrix_c(){thread_.free();}
//
//	double *operator()(int i){return thread_.data[i];}
//
//	template<typename T, typename U, typename V>
//	void form_gamma(T &eri, T &C_jb_P, T &gamma_ia_P, U &ea, V &ev, int iocc, int jocc){
//
//
//
//		//		double denom = ea[iocc] + ea[jocc];
//		//		for(int a=0; a < nv_; a++){
//		//
//		//			for(int P=0; P < nl_; P++){
//		//			for(int b=0; b< nv_; b++){
//		//
//		//				double de = denom - ev[a] - ev[b];
//		//				double t_ij_ab = ( 2*eri(a,b) - eri(b,a) )/de;
//		//
//		//
//		//					gamma_ia_P(a,P) += t_ij_ab * C_jb_P(b,P);
//		//
//		//
//		//			}//b
//		//
//		//			}//P
//		//			}//a
//
//
//		Eigen::MatrixXd tij(nv_,nv_);
//		double denom = ea[iocc] + ea[jocc];
//		for(int a=0; a < nv_; a++){
//
//			for(int b=0; b< nv_; b++){
//
//				double de = denom - ev[a] - ev[b];
//				tij(a,b) = ( 2*eri(a,b) - eri(b,a) )/de;
//
//
//				//					gamma_ia_P(a,P) += t_ij_ab * C_jb_P(b,P);
//
//
//			}//b
//		}//a
//
//		gamma_ia_P += tij*C_jb_P;
//
//
//
//
//	}//form_gamma
//
//	template<typename T>
//	void transform_gamma(const T &gamma_ia_P, const T &virt_coeff_mat, T &gamma_inu_P){
//
//		//virt_coeff_mat(nv,N)
//
//		gamma_inu_P = virt_coeff_mat.transpose()*gamma_ia_P;
//
//
//	}//transform_gamma
//
//
//};//gamma_matrix_c



}//namespace detail
}//namespace rimp2
}//namespace cchem


#endif /* SRC_RI_MP2_GRADIENT_HELPERS_HPP_ */
