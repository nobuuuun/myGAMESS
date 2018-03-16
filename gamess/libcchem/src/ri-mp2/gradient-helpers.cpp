/*
 * gradient-helper.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: luke
 */


#include <gradient-helpers.hpp>








void ::cchem::rimp2_gradient::detail::wmat_vovo_functor::operator()
(MapMatrixXd &eri,
		MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
		MapMatrixXd &bba, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
		bool &loop_restriction){

	double *ptr_wmat_thread = new double[no_*no_];
	MapMatrixXd wmat_thread(ptr_wmat_thread, no_ ,no_);
	memset(ptr_wmat_thread, 0, no_*no_*sizeof(double) );


#pragma omp for
	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
		size_t a1off = iocc-a1_start;

		//		for(size_t jocc = a2_start; jocc < a2_stop; jocc++){
		//			size_t a2off = jocc-a2_start;

		int a2_start_prime = (loop_restriction ? iocc : a2_start);

		for(size_t jocc = a2_start_prime; jocc < a2_stop; jocc++){
			size_t a2off = jocc-a2_start;

			eri = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
					bba.block(a2off*NA2BRows,0,NA2BRows,NA2BCols).transpose();

			//			for(int a=0; a<nv_;a++){
			//				for(int b=0; b<nv_;b++){
			//					wmat_thread(iocc,jocc) += pmat_(a+no_,b+no_)*eri(a,b);
			//				}//b
			//			}//a

			wmat_thread(iocc,jocc) += pmat_.block(no_,no_,nv_,nv_).cwiseProduct(eri.transpose()).sum();


			//			if(iocc != jocc){
			//				for(int a=0; a<nv_;a++){
			//					for(int b=0; b<nv_;b++){
			//						wmat_thread(jocc,iocc) += pmat_(a+no_,b+no_)*eri(b,a);
			//					}//b
			//				}//a
			//			}//(iocc!=jocc)

			if(iocc != jocc)wmat_thread(iocc,jocc) += pmat_.block(no_,no_,nv_,nv_).cwiseProduct(eri).sum();

		}//jocc
	}//iocc

#pragma omp critical
	wmat_.block(0,0,no_,no_) += wmat_thread;

	delete [] ptr_wmat_thread;

}//wmat_vovo_functor


//	inline --this is bad
//used to compute vooo contributions to wmat[III]
void ::cchem::rimp2_gradient::detail::wmat_vooo_functor::operator()
													(MapMatrixXd &eri,
															MapMatrixXd &bia, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
															MapMatrixXd &bba, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
															bool &loop_restriction){

	double *ptr_wmat_thread = new double[no_*no_];
	MapMatrixXd wmat_thread(ptr_wmat_thread, no_ ,no_);
	memset(ptr_wmat_thread, 0, no_*no_*sizeof(double) );

#pragma omp for collapse(2)
	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
		for(size_t jocc = a3_start; jocc < a3_stop; jocc++){

			size_t a1off = iocc-a1_start;
			size_t a3off = jocc-a3_start;

			eri = bia.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
					bba.block(a3off*NA3BRows,0,NA3BRows,NA3BCols).transpose();

			//			for(int kocc = 0; kocc <= jocc; kocc++){
			//				for(int a = 0; a < nv_; a++){
			//
			////	old				wmat_thread(jocc,kocc) -= 4*pmat_(a+no_,iocc)*eri(a,kocc); //1
			////					wmat_thread(jocc,kocc) -= 4*pmat_(iocc,a+no_)*eri(a,kocc); //1
			//
			////old					wmat_thread(jocc,iocc) += pmat_(kocc,a+no_)*eri(a,kocc); //2
			////					wmat_thread(jocc,iocc) += pmat_(a+no_,kocc)*eri(a,kocc);//2
			//
			////old					wmat_thread(kocc,iocc) += pmat_(jocc,a+no_)*eri(a,kocc);//3
			////					wmat_thread(iocc,kocc) += pmat_(jocc,a+no_)*eri(a,kocc);//3
			//
			//					if(jocc != kocc){
			////						wmat_thread(kocc,jocc) -= 4*pmat_(a+no_,iocc)*( eri(a,kocc) ); //1
			////						wmat_thread(iocc,jocc) += pmat_(kocc,a+no_)*eri(a,kocc);//2
			////						wmat_thread(iocc,kocc) += pmat_(jocc,a+no_)*eri(a,kocc);//3
			//					}//(iocc != jocc)
			//
			//
			//				}//a
			//			}//kocc

			//1
			wmat_thread.row(jocc) -= 4.0*pmat_.block(0,no_,no_,nv_).row(iocc)*eri;


			//2
			//			wmat_thread(iocc,jocc) += pmat_.block(no_,0,nv_,jocc+1).cwiseProduct(eri.block(0,0,nv_,jocc+1)).sum();
			//			wmat_thread(iocc,jocc) += pmat_.block(no_,0,nv_,jocc).cwiseProduct(eri.block(0,0,nv_,jocc)).sum();
			wmat_thread(iocc,jocc) += pmat_.block(no_,jocc,nv_,1).cwiseProduct(eri.block(0,jocc,nv_,1)).sum();
			wmat_thread(iocc,jocc) += 2.0*pmat_.block(no_,0,nv_,jocc).cwiseProduct(eri.block(0,0,nv_,jocc)).sum();



			//3
			//			wmat_thread.block(iocc ,0, 1,jocc+1) += pmat_.block(jocc,no_,1,nv_)*eri.block(0,0,nv_,jocc+1);
			//			wmat_thread.block(iocc ,0, 1,jocc) += pmat_.block(jocc,no_,1,nv_)*eri.block(0,0,nv_,jocc);
			wmat_thread.block(iocc ,jocc, 1,1) += pmat_.block(jocc,no_,1,nv_)*eri.block(0,jocc,nv_,1);
			wmat_thread.block(iocc ,0, 1,jocc) += 2.0*pmat_.block(jocc,no_,1,nv_)*eri.block(0,0,nv_,jocc);


		}//jocc
	}//iocc

#pragma omp critical
	wmat_.block(0,0,no_,no_) += wmat_thread;

	delete [] ptr_wmat_thread;
}//wmat_vooo_functor




//used to compute vvoo contributions to wmat[III]
void ::cchem::rimp2_gradient::detail::wmat_vvoo_functor::operator()
													(MapMatrixXd &eri,
															MapMatrixXd &bij, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
															MapMatrixXd &bab, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols,
															bool &loop_restriction){

	//	void * awmat_thread = NULL;
	//	posix_memalign(&awmat_thread, 16, no_*no_*sizeof(double) );
	//	double *ptr_wmat_thread = new(awmat_thread) double[no_*no_];
	//	MapMatrixXd wmat_thread(ptr_wmat_thread, no_ ,no_);
	//	wmat_thread.setZero();

#pragma omp for schedule(dynamic)
	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){

		for(size_t avir = a2_start; avir < a2_stop; avir++){

			size_t a1off = iocc-a1_start;
			size_t a2off = avir-a2_start;

			eri = bij.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
					bab.block(a2off*NA2BRows,0,NA2BRows,NA2BCols).transpose();

			//			for(int jocc = 0; jocc < no_; jocc++){
			//				for(int b = 0; b < nv_; b++){
			//					wmat_(iocc,jocc) -= pmat_(avir+no_,b+no_)*( 2*eri(jocc,b) );
			//				}//b
			//			}//jocc

			//#pragma omp critical
			//						wmat_.row(iocc) -= 2.0*pmat_.row(no_+avir)*eri.transpose();
			//			wmat_.block(iocc,0,1,no_) -= 2.0*pmat_.row(no_+avir)*eri.transpose();
			//			wmat_thread.row(iocc) -= 2.0*pmat_.block(no_+avir,no_,1,nv_)*eri.transpose();
			//			wmat_thread.block(iocc,0,1,no_) -= 2.0*pmat_.block(no_+avir,no_,1,nv_)*eri.transpose();

			wmat_.block(iocc,0,1,no_) -= 2.0*pmat_.block(no_+avir,no_,1,nv_)*eri.transpose();


		}//avir
	}//iocc

	//#pragma omp barrier

	//#pragma omp critical
	//	wmat_.block(0,0,no_,no_) += wmat_thread;

	//	delete ptr_wmat_thread;


}//wmat_vvoo_functor



//used to compute oooo contributions to wmat[III]
void ::cchem::rimp2_gradient::detail::wmat_oooo_functor::operator()
													(MapMatrixXd &eri,
															MapMatrixXd &bij, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
															MapMatrixXd &bkl, size_t &a3_start, size_t &a3_stop, size_t &NA3BRows, size_t &NA3BCols,
															bool &loop_restriction){

	double *ptr_wmat_thread = new double[no_*no_];
	MapMatrixXd wmat_thread(ptr_wmat_thread, no_ ,no_);
	memset(ptr_wmat_thread, 0, no_*no_*sizeof(double) );

#pragma omp for
	for(size_t iocc = a1_start; iocc < a1_stop; iocc++){
		size_t a1off = iocc-a1_start;
		//		for(size_t jocc = a3_start; jocc < a3_stop; jocc++){
		//			size_t a2off = jocc-a3_start;

		int a3_start_prime = (loop_restriction ? iocc : a3_start);

		for(size_t jocc = a3_start_prime; jocc < a3_stop; jocc++){
			size_t a2off = jocc-a3_start;

			eri = bij.block(a1off*NA1BRows,0,NA1BRows,NA1BCols)*
					bkl.block(a2off*NA3BRows,0,NA3BRows,NA3BCols).transpose();

			for(int kocc = 0; kocc < no_; kocc++){
				for(int locc = 0; locc < no_; locc++){
					wmat_thread(iocc,kocc) -= pmat_(jocc,locc)*( 2*eri(kocc,locc) );
					wmat_thread(iocc,jocc) += pmat_(kocc,locc)*eri(kocc,locc);
				}//locc
			}//jocc

			if(iocc != jocc){ //i have not tested the index switch , i may have to change it agains
				for(int kocc = 0; kocc < no_; kocc++){
					for(int locc = 0; locc < no_; locc++){
						wmat_thread(jocc,kocc) -= pmat_(iocc,locc)*( 2*eri(locc,kocc) );
						wmat_thread(jocc,iocc) += pmat_(kocc,locc)*eri(locc,kocc);
					}//locc
				}//kocc
			}//(iocc!=jocc)

		}//jocc

	}//iocc

#pragma omp critical
	wmat_.block(0,0,no_,no_) += wmat_thread;

	delete [] ptr_wmat_thread;

}//wmat_oooo_functor










void ::cchem::rimp2_gradient::detail::transform_functor::operator()
					(MapMatrixXd &bij, size_t &a1_start, size_t &a1_stop, size_t &NA1BRows, size_t &NA1BCols,
							MapMatrixXd &bkl, size_t &a2_start, size_t &a2_stop, size_t &NA2BRows, size_t &NA2BCols)
{
#pragma omp for
	for(size_t a = a1_start; a < a1_stop; a++){
		size_t aoff = a-a1_start;

		bij.block(aoff*NA1BRows,0,NA1BRows,NA1BCols) =
				bkl.block(aoff*NA1BRows,0,NA1BRows,NA1BCols)*matrix_;

	}//a
}//transform_functor


void ::cchem::rimp2_gradient::detail::transform_functor_new::operator()
					(double *ptr_A1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
				                     size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
					 double *ptr_A2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
							         size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols)
{
size_t A1_Blocks = A1_BRStop - A1_BRStart;
//Eigen	MapMatrixXd Integrals(ptr_A1, A1_Blocks*NA1BRows , NA1BCols );

size_t A2_Blocks = A2_BRStop - A2_BRStart;
//Eigen	MapMatrixXd Coeff(ptr_A2, A2_Blocks*NA2BRows , NA2BCols );


#pragma omp for
	for(size_t a = A1_BRStart; a < A1_BRStop; a++){
		size_t aoff = a - A1_BRStart;

//Eigen		Coeff.block(aoff*NA1BRows,0,NA1BRows,NA1BCols) =
//Eigen				Integrals.block(aoff*NA1BRows,0,NA1BRows,NA1BCols)*matrix_;

		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
				NA1BRows,NA1BCols,NA1BCols,
				1.0, &ptr_A1[aoff*NA1BRows], A1_Blocks*NA1BRows,
				matrix_, NA1BCols,
				0.0, &ptr_A2[aoff*NA1BRows], A2_Blocks*NA2BRows );


	}//a


//
//	int nThreads = omp_get_num_threads();
//	int myThreadId = omp_get_thread_num();
//
//	size_t myChunkStart;
//	size_t myChunkStop;
//
////	if(myThreadId == 0) {
//		myChunkStart = myThreadId*A1_Blocks*NA1BRows/nThreads;
//		myChunkStop = (myThreadId+1)*A1_Blocks*NA1BRows/nThreads;// + (A1_Blocks*NA1BRows/nThreads > myThreadId);
////	}
////		std::cout << "x " << myChunkStart << " " << myChunkStop << std::endl;
// 		if(myThreadId == nThreads-1)myChunkStop = A1_Blocks*NA1BRows;
////		if( (A1_Blocks*NA1BRows)%nThreads > 0 )
//		size_t width = myChunkStop - myChunkStart;
//
//		std::cout <<"x " << myChunkStart << " " << myChunkStop << " " << width << std::endl;
//
////		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
////				width,NA1BCols,NA1BCols,
////				1.0, &ptr_A1[myChunkStart], A1_Blocks*NA1BRows,
////				matrix_, NA1BCols,
////				0.0, &ptr_A2[myChunkStart], A2_Blocks*NA2BRows );


}//transform_functor_new






void ::cchem::rimp2_gradient::detail::transform_functor_new2::operator()
    (double *ptr_A1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
     size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
     double *ptr_A2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
     size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols)
{
#pragma omp master
    timerTest.reset();
    size_t A1_Blocks = A1_BRStop - A1_BRStart;
//Eigen	MapMatrixXd Integrals(ptr_A1, A1_Blocks*NA1BRows , NA1BCols );

    size_t A2_Blocks = A2_BRStop - A2_BRStart;
//Eigen	MapMatrixXd Coeff(ptr_A2, A2_Blocks*NA2BRows , NA2BCols );


//#pragma omp for
//	for(size_t a = A1_BRStart; a < A1_BRStop; a++){
//		size_t aoff = a - A1_BRStart;
//
////Eigen		Coeff.block(aoff*NA1BRows,0,NA1BRows,NA1BCols) =
////Eigen				Integrals.block(aoff*NA1BRows,0,NA1BRows,NA1BCols)*matrix_;
//
//		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
//				NA1BRows,NA1BCols,NA1BCols,
//				1.0, &ptr_A1[aoff*NA1BRows], A1_Blocks*NA1BRows,
//				matrix_, NA1BCols,
//				0.0, &ptr_A2[aoff*NA1BRows], A2_Blocks*NA2BRows );
//
//
//	}//a



	int nThreads = omp_get_num_threads();
	int myThreadId = omp_get_thread_num();

	size_t myChunkStart;
	size_t myChunkStop;

//	if(myThreadId == 0) {
		myChunkStart = myThreadId*A1_Blocks*NA1BRows/nThreads;
		myChunkStop = (myThreadId+1)*A1_Blocks*NA1BRows/nThreads;// + (A1_Blocks*NA1BRows/nThreads > myThreadId);
//	}
//		std::cout << "x " << myChunkStart << " " << myChunkStop << std::endl;
 		if(myThreadId == nThreads-1)myChunkStop = A1_Blocks*NA1BRows;
//		if( (A1_Blocks*NA1BRows)%nThreads > 0 )
		size_t width = myChunkStop - myChunkStart;

//		std::cout <<"x " << myChunkStart << " " << myChunkStop << " " << width << std::endl;

		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				width,NA1BCols,NA1BCols,
				1.0, &ptr_A1[myChunkStart], A1_Blocks*NA1BRows,
				matrix_, NA1BCols,
				0.0, &ptr_A2[myChunkStart], A2_Blocks*NA2BRows );


#pragma omp barrier
		
#pragma omp master 
		CPUTime += timerTest;

}//transform_functor_new2




// #if HAVE_CUBLAS

// void ::cchem::rimp2_gradient::detail::transformFunctorGPU::operator()
//     (double *ptr_A1, size_t &A1_BRStart, size_t &A1_BRStop, size_t &NA1BRows,
//      size_t &A1_BCStart, size_t &A1_BCStop, size_t &NA1BCols,
//      double *ptr_A2, size_t &A2_BRStart, size_t &A2_BRStop, size_t &NA2BRows,
//      size_t &A2_BCStart, size_t &A2_BCStop, size_t &NA2BCols)
// {

// size_t A1_Blocks = A1_BRStop - A1_BRStart;
// size_t A2_Blocks = A2_BRStop - A2_BRStart;

// // 	int nThreads = omp_get_num_threads();
// // 	int myThreadId = omp_get_thread_num();
	
// // 	size_t myChunkStart = myThreadId*A1_Blocks*NA1BRows/nThreads;
// // 	size_t myChunkStop = (myThreadId+1)*A1_Blocks*NA1BRows/nThreads;// + (A1_Blocks*NA1BRows/nThre

// // 	if(myThreadId == nThreads-1)myChunkStop = A1_Blocks*NA1BRows;

// // 	size_t width = myChunkStop - myChunkStart;

// // 	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
// // 		    width,NA1BCols,NA1BCols,
// // 		    1.0, &ptr_A1[myChunkStart], A1_Blocks*NA1BRows,
// // 		    matrix_, NA1BCols,
// // 		    0.0, &ptr_A2[myChunkStart], A2_Blocks*NA2BRows );

// // #pragma omp barrier



//  const int maxRowsPerDevice = ((rows_/nDevices_) + rows_%nDevices_); // * NA1BRows; 


//  //    std::cout << "rows_ " << rows_ << std::endl;
//          for(int iDevice = 0; iDevice < nDevices_; iDevice++){

// 	     //switch device
// 	     //		std::cout << "iDevice: " << iDevice << std::endl;
// 		GPUerrchk(cudaSetDevice( iDevice )); // iDevice ));

// 		const size_t hostRows = A1_Blocks*NA1BRows;
		
// 		const size_t myChunkStart = iDevice * maxRowsPerDevice;
	
// 		const size_t myChunkStop =  std::min( (size_t) ((iDevice+1)*maxRowsPerDevice) , hostRows );

// 		const size_t myChunkLength = myChunkStop-myChunkStart;

// 		// std::cout << A1_Blocks << " " << NA1BRows << std::endl;
// 		// std::cout<< iDevice <<  " ::: "  << myChunkStart << " " << myChunkStop << " " <<
// 		//     maxRowsPerDevice << " " << hostRows << " " << myChunkLength << std::endl;


// //		size_t dims[2] = {  no*(nv+ns), nl };
// //		size_t chunk[2] = {  nvs, nl };
// 	GPUerrchk( cudaMemcpy2DAsync(pDeviceVectorBare_[iDevice], 
// 				myChunkLength*sizeof(double),
// 				&ptr_A1[myChunkStart],
// 				hostRows*sizeof(double),
// 				myChunkLength*sizeof(double), //width
// 				nl_,      //height
// 				cudaMemcpyHostToDevice) );

// 	double alpha = 1.0;
// 	double beta = 0.0;
// 	cublasDgemm(dVHandles_[iDevice], CUBLAS_OP_N, CUBLAS_OP_N, 
// 		    myChunkLength, nl_, nl_, 
// 		    &alpha, pDeviceVectorBare_[iDevice], myChunkLength, 
// 		    pDeviceVectorMetric_[iDevice],nl_,
// 		    &beta, pDeviceVectorDressed_[iDevice], myChunkLength);

// 	GPUerrchk( cudaMemcpy2DAsync(&ptr_A2[myChunkStart],
// 				hostRows*sizeof(double),
// 				pDeviceVectorDressed_[iDevice], 
// 				myChunkLength*sizeof(double),
// 				myChunkLength*sizeof(double), //width
// 				nl_,      //height
// 				cudaMemcpyDeviceToHost) );

//     }//iDevice


//      for (int iDevice = 0; iDevice < nDevices_; iDevice++){
// 	 GPUerrchk( cudaSetDevice( iDevice ) );
// 	 GPUerrchk( cudaDeviceSynchronize() );
//      }//iDevice

    
// #pragma omp barrier

// }//transformFunctorGPU::operator()



// #endif 








void ::cchem::rimp2_gradient::detail::nuclear_deriv(Molecule molecule, double *ptr_eg_global){
	Eigen::MatrixXd drg(molecule.size(),molecule.size());
	drg.setZero();
	MapMatrixXd eg_global(ptr_eg_global,molecule.size(),3);
	int nat = molecule.size();
	for (int k = 1; k < nat; k++){
		drg(k,k) = (double)0;
		int km1 = k-1;
		for(int l = 0; l < k; l++){
			double rkl = (double)0;
			for (int i = 0; i < 3; i++){
				rkl += (molecule(k,i)-molecule(l,i)) *(molecule(k,i)-molecule(l,i));
			}//i
			drg(k,l) = -1/rkl;
			drg(l,k) = std::sqrt(rkl);
		}//l
	}//i
	Eigen::MatrixXd eg(nat,3); eg.setZero();
	for (int kk = 0; kk < 3; kk++){
		for (int k = 1; k < nat; k++) {
			double zak = (double)molecule.get_atom(k).Z();
			for(int l= 0; l < k; l++){
				double zal = (double)molecule.get_atom(l).Z();
				double pkl = (molecule(k,kk)-molecule(l,kk))/drg(l,k);
				eg_global(k,kk) += pkl*drg(k,l)*zak*zal;
			}//l
		}//k
		for (int k = 0; k < nat-1; k++){
			double zak = (double)molecule.get_atom(k).Z();
			for (int l = k+1; l<nat;l++){
				double zal = (double)molecule.get_atom(l).Z();
				double pkl = (molecule(k,kk)-molecule(l,kk))/drg(k,l);
				eg_global(k,kk) += pkl*drg(l,k)*zak*zal;
			}//l
		}//k
	}//kk

};





void ::cchem::rimp2_gradient::detail::backtransform::operator()(MapMatrixXd &wmat){

	BOOST_AUTO(const &C, C_.get());
	BOOST_AUTO(const &Cv, Cv_.get());

	size_t nbf = C.size2();
	size_t nmo =  C.size1() + Cv.size1();

	size_t no = C.size1();
	size_t nv = Cv.size1();

	MapMatrixXd temp_wmat(thread_.data[0], nbf, nbf );
	temp_wmat.setZero();

	ConstMapMatrixXd virt_coeff_mat(  Cv.data().begin(), Cv.size1(), Cv.size2() );
	ConstMapMatrixXd occ_coeff_mat(  C.data().begin(), C.size1(), C.size2() );

	//  1/2 back transform:
	//transformation is done in two steps since the virtual and occupied LCAO coefficients
	//   are seperated (e.g. C and Cv)
	temp_wmat.block(0,0,nbf,nmo) =  occ_coeff_mat.transpose()*wmat.block(0,0,no,nmo);
	temp_wmat.block(0,0,nbf,nmo) += virt_coeff_mat.transpose()*wmat.block(no,0,nv,nmo);

	// 2/2 back transform:
	//transformation is done in two steps since the virtual and occupied LCAO coefficients
	//   are seperated (e.g. C and Cv)
	//		wmat = occ_coeff_mat.transpose()*temp_wmat.transpose().block(0,0,nmo,nbf);
	//		wmat += virt_coeff_mat.transpose()*temp_wmat.transpose().block(no,0,nmo,nbf);

	wmat = occ_coeff_mat.transpose()*temp_wmat.transpose().block(0,0,no,nbf);
	wmat += virt_coeff_mat.transpose()*temp_wmat.transpose().block(no,0,nv,nbf);

	//symmeterize wmat
	wmat += wmat.transpose().eval();
	wmat = wmat/(double)2;

} //operator()

