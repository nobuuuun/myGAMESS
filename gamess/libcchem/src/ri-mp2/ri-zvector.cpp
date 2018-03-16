/*
 * ri-zvector.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: luke
 */




#include <Eigen/Dense>

#include <vector>
#include <iostream>

#include <ri-zvector.hpp>

void ::cchem::rimp2_gradient::detail::ZVector::solve(MapMatrixXd &lag_mo, MapMatrixXd &pmat,
		MapMatrixXd &zai){

	BOOST_AUTO(const &ea, ea_.get());
	BOOST_AUTO(const &ev, ev_.get());


	lag_mo *= (double)(-1);

        if(pe_.rank() == 0)
	    std::cout << "ITERATION      TESTER      EXPANSION VECTOR NORM      ITERATION TIME"
                      << std::endl
                      << "===================================================================="
                      << std::endl << std::endl;

	
	double *ptr_u = new double[maxitc_*no_*nv_];

	std::vector< MapMatrixXd > u;
	for (size_t iter = 0; iter < maxitc_; iter++){
		u.push_back( MapMatrixXd(ptr_u + iter*no_*nv_, nv_, no_) );
		u[iter].setZero();
	}//iter


	double *ptr_unxt = new double[no_*nv_];
	MapMatrixXd unxt(ptr_unxt,nv_,no_);

	double *ptr_rhs = new double[no_*nv_];
	MapMatrixXd rhs(ptr_rhs,nv_,no_);

	double *ptr_prhs = new double[no_*nv_];
	MapMatrixXd prhs(ptr_prhs,nv_,no_);

	int *ptr_ipivot = new int[maxitc_];
	MapVectorXi ipivot(ptr_ipivot,maxitc_);


	//			cchem::rimp2_gradient::detail::ZVector zvector(nd_,ns_,nv_,nl_,maxitc);
	//				Eigen::VectorXi ipivot(maxitc_); //for lu pivot

	//get zeroth order estimate
	this->sym_eig(lag_mo, ea, ev, u); //B dotted here

	//start iterations
	utility::timer timer;
	for(int z_iter = 0; z_iter < maxitc_; z_iter++){
	    timer.reset();
		//generate orbital hessian trial vector product (need OV, OO, and VV)
		//build unxt
		unxt.setZero();

//		this->build_orbital_hessian_old(u[z_iter], unxt);

		this->build_orbital_hessian(u[z_iter], unxt);

		//modify unxt
		this->sym_eig(ea, ev, unxt, z_iter);

		this->build_uau(z_iter, u, unxt);

		this->form_alpha(z_iter, u);

		this->lu(z_iter, ipivot);

		this->lus(z_iter, ipivot);

		prhs = rhs;
		this->form_new_solution(z_iter, u, rhs);

		//check for convergence of solution vectors
		if(z_iter > 0){
			prhs -= rhs;
			double gnorm = std::sqrt(  prhs.cwiseProduct(prhs).sum()/(no_*nv_) );
			double gmax = std::max(gnorm,0.0);
                        if(pe_.rank() == 0)std::cout << std::setw(4) << z_iter+1
                                                     << "        " 
						     << std::fixed
                                                     << std::setprecision(10)
                                                     << gnorm << std::endl;
			if(gmax < uconv_){
			    if(pe_.rank() == 0)std::cout << std::endl
							 << "      CONVERGED - WAVEFUNCTIONS STATIONARY"
							 << std::endl << std::endl;
				break;
			}//(gmax < uconv)
		}//(z_iter > 0)

		this->update_unxt(z_iter, prhs, unxt, u);

		//check for convergence of next expansion vector
		double gnorm = std::sqrt(  unxt.cwiseProduct(unxt).sum()/(no_*nv_) );

                if(pe_.rank() == 0)std::cout  << std::setprecision(15)
                                              << "                             "
					      << std::fixed
                                              << gnorm << "       "
                                              << timer << std::endl;

		double gmax = std::max(gnorm,0.0);
		if(gmax < small_){
		    if(pe_.rank() == 0)std::cout << std::endl
						 << "      CONVERGED - NEW EXPANSION VECTOR NEGLIGIBLE"
						 << std::endl << std::endl;
			break;
		}//(gmax < small)

		//add unxt to u[z_iter+1]
		u[z_iter+1] = unxt;

	}//z_iter

	//set virtual-occupied part of P matrix
	pmat.block(no_,0,nv_,no_) = rhs;
	pmat.block(0,no_,no_,nv_) = rhs.transpose();
//	zai = rhs; //temp

	zai.block(no_,0,nv_,no_) = rhs;
	zai.block(0,no_,no_,nv_) = rhs.transpose();


	delete [] ptr_ipivot;
	delete [] ptr_u;
	delete [] ptr_unxt;
	delete [] ptr_prhs;
	delete [] ptr_rhs;

}//::cchem::rimp2_gradient::detail::RIZVector::solve



//void ::cchem::rimp2_gradient::detail::RIZVector::build_orbital_hessian(MapMatrixXd &u, MapMatrixXd &unxt){
//
//	double *ptr_bia = thread_.data[0?]; //these are not existent now
//	double *ptr_bjb = thread_.data[1?];
//	double *ptr_eri = thread_.data[2?];
//
//	MapMatrixXd bia(ptr_bia,nv_,nl_);
//	MapMatrixXd bjb(ptr_bjb,nv_,nl_);
//	MapMatrixXd eri(ptr_eri,nv_,nv_);
//
//	for(int iocc = 0; iocc < no_; iocc++){
//		size_t start[] = { iocc*nv_, 0 };
//		size_t finish[] = { (iocc+1)*nv_, nl_ };
//		CIAQ_->get(ptr_bia, start, finish); //bja
//
////		for(int jocc = 0; jocc < no_; jocc++){
//			for(int jocc = iocc; jocc < no_; jocc++){ //avoid some I/O
//			size_t start2[] = { jocc*nv_, 0 };
//			size_t finish2[] = { (jocc+1)*nv_, nl_ };
//			BARE_IA_->get(ptr_bjb, start2, finish2); //bia
//
//			eri = bia*bjb.transpose(); //(ia|jb)
//
//			for(int a = 0; a < nv_; a++){
//				for(int b = 0; b < nv_; b++){
//					unxt(a,iocc) += 4.0*eri(a,b)*u(b,jocc);
//					unxt(a,iocc) -= eri(b,a)*u(b,jocc);
//				}//b
//			}//a
//
//			if(iocc != jocc){
//				for(int a = 0; a < nv_; a++){
//					for(int b = 0; b < nv_; b++){
//						unxt(b,jocc) += 4.0*eri(a,b)*u(a,iocc);
//						unxt(b,jocc) -= eri(b,a)*u(a,iocc);
//					}//b
//				}//a
//			}//(iocc != jocc)
//
//		}//jocc
//	}//iocc
//
//
//	//just reuse buffers/pointers
//	MapMatrixXd bij(ptr_bia,no_,nl_);
//	MapMatrixXd bab(ptr_bjb,nv_,nl_);
//	MapMatrixXd eri2(ptr_eri,no_,nv_);
//
//	for(int iocc = 0; iocc < no_; iocc++){
//		size_t start[] = { iocc*no_, 0 };
//		size_t finish[] = { (iocc+1)*no_, nl_ };
//		CIJQ_->get(ptr_bia, start, finish); //bij
//
//		for(int a = 0; a < nv_; a++){
//			size_t start2[] = { a*nv_, 0 };
//			size_t finish2[] = { (a+1)*nv_, nl_ };
//			BARE_AB_->get(ptr_bjb, start2, finish2); //bab
//
//			eri2 = bij*bab.transpose();
//
//			for(int jocc =0; jocc < no_; jocc++){
//				for(int b = 0; b < nv_; b++){
//					unxt(a,iocc) -= eri2(jocc,b)*u(b,jocc);
//				}//b
//			}//jocc
//
//		}//a
//	}//iocc
//
//};//build_orbital_hessian








//void ::cchem::rimp2_gradient::detail::RIZVector::build_orbital_hessian_old(MapMatrixXd &u, MapMatrixXd &unxt){
//
//	//	double *ptr_bia = thread_.data[0?]; these are nonexistent
//	//	double *ptr_bjb = thread_.data[1?];
//	//	double *ptr_eri = thread_.data[2?];
//
//	void * aptr_bia = NULL;
//	posix_memalign(&aptr_bia, 16, MaxNCIA_*nv_*nl_*sizeof(double) );
//	double *ptr_bia = new(aptr_bia) double[MaxNCIA_*nv_*nl_];
//
//	void * aptr_bjb = NULL;
//	posix_memalign(&aptr_bjb, 16, MaxNCIA_*nv_*nl_*sizeof(double) );
//	double *ptr_bjb = new(aptr_bjb) double[MaxNCIA_*nv_*nl_];
//
//	//	void * aptr_eri = NULL;
//	//	posix_memalign(&aptr_eri, 16, nv_*nv_*sizeof(double) );
//	//	double *ptr_eri = new(aptr_eri) double[nv_*nv_];
//
//	//	MapMatrixXd bia(ptr_bia,nv_,nl_);
//	//	MapMatrixXd bjb(ptr_bjb,nv_,nl_);
//	//	MapMatrixXd eri(ptr_eri,nv_,nv_);
//	//
//	//	for(int iocc = 0; iocc < no_; iocc++){
//	//		size_t start[] = { iocc*nv_, 0 };
//	//		size_t finish[] = { (iocc+1)*nv_, nl_ };
//	//		CIAQ_->get(ptr_bia, start, finish); //bja
//	//
//	//		//		for(int jocc = 0; jocc < no_; jocc++){
//	//		for(int jocc = iocc; jocc < no_; jocc++){ //avoid some I/O
//	//			size_t start2[] = { jocc*nv_, 0 };
//	//			size_t finish2[] = { (jocc+1)*nv_, nl_ };
//	//			BARE_IA_->get(ptr_bjb, start2, finish2); //bia
//	//
//	//			eri = bia*bjb.transpose(); //(ia|jb)
//	//
//	//			for(int a = 0; a < nv_; a++){
//	//				for(int b = 0; b < nv_; b++){
//	//					unxt(a,iocc) += 4.0*eri(a,b)*u(b,jocc);
//	//					unxt(a,iocc) -= eri(b,a)*u(b,jocc);
//	//				}//b
//	//			}//a
//	//
//	//			if(iocc != jocc){
//	//				for(int a = 0; a < nv_; a++){
//	//					for(int b = 0; b < nv_; b++){
//	//						unxt(b,jocc) += 4.0*eri(a,b)*u(a,iocc);
//	//						unxt(b,jocc) -= eri(b,a)*u(a,iocc);
//	//					}//b
//	//				}//a
//	//			}//(iocc != jocc)
//	//
//	//		}//jocc
//	//	}//iocc
//
//
//
//
//
//
////	omp_set_num_threads(8);
//#pragma omp parallel
//	{
//
////		std::cout << omp_get_num_threads() << std::endl;
//		//the following are thread buffers
//		void * aptr_eri = NULL;
//		posix_memalign(&aptr_eri, 16, nv_*nv_*sizeof(double) );
//		double *ptr_eri = new(aptr_eri) double[nv_*nv_];
//		MapMatrixXd eri(ptr_eri,nv_,nv_);
//
//		void * aptr_unxt_thread = NULL;
//		posix_memalign(&aptr_unxt_thread, 16, nv_*no_*sizeof(double) );
//		double *ptr_unxt_thread = new(aptr_unxt_thread) double[nv_*no_];
//		MapMatrixXd unxt_thread(ptr_unxt_thread,nv_,no_);
//		unxt_thread.setZero();
//
//		for(size_t readc = 0; readc < NCIAReads_; readc++){
//#pragma omp single
//			CIAQ_->get(ptr_bia,  CIAReadPos_[readc].first, CIAReadPos_[readc].second);
//			size_t c_start = CIABlocks_[readc].first;
//			size_t c_stop = CIABlocks_[readc].second;
//			size_t c_length = c_stop - c_start;
//			MapMatrixXd bia2(ptr_bia,c_length*nv_,nl_);
//
//			for(size_t readia=readc; readia< NCIAReads_; readia++){ //avoid some I/O:
//#pragma omp single
//				BARE_IA_->get(ptr_bjb, CIAReadPos_[readia].first, CIAReadPos_[readia].second);
//				size_t ia_start = CIABlocks_[readia].first;
//				size_t ia_stop = CIABlocks_[readia].second;
//				size_t ia_length = ia_stop - ia_start;
//				MapMatrixXd bjb2(ptr_bjb,ia_length*nv_,nl_);
//
//#pragma omp for schedule(dynamic,1)
//				for(size_t iocc = c_start; iocc < c_stop; iocc++ ){
//					size_t ioff = iocc - c_start;
//
//					int jstart = (readc == readia ? iocc : ia_start );
//
//					for(int jocc = jstart; jocc < ia_stop; jocc++){ // jocc >= iocc
//						size_t joff = jocc - ia_start;
//
//						eri = bia2.block(ioff*nv_,0,nv_,nl_)
//							*bjb2.block(joff*nv_,0,nv_,nl_).transpose(); //(ia|jb)
//
//						for(int a = 0; a < nv_; a++){
//							for(int b = 0; b < nv_; b++){
////								unxt(a,iocc) += 4.0*eri(a,b)*u(b,jocc);
////								unxt(a,iocc) -= eri(b,a)*u(b,jocc);
//								unxt_thread(a,iocc) += 4.0*eri(a,b)*u(b,jocc);
//								unxt_thread(a,iocc) -= eri(b,a)*u(b,jocc);
//							}//b
//						}//a
//
//						if(iocc != jocc){
//							for(int a = 0; a < nv_; a++){
//								for(int b = 0; b < nv_; b++){
////									unxt(b,jocc) += 4.0*eri(a,b)*u(a,iocc);
////									unxt(b,jocc) -= eri(b,a)*u(a,iocc);
//									unxt_thread(b,jocc) += 4.0*eri(a,b)*u(a,iocc);
//									unxt_thread(b,jocc) -= eri(b,a)*u(a,iocc);
//								}//b
//							}//a
//						}//(iocc != jocc)
//
//					}//jocc
//				}//iocc
//
//			}//readia
//		}//readc
//
//#pragma omp critical
//		unxt += unxt_thread;
//
//		delete ptr_unxt_thread;
//		delete ptr_eri;
//
//	}//omp parallel
//
//	delete ptr_bia;
//	delete ptr_bjb;
//
//
//
//	//		//just reuse buffers/pointers
//	//		MapMatrixXd bij(ptr_bia,no_,nl_);
//	//		MapMatrixXd bab(ptr_bjb,nv_,nl_);
//	//		MapMatrixXd eri2(ptr_eri,no_,nv_);
//	//
//	//		for(int iocc = 0; iocc < no_; iocc++){
//	//			size_t start[] = { iocc*no_, 0 };
//	//			size_t finish[] = { (iocc+1)*no_, nl_ };
//	//			CIJQ_->get(ptr_bia, start, finish); //bij
//	//
//	//			for(int a = 0; a < nv_; a++){
//	//				size_t start2[] = { a*nv_, 0 };
//	//				size_t finish2[] = { (a+1)*nv_, nl_ };
//	//				BARE_AB_->get(ptr_bjb, start2, finish2); //bab
//	//
//	//				eri2 = bij*bab.transpose();
//	//
//	//				for(int jocc =0; jocc < no_; jocc++){
//	//					for(int b = 0; b < nv_; b++){
//	//						unxt(a,iocc) -= eri2(jocc,b)*u(b,jocc);
//	//					}//b
//	//				}//jocc
//	//
//	//			}//a
//	//		}//iocc
//
//
//
//	//the adjusted size of the AB buffer is MaxNAB_*nv_*nl_
//	//the adjusted size of the C buffer is MaxNC_*no_*nl_
//
//	//new buffers for now:
//	void * aabblock = NULL;
//	posix_memalign(&aabblock, 16, MaxNAB_*nv_*nl_*sizeof(double) );
//	double *ptr_abblock = new(aabblock) double[MaxNAB_*nv_*nl_];
//
//	void * acblock = NULL;
//	posix_memalign(&acblock, 16, MaxNC_*no_*nl_*sizeof(double) );
//	double *ptr_cblock = new(acblock) double[MaxNC_*no_*nl_];
//
////	omp_set_num_threads(8);
//#pragma omp parallel
//	{
//
//		void * aptr_eri = NULL;
//		posix_memalign(&aptr_eri, 16, no_*nv_*sizeof(double) );
//		double *ptr_eri2 = new(aptr_eri) double[no_*nv_];
//		MapMatrixXd eri2(ptr_eri2,no_,nv_);
//
//		for(size_t readc = 0; readc < NCReads_; readc++){
//#pragma omp single
//			CIJQ_->get(ptr_cblock, CReadPos_[readc].first, CReadPos_[readc].second);
//			size_t i_start = CBlocks_[readc].first;
//			size_t i_stop = CBlocks_[readc].second;
//			size_t i_length = i_stop - i_start;
//			MapMatrixXd bij2(ptr_cblock,i_length*no_,nl_);
//
//			for(size_t readab = 0; readab < NABReads_; readab++){
//#pragma omp single
//				BARE_AB_->get(ptr_abblock, ABReadPos_[readab].first, ABReadPos_[readab].second);
//				size_t a_start = ABBlocks_[readab].first;
//				size_t a_stop = ABBlocks_[readab].second;
//				size_t a_length = a_stop - a_start;
//				MapMatrixXd bab2(ptr_abblock,a_length*nv_,nl_);
//
//#pragma omp barrier
//
//#pragma omp for schedule(dynamic,1)
//				for(size_t iocc = i_start; iocc < i_stop; iocc++ ){
//					size_t ioff = iocc - i_start;
//
//					for(size_t a = a_start; a < a_stop; a++){
//						size_t offset = a - a_start;
//
//						eri2 = bij2.block(ioff*no_,0,no_,nl_)*
//								bab2.block(offset*nv_,0,nv_,nl_).transpose();
//
//						for(int jocc =0; jocc < no_; jocc++){
//							for(int b = 0; b < nv_; b++){
//								unxt(a,iocc) -= eri2(jocc,b)*u(b,jocc);
//							}//b
//						}//jocc
//
//					}//a
//				}//iocc
//
//			}//read
//		}//readc
//
//		delete ptr_eri2;
//
//	}//omp parallel
////	omp_set_num_threads(1);
//
//	delete ptr_abblock;
//	delete ptr_cblock;
//
//
//	//	//just reuse buffers/pointers
//	//	MapMatrixXd bij(ptr_bia,no_,nl_); //C_ij^Q
//	//	MapMatrixXd bab(ptr_bjb,nv_*nv_,nl_); // bare_ab this buffer may be taxing
//	//
//	//	size_t start_test[] = { 0, 0 };
//	//	size_t finish_test[] = { no_*no_, nl_ };
//	//	C->get(ptr_bjb, start_test, finish_test);
//	//
//	//	omp_set_num_threads(4);
//	//#pragma omp parallel
//	//	{
//	//
//	//		this->allocate_thread_memory();
//	//
//	//		double *ptr_eri = this->get_pointer(8);
//	//		MapMatrixXd eri2(ptr_eri,no_,nv_);
//	//
//	//
//	//
//	//		for(int iocc = 0; iocc < no_; iocc++){
//	//			size_t start[] = { iocc*no_, 0 };
//	//			size_t finish[] = { (iocc+1)*no_, nl_ };
//	//
//	//#pragma omp single
//	//			CIJQ_->get(ptr_bia, start, finish); //bij
//	//
//	//#pragma omp for
//	//			for(int a = 0; a < nv_; a++){
//	//
//	//				eri2 = bij*bab.block(nv_*a,0,nv_,nl_).transpose();
//	//
//	//				for(int jocc =0; jocc < no_; jocc++){
//	//					for(int b = 0; b < nv_; b++){
//	//						unxt(a,iocc) -= eri2(jocc,b)*u(b,jocc);
//	//					}//b
//	//				}//jocc
//	//
//	//			}//a
//	//
//	//		}//iocc
//	//
//	//	}//omp parallel
//	//
//	//	omp_set_num_threads(1);
//	//
//
//
//};//build_orbital_hessian




//void ::cchem::rimp2_gradient::detail::RIZVector::set_up_block_information(){
//
//	//we need to see how many blocks of CIJQ and BARE_AB can
//	// be loaded into the buffers at the same time
//	size_t NBytes = NMB_*1024*1024;
//
//	size_t ABBlockSize = nv_*nl_*sizeof(double);
//	size_t CBlockSize = no_*nl_*sizeof(double);
//
//	//we'll need at least one BARE_AB block
//	MaxNAB_ = 1;
//	MaxNC_ = 0;
//	NBytes -= ABBlockSize;
//
//	//see how many CijQ blocks we can fit into the buffer
//	while(NBytes > CBlockSize){
//		MaxNC_++;
//		NBytes -= CBlockSize;
//		if(MaxNC_ == no_)break;
//	}
//
//	//see how many BARE_AB blocks we can fit into the buffer
//	while(NBytes > ABBlockSize){
//		MaxNAB_++;
//		NBytes -= ABBlockSize;
//		if( MaxNAB_ == nv_)break;
//	}
//
//	//figure out how many reads will need to occur
//	NABReads_ = nv_/MaxNAB_ + (nv_%MaxNAB_ > 0);
//	NCReads_ = no_/MaxNC_ + (no_%MaxNC_ > 0);
//
//	//		std::cout << std::endl << "we can fit " <<MaxNAB_ << " NAB vecotors "<< std::endl;
//	//		std::cout << "we can fit " <<MaxNC_ << " C vectors "<< std::endl;
//	//		std::cout << "using " << (NMB_*1024*1024-NBytes)/1024/1024 << " MB of " <<  NMB_ << " available Bytes" << std::endl;
//	//		std::cout << "there are " << NBytes/1024/1024 << " MB still available" << std::endl<< std::endl;
//
//	if(NABReads_ == 1){
//		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//		temp.first.push_back(0); temp.first.push_back(0);
//		temp.second.push_back(nv_*nv_); temp.second.push_back(nl_);
//		ABReadPos_.push_back(temp);
//
//		ABBlocks_.push_back( std::pair<size_t,size_t> ( 0,nv_) );
//
//	}else{
//		for(int i =0; i<NABReads_-1; i++){
//			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//			temp.first.push_back( (i*MaxNAB_)*nv_); temp.first.push_back(0);
//			temp.second.push_back( ((i+1)*MaxNAB_)*nv_); temp.second.push_back(nl_);
//			ABReadPos_.push_back(temp);
//
//			//				ABBlocks.push_back(MaxNAB*(i+1) );
//			ABBlocks_.push_back( std::pair<size_t,size_t> ( i*MaxNAB_ ,(i+1)*MaxNAB_) );
//		}//i
//		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//		temp.first.push_back( (NABReads_-1)*MaxNAB_ *nv_); temp.first.push_back(0);
//		temp.second.push_back( nv_*nv_ ); temp.second.push_back(nl_);
//		ABReadPos_.push_back(temp);
//
//		//			ABBlocks.push_back(nv_);
//		ABBlocks_.push_back( std::pair<size_t,size_t> ( (NABReads_-1)*MaxNAB_ ,nv_) );
//	}//(NABReads_ == 1)
//
//	if(NCReads_ == 1){
//		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//		temp.first.push_back(0); temp.first.push_back(0);
//		temp.second.push_back(no_*no_); temp.second.push_back(nl_);
//		CReadPos_.push_back(temp);
//
//		CBlocks_.push_back( std::pair<size_t,size_t> ( 0,no_ ) );
//	}else{
//		for(int i =0; i<NCReads_-1; i++){
//			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//			temp.first.push_back( (i*MaxNC_)*no_); temp.first.push_back(0);
//			temp.second.push_back( ((i+1)*MaxNC_)*no_); temp.second.push_back(nl_);
//			CReadPos_.push_back(temp);
//
//			CBlocks_.push_back( std::pair<size_t,size_t> ( i*MaxNC_,(i+1)*MaxNC_ ) );
//
//		}//i
//		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//		temp.first.push_back( (NCReads_-1)*MaxNC_ *no_); temp.first.push_back(0);
//		temp.second.push_back( no_*no_ ); temp.second.push_back(nl_);
//		CReadPos_.push_back(temp);
//
//		CBlocks_.push_back( std::pair<size_t,size_t> ( (NCReads_-1)*MaxNC_, no_ ) );
//
//	}//(NCReads_ == 1)
//
//
//	//		std::cout << "need to perform " << NABReads_ << " NABReads_" << std::endl;
//	//		std::cout << "need to perform " << NCReads_ << " NCReads_" << std::endl;
//	//		std::cout << "MaxNC_ " << MaxNC_ << std::endl;
//	//		std::cout << "AB reads look like: " << std::endl;
//	//		for(int i = 0; i < NABReads_; i++){
//	//			std::cout << "{ " << ABReadPos_[i].first[0] << " " << ABReadPos_[i].first[1] << " }" << std::endl;
//	//			std::cout << "{ "<< ABReadPos_[i].second[0] << " " << ABReadPos_[i].second[1] << " }" << std::endl;
//	//		}
//	//
//	//		std::cout <<std::endl << "C reads look like: " << std::endl;
//	//		for(int i = 0; i < NCReads_; i++){
//	//			std::cout << "{ "<< CReadPos_[i].first[0] << " " << CReadPos_[i].first[1] << " }" << std::endl;
//	//			std::cout << "{ "<< CReadPos_[i].second[0] << " " << CReadPos_[i].second[1] << " }" << std::endl;
//	//		}
//	//
//	//
//	//		std::cout << "AB blocks look like: " << std::endl;
//	//		for(int i = 0; i < NABReads_; i++){
//	//			std::cout << ABBlocks_[i].first << " : " << ABBlocks_[i].second << std::endl;
//	//		}
//	//		std::cout <<std::endl << "C blocks  look like: " << std::endl;
//	//		for(int i = 0; i < NCReads_; i++){
//	//			std::cout << CBlocks_[i].first << " : " << CBlocks_[i].second << std::endl;
//	//		}
//
//
//	//we need to see how many blocks of CIJQ and BARE_AB can
//	// be loaded into the buffers at the same time (reset NBytes)
//	NBytes = NMB_*1024*1024;
//
//	size_t CIABlockSize = nv_*nl_*sizeof(double);
//
//	//we'll need at least one CijQ and one BARE_IA block (multiply by two since they are teh same length)
//	MaxNCIA_ = 1;
//	NBytes -= 2*CIABlockSize;
//
//
//	//see how many CijQ/BARE_IA blocks we can fit into the buffer
//	while(NBytes > 2*CIABlockSize){
//		MaxNCIA_++;
//		NBytes -= 2*CIABlockSize;
//		if(MaxNCIA_ == no_)break;
//	}
//
//	//figure out how many reads will need to occur
//	NCIAReads_ = no_/MaxNCIA_ + (no_%MaxNCIA_ > 0);
//
//	//		std::cout << std::endl << "we can fit " <<MaxNCIA_ << " NCIA and NIA vectors "<< std::endl;
//	//		std::cout << "using " << (NMB_*1024*1024-NBytes)/1024/1024 << " MB of " <<  NMB_ << " available Bytes" << std::endl;
//	//		std::cout << "there are " << NBytes/1024/1024 << " MB still available" << std::endl<< std::endl;
//
//
//
//
//	if(NCIAReads_ == 1){
//		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//		temp.first.push_back(0); temp.first.push_back(0);
//		temp.second.push_back(no_*nv_); temp.second.push_back(nl_);
//		CIAReadPos_.push_back(temp);
//
//		CIABlocks_.push_back( std::pair<size_t,size_t> ( 0,no_) );
//
//	}else{
//		for(int i =0; i<NCIAReads_-1; i++){
//			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//			temp.first.push_back( (i*MaxNCIA_)*nv_); temp.first.push_back(0);
//			temp.second.push_back( ((i+1)*MaxNCIA_)*nv_); temp.second.push_back(nl_);
//			CIAReadPos_.push_back(temp);
//
//			CIABlocks_.push_back( std::pair<size_t,size_t> ( i*MaxNCIA_ ,(i+1)*MaxNCIA_) );
//		}//i
//		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
//		temp.first.push_back( (NCIAReads_-1)*MaxNCIA_ *nv_); temp.first.push_back(0);
//		temp.second.push_back( no_*nv_ ); temp.second.push_back(nl_);
//		CIAReadPos_.push_back(temp);
//
//		CIABlocks_.push_back( std::pair<size_t,size_t> ( (NCIAReads_-1)*MaxNCIA_ ,no_) );
//	}//(NCIAReads_ == 1)
//
//
//
//	//	std::cout << "need to perform " << NCIAReads_ << " NCIAReads_" << std::endl;
//	//	std::cout << "MaxNCIA_ " << MaxNCIA_ << std::endl;
//	//	std::cout << "CIA reads look like: " << std::endl;
//	//	for(int i = 0; i < NCIAReads_; i++){
//	//		std::cout << "{ " << CIAReadPos_[i].first[0] << " " << CIAReadPos_[i].first[1] << " }" << std::endl;
//	//		std::cout << "{ "<< CIAReadPos_[i].second[0] << " " << CIAReadPos_[i].second[1] << " }" << std::endl;
//	//	}
//	//
//	//	std::cout << "CIA blocks look like: " << std::endl;
//	//	for(int i = 0; i < NCIAReads_; i++){
//	//		std::cout << CIABlocks_[i].first << " : " << CIABlocks_[i].second << std::endl;
//	//	}
//
//
//	//if any of the following are zero, there is likely not enough memory allowed for allocation
//	if(MaxNCIA_ == 0 || MaxNAB_ == 0 || MaxNC_ == 0 ){
//		std::cout <<"Error: not enough memory to solve ZVector equations. See file: "
//				<< __FILE__ <<" at line: "<< __LINE__ << std::endl;
//		exit (EXIT_FAILURE);
//	}
//
//
//
//};//set_up_block_information





//void ::cchem::rimp2_gradient::detail::ZVector::lus(const int iter, Eigen::VectorXi &ipivot){
void ::cchem::rimp2_gradient::detail::ZVector::lus(const int iter, MapVectorXi &ipivot){
	MapMatrixXd alpha(thread_.data[4],iter+1,iter+1);
	MapVectorXd cc(thread_.data[2],iter+1);
	MapVectorXd b(thread_.data[0],iter+1); //B
	for (int k =0; k < iter+1; k++) {
		cc(k) = b(k);
	}//k
	//swap coefficient according to ipivot
	for (int k = 0; k < iter+1; k++){
		int l = ipivot(k);
		double xk = cc(l);
		cc(l) = cc(k);
		cc(k) = xk;
	}//k
	for(int k = 0; k < iter+1; k++){
		double xk = cc(k)*alpha(k,k);
		if ( xk != 0.0 ){
			for(int i = k+1; i < iter+1; i++){
				cc(i) -= alpha(i,k)*xk;
			}//i
		}//( xk != 0.0)
		cc(k) = xk;
	}//k
	for(int k = iter; k >= 0; k--){
		double xk = cc(k);
		if( xk != 0.0 ){
			for( int i=0; i < k; i++){
				cc(i) += alpha(i,k)*xk;
			}// i
		}//( xk != 0.0 )
	}//k
};//lus


//	void ::cchem::rimp2_gradient::detail::ZVector::lu(int iter, Eigen::VectorXi &ipivot){
void ::cchem::rimp2_gradient::detail::ZVector::lu(int iter, MapVectorXi &ipivot){
	MapMatrixXd alpha(thread_.data[4],iter+1,iter+1);
	for(int j = 0; j < iter+1; j++){
		//smxpy here : y.block <- a.block*x.block
		if(j > 0){
			alpha.block(j,j,iter+1-j,1) += alpha.block(j,0,iter+1-j,j)*
					alpha.block(0,j,j,1);
		}//(j > 0
		//find pivot
		double temp = 0.0;
		int k;
		for(int i = j; i < iter+1; i++){
			if( std::abs( alpha(j,i) ) > temp ){
				temp = std::abs( alpha(j,i) );
				k = i;
			}//endif
		}//i
		ipivot(j) = k;
		if( temp == 0.0 ){
			std::cout << "LU: singular matrix!!" << std::endl;
			//throw exception here!!!
		}// ( temp == 0.0)
		//swap rows
		alpha.row(j).swap(alpha.row( k ));
		//invert diagonal element
		alpha(j,j) = 1.0/alpha(j,j);
		//sxmpy here : y.block <- x.block*a.block
		if( j > 0 && j < iter){
			alpha.block(j,j+1,1,iter-j) += alpha.block(j,0,1,j)*
					alpha.block(0,j+1,j,iter-j);
		}//( j > 0 && j < iter)
		//scale part of row j
		temp = -alpha(j,j);
		for (int i = j+1; i < iter+1; i++){
			alpha(j,i) *= temp;
		}//i
	}//j
};//lu


void ::cchem::rimp2_gradient::detail::ZVector::form_alpha(int iter, std::vector<MapMatrixXd> &u){
	MapMatrixXd uau(thread_.data[1],maxitc_,maxitc_);
	MapVectorXd uu(thread_.data[3],iter+1);
	MapMatrixXd alpha(thread_.data[4],iter+1,iter+1);
	uu(iter) = u[iter].cwiseProduct(u[iter]).sum();
	for (int i = 0; i<iter+1; i++){
		for (int j = 0; j < iter+1; j++){
			alpha(i,j) = uau(i,j);
		}//j
		alpha(i,i) += uu(i);
	}//i
};//form_alpha


void ::cchem::rimp2_gradient::detail::ZVector::build_uau(int iter, std::vector<MapMatrixXd> &u, MapMatrixXd &unxt){
	MapMatrixXd uau(thread_.data[1],maxitc_,maxitc_);
	uau(iter,iter) = u[iter].cwiseProduct(unxt).sum();
	if(iter > 0){
		for(int j = 0; j < iter; j++){
			//for lower+diagonal part of uau
			uau(iter,j) = u[j+1].cwiseProduct(u[iter]).sum();
			//form upper part of uau
			uau(j,iter) = u[j].cwiseProduct(unxt).sum();
		}//j
	}//(iter > 0)
};//lower_uau


void ::cchem::rimp2_gradient::detail::ZVector::form_new_solution(int iter, std::vector<MapMatrixXd> &u, MapMatrixXd &rhs){
	MapVectorXd cc(thread_.data[2],iter+1);
	rhs.setZero();
	for (int i = 0; i< iter+1; i++){
		rhs += u[i]*cc(i);
	}//i
};//form_new_solution


void ::cchem::rimp2_gradient::detail::ZVector::update_unxt(int iter, MapMatrixXd &prhs, MapMatrixXd &unxt, std::vector<MapMatrixXd> &u){
	MapMatrixXd uau(thread_.data[1],maxitc_,maxitc_);
	MapVectorXd uu(thread_.data[3],iter+1);
	prhs.setZero();
	for (int i = 0; i< iter+1; i++){
		double fac = -uau(i,iter)/uu(i);
		prhs += u[i]*fac;
	}//i
	unxt += prhs;
};//update_unxt


void ::cchem::rimp2_gradient::detail::ZVector::sym_eig(MapMatrixXd &lag, const vector_range &ea, const vector_range &ev, std::vector<MapMatrixXd> &u){
	MapVectorXd b(thread_.data[0],maxitc_); //B
	b.setZero();
	for(int a = 0; a < nv_; a++){
		for(int i = 0; i < no_; i++){
			u[0](a,i) = lag(a,i)/(ev[a]-ea[i]);
		}//i
	}//a
	b(0) = u[0].cwiseProduct(u[0]).sum(); //like matrix dot product
};//sym_eig



void ::cchem::rimp2_gradient::detail::ZVector::sym_eig(const vector_range &ea, const vector_range &ev, MapMatrixXd &unxt, int iter){
	for(int a = 0; a < nv_; a++){
		for(int i = 0; i < no_; i++){
			unxt(a,i) /= (ev[a]-ea[i]);
		}//i
	}//a
};//sym_eig












void ::cchem::rimp2_gradient::detail::JK_RIZVector::build_orbital_hessian(MapMatrixXd &u, MapMatrixXd &unxt)
	{

		utility::timer timer;
		timer.reset();

		utility::timer timerTest;

		utility::timer timer2;

		BOOST_AUTO(const &Ca, Ca_.get());
		BOOST_AUTO(const &Cv, Cv_.get());

		typedef const Eigen::Map<const Eigen::MatrixXd,Eigen::AutoAlign> ConstMapMatrixXd;
		ConstMapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() );
		ConstMapMatrixXd virt_coeff_mat(  Cv.data().begin(), Cv.size1(), Cv.size2() );

		Eigen::MatrixXd Cocc(Ca.size1(), Ca.size2()); //fix this later so we do not need the additional Cocc array
		Cocc = occ_coeff_mat;

		//construct pseudo density

		//update "Right" side (this us the only thing to change each iteration)
		cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
				no_,N_,nv_,
				1.0, u.data(), nv_,
				virt_coeff_mat.data(), nv_,
				0.0, pRight_[0], no_);

		//build pseudo density
		cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
				N_, N_, no_,
				1.0, pLeft_[0], no_,
				pRight_[0], no_,
				0.0, pD_[0], N_);

		size_t num_mn = N_*N_; //no symmetry...


		memset(pJ_[0], 0, N_*N_*sizeof(double) );
		memset(pK_[0], 0, N_*N_*sizeof(double) );

		timerTest.reset();




		std::vector<int> dims;
		dims.push_back(no_); dims.push_back(0);

		timerTest.reset();

		const size_t mnSize = N_*(N_+1)/2;
		size_t NMB = 500;
		NMB =3000;

		#if HAVE_CUBLAS
		size_t GPUMem;
		if(pGPU_ != NULL )
		    GPUMem = ::cchem
			::rimp2_gradient
			::detail
			::GPUBuildJKMatrixFunctor
			::GPUMemLimit(nl_, dims, N_, mnSize, pGPU_->getNumStreams(),1 );
		size_t NMBDeviceTemp = GPUMem;
		for(size_t ip = 0; ip < pe_.size(); ip++){
		    if(ip == pe_.rank()){
			pe_.broadcast(&GPUMem, 1, ip);
		    }else{
			pe_.broadcast(&NMBDeviceTemp, 1, ip);
		    }
		    if(NMBDeviceTemp < GPUMem) GPUMem = NMBDeviceTemp;
		} //ip
		
		//if(pe_.rank() == 0)std::cout << "GPUMem: " << GPUMem << std::endl;
		NMB = GPUMem;

		#endif


		bool ReadRow = 0; bool debug = 0; bool printArray = 0;// pe_.rank() ==0;
		::cchem::rimp2_gradient::detail::DFAsync 
		      JMatrix(NMB,
			      LM1_MN_SYM_, 1, 0, nl_, 0,
			      ReadRow, debug, pe_, printArray);


		//JK on CPU
		if(pGPU_ == NULL)
		    {
		    	::cchem::rimp2_gradient::detail::buildJKMatrixFunctor 
		    	    Functor(mnSize, N_, pJTemp_, pD_, pJ_, pDmn_, pDAux_, 
		    		    pK_, pLeft_, pRight_, dims, JMatrix.getMaxNAB(0) );
		    	JMatrix.DoOpAsync_R_NODE_PARALLEL( Functor );
		    }

		#if HAVE_CUBLAS
		//JK with GPU
		if(pGPU_ != NULL)
		    {//scope
			int myRank = pe_.rank();
			timer2.reset();
			//int nThreads = omp_get_max_threads();
			//omp_set_num_threads(1);

			{
			    ::cchem::rimp2_gradient::detail::GPUBuildJKMatrixFunctor
				Functor(mnSize, N_, no_, pJTemp_, pD_, pJ_, pDmn_, pDAux_,
					pK_, pLeft_, pRight_, dims, JMatrix.getMaxNAB(0), pGPU_,
					myRank);
			    //if(pe_.rank() == 0)
			    //std::cout  << "  build time "<< timer2 << std::endl;
			    //if(pe_.rank() == 0)JMatrix.DoOpAsync_R( Functor );
			    JMatrix.DoOpAsync_R_NODE_PARALLEL( Functor );
			    //JMatrix.DoOpAsync_R_NODE_PARALLEL_GPU( Functor, pGPU_ );
			}
			//omp_set_num_threads(nThreads);			

			//std::cout  << "  JK matrix with I/O: "<< timerTest << std::endl;
		    }//scope
		#endif

		pe_.reduce("+",pJ_[0], (size_t)N_*N_); 
		pe_.reduce("+",pK_[0], (size_t)N_*N_); 

		//4(ia|jb)P_jb = C_im J_mn C_na
//		unxt = virt_coeff_mat*J_[0]*occ_coeff_mat.transpose();
//		unxt *= 4.0;
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
				N_,no_,N_,
				1.0, pJ_[0], N_,
				Ca.data().begin(), no_,
				0.0, pHalf_, N_);

		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				nv_,no_,N_,
				4.0, Cv.data().begin(), nv_,
				pHalf_, N_,
				0.0, unxt.data(), nv_);


		// -(ij|ab)P_jb = C_im K_mn C_ra
//		unxt -= virt_coeff_mat*K_[0]*occ_coeff_mat.transpose();
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
				N_,no_,N_,
				1.0, pK_[0], N_,
				Ca.data().begin(), no_,
				0.0, pHalf_, N_);

		cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,
				N_,no_,N_,
				1.0, pK_[0], N_,
				Ca.data().begin(), no_,
				1.0, pHalf_, N_);

//		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
//				nv_,no_,N_,
//				-1.0, virt_coeff_mat.data(), nv_,
//				pHalf_, N_,
//				1.0, unxt.data(), nv_);

		// -(ib|ja)P_jb = C_in K_nm C_ma
//		unxt -= virt_coeff_mat*K_[0].transpose()*occ_coeff_mat.transpose();

//		cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,
//				N_,no_,N_,
//				1.0, K_[0].data(), N_,
//				occ_coeff_mat.data(), no_,
//				0.0, pHalf_, N_);

		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
				nv_,no_,N_,
				-1.0, Cv.data().begin(), nv_,
				pHalf_, N_,
				1.0, unxt.data(), nv_);

		//if(pe_.rank() == 0)
		    //std::cout << "  Elec. Hessian time: " << timer << std::endl;

	}//build_orbtial_hessian

