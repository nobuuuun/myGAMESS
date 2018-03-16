/*
 * metric.cpp
 *
 *  Created on: Nov 4, 2015
 *      Author: luke
 */

#include "metric.hpp"
#include "ri-mp2/ri-integrals.hpp"
#include <math.hpp>

void cchem::rimp2_gradient::detail::Metric::
BuildTwoCenterEri(std::vector<MapMatrixXd> &sph_c){
    
    MapMatrixXd twoc_eri(data_twoc_, nl_, nl_);

    pe_.task().reset(); 
    if(pe_.rank() != 0) goto skip1;

#pragma omp parallel
    if (pe_.node().rank() == 0) {
	
	BOOST_AUTO(const &auxbasis, auxbasis_.get());
	BOOST_AUTO(const &auxshells, auxbasis.shells());
	
	//for possible cartesian to spherical transformation (allocated for each thread)
	double *data_cart_to_sph = NULL;
	if(spherical_){
	    data_cart_to_sph= new double [auxbasis.max().size()*auxbasis.max().size()];
	}
	//spherical_end
	
	detail::Thread::Task<Parallel::Task&> task(pe_.task());
#pragma omp barrier
	
	cchem::ri::AuxiliaryTwoCenterInt< ::rysq::TwoCenterEri > auxiliary_eri(boost::cref(auxbasis));

	while (++task < auxshells.size()) {
	    const Basis::Shell &S = auxbasis.shells().at(task);
	    size_t s = task;
	    
	    for(size_t q= 0; q <= task; ++q){
		
		const Basis::Shell &Q = auxbasis.shells().at(q);
		
		MapMatrixXd twoc_batch(auxiliary_eri(Q,S), S.size(), Q.size());
		
		if(spherical_){
		    
		    MapMatrixXd transformed_batch(data_cart_to_sph,S.sphsize(), Q.sphsize());
		    transformed_batch = sph_c[S.Lmin()].block(0,0,S.size(),S.sphsize()).transpose()*twoc_batch*sph_c[Q.Lmin()].block(0,0,Q.size(),Q.sphsize());
		    twoc_eri.block(S.sphstart(),Q.sphstart(),S.sphsize(),Q.sphsize()) = transformed_batch;
		    twoc_eri.block(Q.sphstart(),S.sphstart(),Q.sphsize(),S.sphsize()) = transformed_batch.matrix().transpose();
		    
		}else{
		    
		    twoc_eri.block(S.start(),Q.start(),S.size(),Q.size()) = twoc_batch;
		    twoc_eri.block(Q.start(),S.start(),Q.size(),S.size()) = twoc_batch.matrix().transpose();
		    
		}//(spherical_)

		
	    } //q
	    
	} //task

	delete [] data_cart_to_sph;

    }//(pe_.node().rank() == 0)

 skip1:;

    pe_.broadcast(data_twoc_, nl_*nl_, 0);
    
    return;
};

void ::cchem::rimp2_gradient::detail::Metric::DecomposeLDLT(){
    MapMatrixXd twoc_eri(data_twoc_, nl_, nl_);
    MapMatrixXd twoc_pqm1(data_pqm1_, nl_, nl_);
    //--------LDLT from eigen--------
    twoc_pqm1 = twoc_eri.ldlt().solve(Eigen::MatrixXd::Identity(nl_,nl_));
    return;
}

void ::cchem::rimp2_gradient::detail::Metric::DecomposeLLT(){
    MapMatrixXd twoc_eri(data_twoc_, nl_, nl_);
    MapMatrixXd twoc_pqm1(data_pqm1_, nl_, nl_);
    //--------LLT from eigen--------
    Eigen::LLT<Eigen::MatrixXd> cho_twoc_eri(twoc_eri); //.solver();
    twoc_pqm1 = cho_twoc_eri.matrixL();
    return;
};

void ::cchem::rimp2_gradient::detail::Metric::InvertL(){
    MapMatrixXd twoc_pqm1(data_pqm1_, nl_, nl_);
    twoc_pqm1 = twoc_pqm1.inverse();
    return;
};










void ::cchem::rimp2_gradient::detail::Metric::
pseudoInvertedSqrt(const double &tol, double *pPqm12, double *pTwoC, bool doPrint){
    
    if(doPrint)std::cout << "Performing Pseudo-Inversion" << std::endl;
    MapMatrixXd TwoC(pTwoC,nl_,nl_);
    MapMatrixXd Pqm12(pPqm12,nl_,nl_);

    Pqm12 = TwoC; //lm1 is temp here

    Eigen::VectorXd evalues(nl_);
    Eigen::MatrixXd evectors(nl_,nl_);
    Eigen::MatrixXd evectorsCopy(nl_,nl_);

#if !HAVE_MKL
    //get eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Pqm12);
    if(eigensolver.info())
	std::cout << "!!!something went wrong with eigensolver!!! " 
		  << __FILE__ << ":" << __LINE__ << std::endl;
    evalues  = eigensolver.eigenvalues();
    evectors = eigensolver.eigenvectors();
#endif

    
#if HAVE_MKL
    evectors = TwoC;
    int info;
    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L' , nl_, evectors.data(),  nl_, evalues.data() ); 
#endif

// #if HAVE_OPENBLAS
//     evectors = twoc_eri;
//     int info;
//     LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L' , nl_, evectors.data(),  nl_, evalues.data() ); 
// #endif


    evectorsCopy = evectors;

    double maxJ = evalues(nl_-1);

    //eigenvectors with eigenvalues greater than the threshold 'tol' contribute
    int sig = 0;
    for (int ind=0; ind < nl_; ind++){
	if( (evalues(ind) / maxJ) < tol || evalues(ind) <= 0.0 ){
	    sig++;
	    evalues(ind) = 0;
	}else{
	    evalues(ind) = 1/ ( sqrt( evalues(ind) ) );
	}//logic
	for(int jnd = 0; jnd < nl_; jnd++){
	    evectors(jnd,ind) = evectors(jnd,ind)*evalues(ind);
	}//tolerance
	
    }//ind
    if(doPrint)std::cout << "   removed " << sig << " eigenvalues from pseudo-inversion" << std::endl;
    
    //construct (P|Q)^-0.5
    Pqm12 = evectorsCopy*evectors.transpose();

} //pseudoInversion





void ::cchem::rimp2_gradient::detail::Metric::
choleskyInversion(double *data_pqm1, double *data_twoc, bool doPrint){

    MapMatrixXd twoc_eri(data_twoc,nl_,nl_);

#ifdef HAVE_ATLAS
    clapack_dpotrf(CblasColMajor,CblasUpper,
		   (int)nl_,data_twoc,(int)nl_);
    
    clapack_dtrtri(CblasColMajor,CblasUpper,CblasNonUnit,
		   (int)nl_,data_twoc,(int)nl_);
#endif
    
#ifdef HAVE_MKL
    LAPACKE_dpotrf(LAPACK_COL_MAJOR,'U',
		   (int)nl_,data_twoc,(int)nl_);
    
    LAPACKE_dtrtri(LAPACK_COL_MAJOR,'U','N',
		   (int)nl_,data_twoc,(int)nl_);
#endif

#ifdef HAVE_OPENBLAS
            // LAPACKE_dpotrf(LAPACK_COL_MAJOR,'U',
            // 		   (int)nl_,data_twoc,(int)nl_);
    
            // LAPACKE_dtrtri(LAPACK_COL_MAJOR,'U','N',
    	    // 	       (int)nl_,data_twoc,(int)nl_);

    //have problems with openblas lapack --> use eigen
    Eigen::MatrixXd twoc_pqm1(nl_, nl_);
    twoc_pqm1.setZero();
    //--------LLT from eigen--------
    Eigen::LLT<Eigen::MatrixXd> cho_twoc_eri(twoc_eri); //.solver();
    twoc_pqm1 = cho_twoc_eri.matrixL();
    twoc_pqm1 = twoc_pqm1.inverse();
    twoc_eri = twoc_pqm1.transpose();

#endif
    

	// this->pseudoInversion(1e-10,data_pqm1,data_twoc);
	// for(int i = 0; i < 10; i ++ )
	// std::cout << nl_ << std::endl;
	// std::cout << twoc_eri.block(0,0,5,5) << std::endl << std::endl;

	//zap the lower off-diagonal
    twoc_eri = twoc_eri.triangularView<Eigen::Upper>();

	//	std::cout << twoc_eri.block(0,0,5,5) << std::endl << std::endl;
	//	twoc_eri.setZero();
	//	twoc_eri.setIdentity();





    	//form (P|Q)-1 <--- Lm1 * Lm1^T
    	//this is done locally (not on GPU), we use (P|Q)-1 with GPU
#if HAVE_CUBLAS
	// MapMatrixXd pqm1(data_pqm1,nl_,nl_);

	// Eigen::MatrixXd temp(nl_,nl_);
	// temp.setZero();

	// temp =  pqm1*pqm1;
	// pqm1 = temp;
	
	//	pqm1 = twoc_eri*twoc_eri.transpose();

	 // cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
	 // 	    nl_,nl_,nl_,
	 // 	    1.0, data_twoc, nl_,
	 // 	    data_twoc, nl_,
	 // 	    0.0, data_pqm1, nl_);


	// temp = 0.5*pqm1 + 0.5*pqm1.transpose();

	
	// pqm1 = temp;
	
    if(doPrint)std::cout << std::endl
			 << "Building inverse Coulomb metric from Cholesky vectors:"
			 << std::endl << "    (P|Q)^-1 = U^-1 * [U^-1]^T"
			 << std::endl << std::endl;
    
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,
		nl_,nl_,nl_,
		1.0, data_twoc, nl_,
		data_twoc, nl_,
		0.0, data_pqm1, nl_);

#endif

    	//put Lm1 --> data_pqm1
#if !HAVE_CUBLAS

    if(doPrint)std::cout << std::endl
			 << "Using inverse Cholesky decomposed Coulomb metric L^-1"
			 << std::endl << std::endl;
    memcpy(data_pqm1, data_twoc, nl_*nl_*sizeof(double) );

#endif
    

}//choleskyInversion






//--------Direct Eigen Inversion--------
//            twoc_pqm1.setZero();
//            twoc_pqm1 = twoc_eri.inverse(); //inverse for gradient

//            std::cout << "inverted " << std::endl << twoc_pqm1.block(0,0,5,5) << std::endl << std::endl;

//            Eigen::MatrixXd unit(nl,nl);
//            unit = Eigen::MatrixXd::Identity(nl,nl);
//            double error_value = ((twoc_eri*twoc_pqm1) - unit).norm();
//            std::cout << "unit " << std::endl << (twoc_eri*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
//            std::cout << "error_value: " << error_value << std::endl;
//--------------------------------------




//----------------LU factorization with sequential inversion----------------
//            	twoc_pqm1 = twoc_eri;
//                int *IPIV = new int[nl+1];
//
//                int info2 = clapack_dgetrf(CblasColMajor, (int)nl, (int)nl, twoc_pqm1.data(), (int)nl, IPIV);
//                int info3 = clapack_dgetri(CblasColMajor, (int)nl, twoc_pqm1.data(), (int)nl, IPIV);
//
//
//                delete IPIV;
//
//                std::cout << "inverted " << std::endl << twoc_pqm1.block(0,0,5,5) << std::endl << std::endl;
//
//                Eigen::MatrixXd unit(nl,nl);
//                unit = Eigen::MatrixXd::Identity(nl,nl);
//                double error_value = ((twoc_eri*twoc_pqm1) - unit).norm();
//                std::cout << "unit " << std::endl << (twoc_eri*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
//                std::cout << "error_value: " << error_value << std::endl;
//------------------------------------------------



//            //--------Using atlas libaray exclisively to perform cholesky -> invert ( L -> L^-1 )--------
//            twoc_pqm1 = twoc_eri;
//
//            int info4 =  clapack_dpotrf(CblasColMajor, CblasLower, (int)nl, data_real_pqm1, (int)nl);
//
//
////can invert with eigen
//            twoc_pqm1 = twoc_pqm1.triangularView<Eigen::Lower>();
//            twoc_pqm1 = twoc_pqm1.inverse();
//            twoc_pqm1 = twoc_pqm1.transpose()*twoc_pqm1;
//
////            //invert with atlas
////            int info5 = clapack_dpotri(CblasColMajor, CblasLower, (int)nl, data_real_pqm1, (int)nl);
////            std::cout << "info5 " << info5 << std::endl;
////
////            std::cout << "inverted" << std::endl <<twoc_pqm1.block(0,0,5,5) << std::endl << std::endl;
////            //need to symmeterize matrix
////            for (int i =0; i < nl; i++){
////                for (int j=0; j < i; j++){
////                	twoc_pqm1(j,i) = twoc_pqm1(i,j);
////                }//j
////            }//i
//
//                            Eigen::MatrixXd unit(nl,nl);
//                            unit = Eigen::MatrixXd::Identity(nl,nl);
//                            double error_value = ((twoc_eri*twoc_pqm1) - unit).norm();
//                            std::cout << "unit " << std::endl << (twoc_eri*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
//                            std::cout << "error_value: " << error_value << std::endl;
//
//            //----------------------------------------------------------------



//--------here we're eigen SVD--------
// A    = U * D    * V^T
// A^-1 = V * D^-1 * U^T

//            Eigen::MatrixXd tester(nl,nl);
//            twoc_pqm1 = twoc_eri;
//            Eigen::JacobiSVD<Eigen::MatrixXd> svd(twoc_pqm1, Eigen::ComputeThinU | Eigen::ComputeThinV);
//            tester = svd.matrixV()*svd.singularValues().asDiagonal().inverse()*svd.matrixU().transpose();
//            std::cout << (tester*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
//            std::cout << "the singular values are: " << std::endl << svd.singularValues() << std::endl;

////            Eigen::MatrixXd unit(nl,nl);
////            unit = Eigen::MatrixXd::Identity(nl,nl);
////            double error_value = ((tester*twoc_pqm1) - unit).norm();
//////            std::cout << "unit " << std::endl << (tester*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
////            std::cout << "error_value: " << error_value << std::endl;
////            twoc_pqm1 = tester;


////            double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
////            std::cout << "condition number is " << cond << std::endl;
//----------------------------------------






//            //--------This is using the LDLT Eigen routine instead of LLT--------
//
//            Eigen::MatrixXd temp(nl_temp,nl_temp);
//            twoc_pqm1 = twoc_eri;
//            temp = twoc_pqm1.ldlt().solve(Eigen::MatrixXd::Identity(nl_temp,nl_temp));
//            twoc_pqm1 = temp;
//
//            std::cout << "temp" << std::endl <<temp.block(0,0,5,5) << std::endl << std::endl;
//
//            Eigen::MatrixXd unit(nl_temp,nl_temp);
//            unit = Eigen::MatrixXd::Identity(nl_temp,nl_temp);
//            double error_value = ((twoc_eri*twoc_pqm1) - unit).norm();
//            std::cout << "unit " << std::endl << (twoc_eri*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
//            std::cout << "error_value: " << error_value << std::endl;
//
//            //----------------------------------------
//
//
//            //----------------------------------------
//            //find 1-norms
//            double ANorm = 0.0;
//            for(int j=0; j<nl_temp; j++){
//            	double col = 0.0;
//            	for(int i=0; i<nl_temp; i++){
//            		col += std::abs(twoc_eri(i,j)); //1
////            		col += twoc_eri(j,i); //infinity
//
//
//            	}//i
//            	if(col > ANorm) ANorm = col;
//            }//j
//
//
//            double Am1Norm = 0.0;
//            for(int j=0; j<nl_temp; j++){
//            	double col = 0.0;
//            	for(int i=0; i<nl_temp; i++){
//            		col += std::abs(twoc_pqm1(i,j)); //1
////            		col += twoc_pqm1(j,i); //inifinity
//
//
//            	}//i
//            	if(col > Am1Norm) Am1Norm = col;
//            }//j
//
//            std::cout << "condition number " << ANorm*Am1Norm << std::endl;
//            std::cout << "condition number " << twoc_eri.norm()*twoc_pqm1.norm() << std::endl;
//
//            //----------------------------------------






//--------cholesky decomposition from ublas/atlas--------
//        #include "lapack.hpp" //must be declare at head
//        long info;
//        long psize = auxbasis.size();
//        dpotrf_("L", &psize, twoc_eri.data(), &psize, &info);  //this may take some more thought later
//        twoc_eri = twoc_eri.triangularView<Eigen::Lower>();

//        int info;
//        int psize = auxbasis.size();
//        char lower = 'L';
//        info = _clapack_dpotrf_(CblasRowMajor,CblasLower, psize, twoc_eri.data(),psize);
//        twoc_eri = twoc_eri.triangularView<Eigen::Lower>();
//----------------------------------------


//--------default eigen inversion for L matrix--------
//slower for some reason
//            twoc_eri = twoc_eri.inverse(); //eigen-slower?
//----------------------------------------



//--------ublas/atlas inversion for L matrix--------
//(this looks like ublas, should probaly fix this somehow???)
//            long info;  //ublas
//            long psize = auxbasis.size(); //ublas
//            dtrtri_("L","N",&psize,twoc_eri.data(),&psize,&info); //ublas

//            int info = clapack_dtrtri( CblasColMajor,CblasLower,CblasNonUnit,(int)nl_temp, twoc_eri.data(), (int)nl_temp); //atlas

//can't do this since it inverts original matrix (not the L matrix)    dpotri_("L",&psize,twoc_eri.data(),&psize,&info);
//-----------------------------------------

//            std::cout << "inverted" << std::endl << twoc_eri.block(0,0,5,5) << std::endl << std::endl;


//            //--------Here I am generating the (P|Q)^-1 matrix from (L^-1)^T * L --------
//            //alternate: invert twoc_eri to twoc_pqm1 for gradient
//            twoc_pqm1 = twoc_eri.transpose()*twoc_eri;
//            std::cout << twoc_pqm1.block(0,0,5,5) << std::endl << std::endl;
//
//
//            Eigen::MatrixXd unit(nl,nl);
//                        unit = Eigen::MatrixXd::Identity(nl,nl);
//                        double error_value = ((tester*twoc_pqm1) - unit).norm();
//                        std::cout << "unit " << std::endl << (tester*twoc_pqm1).block(0,0,5,5) << std::endl << std::endl;
//                        std::cout << "error_value: " << error_value << std::endl;
////                        twoc_pqm1 = tester;
//
//            //--------------------------------------------------------------------------------

