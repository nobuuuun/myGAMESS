/*
 * DFIntGradient.cpp
 *
 *  Created on: Feb 12, 2016
 *      Author: luke
 */

#include "utility/timer.hpp"

#include <three-center-work.hpp>
#include <DFIntGradient.hpp>


void ::cchem::rimp2_gradient::DFIntGradient::gradient::compute_DFIntGradient(){

	utility::timer timer;

	//factorize Ptot = pscf+pmat into neg/pos factors
	// e.g P(MO,MO) -> P(MO,+/-) -> P(AO,+/-)
	this->setup();

	if(pe_.rank() == 0)std::cout << "Creating Arrays: "
				     << std::endl;
	//BARE_PN
	size_t dims_pn[2] = {  no_*(Npi_+Nni_), nl_ };
	size_t chunk_pn[2] = {  dims_pn[0], 1 };
	rt_.arrays().allocate<double,2>("rimp2.bare_pn", dims_pn, pe_, chunk_pn);
	Array<double> *BARE_PN_ = rt_.arrays().find< Array<double> >("rimp2.bare_pn");
	rt_.cout() << "     " << *BARE_PN_ << std::endl;

	//MNAUX_FORCE
	size_t dims_MNAUX_FORCE[2] = {  N_*N_, nl_ };
	size_t chunk_MNAUX_FORCE[2] = {  dims_MNAUX_FORCE[0], 1 };
	rt_.arrays().allocate<double,2>("rimp2.MNAUX_FORCE", dims_MNAUX_FORCE, pe_, chunk_MNAUX_FORCE);
	Array<double> *MNAUX_FORCE_ = rt_.arrays().find< Array<double> >("rimp2.MNAUX_FORCE");
	rt_.cout() << "     " << *MNAUX_FORCE_ << std::endl
		   << std::endl;

	timer.reset();
	this->build_three_center_terms(BARE_PN_);
	if(pe_.rank() == 0)
	    std::cout << "JK Gradient: Time for 3 center intermediates: " 
		      << timer << std::endl;

	pe_.barrier();
	
	if(pe_.rank() != 0) goto skip1;
	timer.reset();
	this->TwoC_DERI_contribution();
	if(pe_.rank() == 0)
	    std::cout << "JK Gradient: Time for 2 center contributions: " 
		      << timer << std::endl;
 skip1:;
	pe_.barrier();


	timer.reset();
      	this->ThreeC_DERI_contribution(BARE_PN_, MNAUX_FORCE_);
	if(pe_.rank() == 0)
	    std::cout << "JK Gradient: time for 3 center contributions: " 
		      << timer << std::endl;
	pe_.barrier();

	delete MNAUX_FORCE_;
	delete BARE_PN_;
}



void ::cchem::rimp2_gradient::DFIntGradient::gradient::setup(){


	BOOST_AUTO(const &Ca, Ca_.get());
	BOOST_AUTO(const &Cv, Cv_.get());

	ConstMapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() );
	ConstMapMatrixXd virt_coeff_mat(  Cv.data().begin(), Cv.size1(), Cv.size2() );

	Eigen::MatrixXd C(Ca.size1()+Cv.size1(),N_);
	C.block(0,0,Ca.size1(),N_) = occ_coeff_mat;
	C.block(Ca.size1(),0,Cv.size1(),N_) = virt_coeff_mat;

	Eigen::MatrixXd Ptotal(no_+nv_,no_+nv_);
	Ptotal = Pmo_;

	Eigen::MatrixXd PAO(N_,N_);
	PAO = C.transpose()*Ptotal*C;
	P2AO_.push_back(PAO);

	PAO = Pscf_;
	D2AO_.push_back(PAO);

	Ptotal *= double(0.5);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Ptotal);

	Eigen::VectorXd evalues(no_+nv_);
	Eigen::MatrixXd evectors(no_+nv_,no_+nv_);

	evalues  = eigensolver.eigenvalues();
	evectors = eigensolver.eigenvectors();

	for (int i = 0; i < evalues.size(); i++){
		if (evalues(i) >= 0.0){
			Npi_++;
		}else{
			Nni_++;
		}//(evalues(i) >= 0.0)

	}//i

	Eigen::MatrixXd Pos(no_+nv_,Npi_);
	Eigen::MatrixXd Neg(no_+nv_,Nni_);

	int Pcounter = 0;
	int Ncounter = 0;

	for(int i = 0; i < no_+nv_; i++){

		//you can put a cuttoff here... if(fabs(evalue(i) >= delta){ etc...
		//		if(fabs(evalues(i)) > cutoff){
		if(evalues(i) >= 0.0){
			double val = sqrt(fabs(evalues(i)));
			Pos.block(0,Pcounter,no_+nv_,1) = val*evectors.block(0,i,no_+nv_,1);
			Pcounter++;
		}else{
			double val = sqrt(fabs(evalues(i)));
			Neg.block(0,Ncounter,no_+nv_,1) = -val*evectors.block(0,i,no_+nv_,1);
			Ncounter++;
		}//logic
		//		}//logic
	}//i

	Eigen::MatrixXd posAO(Npi_,N_);
	Eigen::MatrixXd negAO(Nni_,N_);

	posAO = Pos.transpose()*C; //(Npi,N) = (Npi,no_+nv_)*(no_+nv_,N)
	negAO = Neg.transpose()*C; //(Nni,N) = (Nni,no_+nv_)*(no_+nv_,N)

	//	this->print_matrix(Pos);
	//	this->print_matrix(Neg);

	posFactor_.push_back(posAO);
	negFactor_.push_back(negAO);

}


void ::cchem::rimp2_gradient::DFIntGradient::gradient::
setupAsyncMNAuxAccess(size_t &NMB, std::vector <std::pair<size_t,size_t> >&AsShWork,
		       std::vector <std::pair<size_t,size_t> >&AsBasWork, size_t &MaxBuffer){

	BOOST_AUTO(const &auxbasis, auxbasis_.get());
	BOOST_AUTO(const &auxshells, auxbasis.shells());

	//find minimum NMB needed

	size_t MaxShSize = auxbasis.max().size();
	if(spherical_) MaxShSize = auxbasis.max_spherical_shell_size();

	size_t  MinNMB = 2*MaxShSize*N_*N_*sizeof(double)/(1024*1024) +1;

	//adjust minimum if needed
	if(MinNMB > NMB)NMB = MinNMB;

	//bytes needed to store each column
	size_t ColCost = N_*N_*sizeof(double);
	size_t memory;

	//figure out column range in terms of shells for each async read
	//  the size of AsShWork is the number of reads needed
	for(size_t shell = 0; shell < auxbasis.shells().size(); shell++){
		const Basis::Shell &L = auxbasis.shells().at(shell);

		size_t ShSize = L.size();
		if(spherical_){
			ShSize = L.sphsize();
		}

		if(AsShWork.size() == 0){
			memory = NMB*1024*1024;
			memory -= 2*ShSize*ColCost;
			AsShWork.push_back ( std::pair<size_t,size_t> (0,1) );

		}else if( (memory > 2*ShSize*ColCost) ){
			memory -= 2*ShSize*ColCost;
			AsShWork.back().second +=1;

		}else{
			memory = NMB*1024*1024;
			memory -= 2*ShSize*ColCost;
			AsShWork.push_back ( std::pair<size_t,size_t> (shell,shell+1) );
		}//logic

	}//shell

	//	set up auxiliary basis indices, find max buffer needed
	//figure out column range in terms of basis function for each async read
	//  this is needed since the MNAux_ array is address as (N_*N*,nl_)
	for(int i = 0; i < AsShWork.size(); i++){
		const Basis::Shell &L1 = auxbasis.shells().at( AsShWork[i].first  );
		const Basis::Shell &L2 = auxbasis.shells().at( AsShWork[i].second-1 );

		size_t diff;
		if(!spherical_){
		AsBasWork.push_back (std::pair<size_t,size_t> ( L1.start(), L2.start() + L2.size() ));
		diff = L2.start() + L2.size() - (L1.start() );
		}
		if(spherical_){
			AsBasWork.push_back (std::pair<size_t,size_t> ( L1.sphstart(), L2.sphstart() + L2.sphsize() ));
			diff = L2.sphstart() + L2.sphsize() - (L1.sphstart() );
		}

		if( diff > MaxBuffer) MaxBuffer = diff;
	}

}//setupAsyncMNAuxAccess









void ::cchem::rimp2_gradient::DFIntGradient::gradient::
build_three_center_terms(Array<double> * BARE_PN_){

    utility::timer timerTest;
    
	BOOST_AUTO(const &Ca, Ca_.get());

	ConstMapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() );

	DAux_.push_back( Eigen::VectorXd(nl_) );
	DAux_[0].setZero();
	PAux_.push_back( Eigen::VectorXd(nl_) );
	PAux_[0].setZero();

	
	pe_.barrier();

	{
	    //refernce denisty
	    double *pDmn = new double[ N_*(N_+1)/2 ];
	    //correlated density
	    double *pPmn = new double[ N_*(N_+1)/2 ];
	    
	    //morph N^2 densities into N*(N+1)/2 densities
	    size_t munu = 0;
	    const double *D = D2AO_[0].data();
	    const double *P = P2AO_[0].data();
	    for(size_t mu = 0; mu < N_; mu++){
		for(size_t nu = 0; nu <=mu ; nu++, munu++){
		    pDmn[munu] = ( mu==nu ? D[mu+nu*N_] : D[mu+nu*N_] + D[nu+mu*N_] );
		    pPmn[munu] = ( mu==nu ? P[mu+nu*N_] : P[mu+nu*N_] + P[nu+mu*N_] );
		}//nu
	    }//mu
	    
	    timerTest.reset();
	    
	    size_t NMB = 1000;
	    NMB=500;
	    NMB=1000;
	    bool ReadRow = 0; bool debug = 1; bool printArray = 0; //(pe_.rank() == 0);
	    ::cchem::rimp2_gradient::detail::DFAsync 
		  MakeDP(NMB,
			 BARE_MN_SYM_, 1, 0, nl_, 0,
			 ReadRow, debug, pe_, printArray);
	    
	    ::cchem::rimp2_gradient::DFIntGradient::gradient::DAuxPAuxFunctor 
		  Functor(pDmn, pPmn, DAux_[0], PAux_[0], N_);
	    
	    //if(pe_.rank() == 0)
	    //MakeDP.DoOpAsync_R( Functor );
	    MakeDP.DoOpAsync_R_NODE_PARALLEL( Functor );
	    
	    if(pe_.rank() == 0)
		std::cout << "     Time for DAuxPAuxFunctor (CPU): " 
			  << timerTest << std::endl;
	    
	    pe_.barrier();
	    
	    delete [] pDmn;
	    delete [] pPmn;
	}
	
	
	pe_.reduce("+", DAux_[0].data(), nl_);
	pe_.reduce("+", PAux_[0].data(), nl_);
	

	DAux_.push_back( Eigen::VectorXd(nl_) );
	//transform : ( (mu nu|Q).dot.Dmn )*(Q|P)^-1
	DAux_[1] = pqm1_*DAux_[0];
	
	PAux_.push_back( Eigen::VectorXd(nl_) );
	//transform : ( (mu nu|Q).dot.Pmn )*(Q|P)^-1
	PAux_[1] = pqm1_*PAux_[0];
	
	///////////////////////////////
	//   BARE_MN_SYM_ --> BARE_PN_
	// 1)  (mu nu|Q)*C --> (i nu|Q)
	// 
	// 2a) (i nu|Q) * PosAO -> (i,Npi|Q)
	// 2b) (i nu|Q) * NegAO -> (i,Nni|Q)
	{
	    timerTest.reset();
	    
	    size_t NMB = 1000;
	    NMB=500;
		NMB=1000;
		bool read_restrict = 1;
		bool ReadRow = 0;
		bool MagicSwap = 0; bool printArray = 0; //(pe_.rank() == 0);
		bool loopRestrict = 0;
		::cchem::rimp2_gradient::detail::DFAsync 
		      MakePN(NMB,
			     BARE_MN_SYM_, 1, 0, nl_, 0,
			     BARE_PN_, 1, 0, nl_, 0,
			     read_restrict,loopRestrict, ReadRow, MagicSwap, pe_, printArray);

		::cchem::rimp2_gradient::DFIntGradient::gradient::PNFunctor
		      Functor(posFactor_[0],negFactor_[0],N_,no_,Npi_,Nni_,occ_coeff_mat);
		
		MakePN.DoOpAsync_R_W_NODE_PARALLEL( Functor );

		if(pe_.rank() == 0)
		std::cout << "     Time for PNFunctor (CPU): " 
			  << timerTest << std::endl; 
	}
	pe_.barrier();

	//build  B_ia^Q : PN*(mu nu| Q) -> PN*(mu nu|Q)(Q|R)-1
	//3/3 transform BARE_PN_
	{
	    timerTest.reset();

		size_t NMB = 1000;
		//NMB=500;
		bool read_restrict = 1;
		bool ReadRow = 1;
		bool MagicSwap = 1;
		bool loopRestrict = 0;
		::cchem::rimp2_gradient::detail::DFAsync 
		      MakeBIAQ(NMB,
			       BARE_PN_, no_, 0, 1, 0,
			       BARE_PN_, no_, 0, 1, 0,
			       read_restrict, loopRestrict, ReadRow, MagicSwap, pe_);

		::cchem::rimp2_gradient::detail::transform_functor_new 
		      Functor( pqm1_.data() );

		MakeBIAQ.DoOpAsync_R_W_NODE_PARALLEL( Functor );


		if(pe_.rank() == 0)
		    std::cout << "     Time for transform_functor (CPU): " 
			      << timerTest << std::endl; 
	}

	pe_.barrier();


	auxForce_.push_back( Eigen::MatrixXd(nl_,nl_) );
	auxForce_[0].setZero();


	//if(pe_.rank()==0)
	{
	    timerTest.reset();
	    size_t NMB = 5000; int tag = 0;

	    bool printDebug = (pe_.rank() == 0);
	    bool read_restrict = 1; bool ReadRow = 0; bool MagicSwap = 0;
	    bool loopRestrict = 1;
	    ::cchem::rimp2_gradient::detail::DFAsync 
		  makeAuxForce(NMB,
			       BARE_PN_, 1, 0, nl_, 0,
			       BARE_PN_, 1, 0, nl_, 0,
			       read_restrict, loopRestrict, ReadRow, MagicSwap, pe_);
	    
	    ::cchem::rimp2_gradient::DFIntGradient::gradient
		  ::auxForceFunctor functor(
					    N_,no_, nl_,Npi_,Nni_, auxForce_[0]
					    );
	    
	    makeAuxForce.DoOpAsync_RR( functor );
	    
	    if(pe_.rank() == 0)
		std::cout << "     Time for auxForceFunctor (CPU): " 
			  << timerTest << std::endl; 
	    
	}

	pe_.reduce("+",auxForce_[0].data(),nl_*nl_);

}//:cchem::rimp2_gradient::DFIntGradient::gradient::build_three_center_terms


void ::cchem::rimp2_gradient::DFIntGradient::gradient::TwoC_DERI_contribution(){

    BOOST_AUTO(const &auxbasis, auxbasis_.get());
    BOOST_AUTO(const &auxshells, auxbasis.shells());

    
    cchem::ri::AuxiliaryTwoCenterInt< ::rysq::TwoCenterDerivativeEri >
	auxiliary_eri(boost::cref(auxbasis));

    for(size_t task = 0; task < auxshells.size(); task++){
	const Basis::Shell &S = auxbasis.shells().at(task);
	
	for(size_t q= 0; q <= task; ++q){
	    const Basis::Shell &Q = auxbasis.shells().at(q);

	    if( S.atom() == Q.atom() )continue;

	    const double perm = (task == q ? 0.5 : 1.0);
	    
	    int satom = S.atom();
	    int qatom = Q.atom();

	    
	    MapMatrixXd twoc_batch(auxiliary_eri(Q,S), Q.size()*S.size(),6);
	    
	    //transform over each coordinate (e.g. d/dAx)
	    double *ptr_deriv = twoc_batch.data();

	    if(spherical_){
		
		Eigen::MatrixXd transformed_batch(Q.sphsize(),S.sphsize());
		for(int nder = 0; nder < 6; nder++){

		    int center = Q.atom();
		    if(nder < 3)center = S.atom();
		    
		    MapMatrixXd coord_batch(ptr_deriv,S.size(),Q.size());
		    transformed_batch = sph_c_[S.Lmin()].block(0,0,S.size(),S.sphsize()).transpose()*coord_batch*sph_c_[Q.Lmin()].block(0,0,Q.size(),Q.sphsize());
		    ptr_deriv += Q.size()*S.size();
		    
		    const int ider = nder %3;
		    for(int st = 0; st < S.sphsize(); st++){
			for(int qt = 0; qt < Q.sphsize(); qt++){
			    
			    const double val = 0.5*(DAux_[1](Q.sphstart()+qt)*PAux_[1](S.sphstart()+st) + DAux_[1](S.sphstart()+st)*PAux_[1](Q.sphstart()+qt) );
			    //coulomb
			    eg_(center,ider) -= val*transformed_batch(st,qt);
			    //exchange
			    eg_(center,ider) += perm*(auxForce_[0](Q.sphstart()+qt,S.sphstart()+st))*transformed_batch(st,qt);
			    
			}//q
		    }//s
		    
		}//ic
	    }else{//spherical
		
		for(int nder = 0; nder < 6; nder++){
		    
		    int center = Q.atom();
		    if(nder < 3)center = S.atom();
		    
		    MapMatrixXd coord_batch(ptr_deriv,S.size(),Q.size());
		    ptr_deriv += Q.size()*S.size();
		    
		    const int ider = nder %3;
		    for(int st = 0; st < S.size(); st++){
			for(int qt = 0; qt < Q.size(); qt++){
			    
			    const double val = 0.5*(DAux_[1](Q.start()+qt)*PAux_[1](S.start()+st) + DAux_[1](S.start()+st)*PAux_[1](Q.start()+qt) );
			    //coulomb
			    eg_(center,ider) -= val*coord_batch(st,qt);
			    //exchange
			    eg_(center,ider) += perm*(auxForce_[0](Q.start()+qt,S.start()+st))*coord_batch(st,qt);
			    
			}//q
		    }//s
		    
		}//ic
		
	    }//spherical

	}//q
	
    }//i

}//TwoC_DERI_contribution




void ::cchem::rimp2_gradient::DFIntGradient::gradient::
ThreeC_DERI_contribution(Array<double> * BARE_PN_,Array<double> * MNAUX_FORCE_){

    BOOST_AUTO(const &Ca, Ca_.get());
    ConstMapMatrixXd occ_coeff_mat(  Ca.data().begin(), Ca.size1(), Ca.size2() );


    //    if(pe_.rank() != 0) goto skip1;
    {

	utility::timer timerTest;
	timerTest.reset();
	size_t NMB = 1000;
	NMB=500;
		NMB=1000;
	//size_t NMB = 4000;
	//size_t NMB = 20;
	bool read_restrict = 1; bool ReadRow = 0; bool MagicSwap = 0;
	bool printDebug = 0; //(pe_.rank() == 0);
	bool loopRestrict = 0;
	::cchem::rimp2_gradient::detail::DFAsync 
	      makeMNAuxForce(NMB,
		      BARE_PN_, 1, 0, nl_, 0,
		      MNAUX_FORCE_,      1, 0, nl_, 0,
		      read_restrict, loopRestrict, ReadRow, MagicSwap, pe_, printDebug);
	
	::cchem::rimp2_gradient::DFIntGradient::gradient::MNAuxForceFuctor
	      Functor(
		      posFactor_[0].data(),negFactor_[0].data(),
		      N_,no_,Npi_,Nni_,
		      Ca.data().begin()
		      ) ;
	
	makeMNAuxForce.DoOpAsync_R_W_NODE_PARALLEL( Functor );
	
	if(pe_.rank() == 0)
	    std::cout << "     Time for MNAuxForceFuctor (CPU): " 
		      << timerTest << std::flush << std::endl; 
    }
    
    pe_.barrier();

    
    BOOST_AUTO(const &auxbasis, auxbasis_.get());
    BOOST_AUTO(const &auxshells, auxbasis.shells());
    BOOST_AUTO(const &basis, basis_.get());
    BOOST_AUTO(const &shells, basis.shells());
    
    BOOST_AUTO(const &Cv, Cv_.get());
    
    
    
    size_t NMB = 1000;
    NMB=500;
    NMB=1000;
    std::vector< std::pair<size_t,size_t> > AsShWork;
    std::vector< std::pair<size_t,size_t> > AsBasWork;
    size_t MaxBuffer = 0;
    this->setupAsyncMNAuxAccess(NMB, AsShWork, AsBasWork, MaxBuffer);
    
    // if(pe_.rank() ==0){
    // 	std::cout << "AsShWork.size() " << AsShWork.size() << std::endl;
    // 	std::cout << "AsBasWork.size() " << AsBasWork.size() << std::endl;
    // 	std::cout << "MaxBuffer " << MaxBuffer << std::endl;
    // }

    const size_t nMNAuxReads = AsBasWork.size();
    
    std::vector<size_t> startMNAuxRead(2);

    
    std::vector<size_t> finishMNAuxRead(2);

    //set up asynchronous object
    ::cchem::ri::async async;

    //asynchronous buffers
    double *pMNAux1 = new double[MaxBuffer*N_*N_];
    double *pMNAux2 = new double[MaxBuffer*N_*N_];

    Eigen::MatrixXd egNodeTemp(natoms_,3);
    egNodeTemp.setZero();


    if(pe_.rank() >= nMNAuxReads)goto nothingToDo;
    startMNAuxRead[0] = 0; // AsBasWork[0].first};
    startMNAuxRead[1] = AsBasWork[pe_.rank()].first;

    finishMNAuxRead[0] = N_*N_; // AsBasWork[0].second};
    finishMNAuxRead[1] = AsBasWork[pe_.rank()].second;


#pragma omp parallel
    {
	cchem::ri::AuxiliaryThreeCenterInt < ::rysq::ThreeCenterDerivativeEri > 
	    auxiliary_eri(
			  boost::cref(auxbasis),
			  boost::cref(basis),
			  boost::cref(Ca),
			  boost::cref(Cv) );
	
	Eigen::MatrixXd eg_temp(natoms_,3);
	eg_temp.setZero();
	
	const size_t MaxObsShell = basis.max().size();    // #functions in largest obs shell
	const size_t MaxAuxShell = auxbasis.max().size(); // #functions in largest aux shell
	double *pTransBatch = NULL;
	if(spherical_)pTransBatch = new double[MaxObsShell*MaxObsShell*MaxAuxShell];

#pragma omp single
	async.get(startMNAuxRead, finishMNAuxRead, *(MNAUX_FORCE_), pMNAux2);

	for(int read = pe_.rank(); read < nMNAuxReads; read += pe_.size() ){

#pragma omp single
	    {
		async.wait();
		std::swap(pMNAux1,pMNAux2);
		
		if(read + pe_.size() < nMNAuxReads){

		    startMNAuxRead[1] = AsBasWork[ read+pe_.size() ].first;
		    finishMNAuxRead[1] = AsBasWork[ read+pe_.size() ].second;

		    async.get(startMNAuxRead, finishMNAuxRead, *(MNAUX_FORCE_), pMNAux2);
		    
		}//( read < nMNAuxReads-1)
	    }//omp single
	    
	    const size_t nl_chunk = AsBasWork[read].second-AsBasWork[read].first;
	    MapMatrixXd MNAux(pMNAux1, N_*N_ , nl_chunk );
	    const size_t offset = AsBasWork[read].first;

#pragma omp for
	    for(size_t task = AsShWork[read].first; task < AsShWork[read].second; task++){
		const Basis::Shell &L = auxbasis.shells().at(task);
		const size_t LSize = L.size();
		const int latom = L.atom();

		for(size_t s = 0; s < basis.shells().size(); s++){
		    const Basis::Shell &S = basis.shells().at(s);
		    const size_t SSize = S.size();
		    const size_t SStart = S.start();
		    const int satom = S.atom();
		    
		    for (size_t q = 0; q <= s; ++q) {
			const Basis::Shell &Q = basis.shells().at(q);
			const size_t QSize = Q.size();
			const size_t SQSize = QSize*SSize;
			const size_t QStart = Q.start();
			const int qatom = Q.atom();
			
			if(latom == satom && latom == qatom)continue;
			
			double perm = 0.5;
			if( q!=s ) perm = 1.0;
			
			//when sorted: s >= q
			//MapMatrixXd threec_batch( auxiliary_eri(Q,S,L), Q.size()*S.size()*L.size(), 9);
			//double * ptr_deriv = threec_batch.data();
			const double * ptr_deriv = auxiliary_eri(Q,S,L);
			
			if(spherical_){
			    //if something is not working, do we have sorted shells???
			    
			    const size_t LSphSize = L.sphsize();
			    const size_t LSphStart = L.sphstart();
			    const double *pLmin = pSphTrans_[L.Lmin()];
			    
			    for(int nder = 0; nder < 9; nder++){
				int center = L.atom();
				if(nder < 6)center = Q.atom();
				if(nder < 3)center = S.atom();
				
				const int ic = (nder % 3);
				
				cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
					    SQSize, LSphSize, LSize,
					    1.0, ptr_deriv, SQSize,
					    pLmin, LSize,
					    0.0, pTransBatch, SQSize );

				double val3 = 0.0;
				for (int lt = 0; lt < LSphSize; lt++) {
				    const double *pInt = &pTransBatch[lt*SQSize];
				    
				    for(int st = SStart, qst=0; st < SStart+SSize; st++){
					for(int qt = QStart; qt < QStart+QSize; qt++,qst++){
					    

					    //coulomb
					    const double val = (PAux_[1](LSphStart+lt)*D2AO_[0]( st, qt) +
								DAux_[1](LSphStart+lt)*P2AO_[0]( st, qt) );

					    val3 += val*(pInt[qst]); //transformed_batch(qst ,lt );

					    //exchange
					    const double val2 = ( MNAux(st*N_ + qt , LSphStart + lt - offset) +
								  MNAux(qt*N_ + st , LSphStart + lt - offset) );
					    val3 -= val2*(pInt[qst]); //transformed_batch(qst ,lt );
					    
					}//qt
				    }//st
				    
				}//lt
				
				eg_temp(center,ic) += perm*val3;
				ptr_deriv += SQSize*LSize;

			    }//nder

			}else{ //(spherical_)


			    const size_t LStart = L.start();

			    for(int nder = 0; nder < 9; nder++){
				int center = L.atom();
				if(nder < 6)center = Q.atom();
				if(nder < 3)center = S.atom();

				const int ic = (nder % 3);

				double val3 = 0.0;
				for (int lt = 0; lt < LSize; lt++) {
				    const double *pInt = &ptr_deriv[lt*SQSize];

				    for(int st = SStart, qst=0; st < SStart+SSize; st++){
					for(int qt = QStart; qt < QStart+QSize; qt++,qst++){


					    //coulomb
					    const double val = (PAux_[1](LStart+lt)*D2AO_[0]( st, qt) +
								DAux_[1](LStart+lt)*P2AO_[0]( st, qt) );

					    val3 += val*(pInt[qst]); //deriv_batch(qst ,lt );

					    //exchange
					    const double val2 = (MNAux(st*N_ + qt , LStart + lt - offset) +
								  MNAux(qt*N_ + st , LStart + lt - offset) );
					    val3 -= val2*(pInt[qst]); //deriv_batch(qst ,lt );

					}//qt
				    }//st
				}//lt

				eg_temp(center,ic) += perm*val3;
				ptr_deriv += Q.size()*S.size()*L.size();

			    }//nder

			} //(spherical_)

		    }//q

		}//s

	    }//task

	}//read

#pragma omp critical
	egNodeTemp += eg_temp;

	delete [] pTransBatch;

    }//omp parallel

 nothingToDo:;

    pe_.reduce("+",egNodeTemp.data(), (size_t)(natoms_*3) );
    eg_ += egNodeTemp;
    

    delete [] pMNAux2;
    delete [] pMNAux1;

}//ThreeC_DERI_contribution
