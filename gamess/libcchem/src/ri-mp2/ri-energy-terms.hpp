/*
 * ri-energy-terms.hpp
 *
 *  Created on: Jun 17, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_MP2_RI_ENERGY_TERMS_HPP_
#define LIBCCHEM_SRC_MP2_RI_ENERGY_TERMS_HPP_

namespace cchem{
namespace rimp2{
namespace detail {


using cchem::Thread;


struct Energy_Term : boost::noncopyable {

private:
	Thread thread_;
	size_t nd_,ns_,nv_,nl_,nvs_;
	int stride_;
	int stride_row_;
	double thread_e_term[7];

public:

	Energy_Term(){};

	Energy_Term(const size_t &nd, const size_t &ns, const size_t &nv,const size_t &nl,
			const int &stride,const int &stride_row)
	:nd_(nd),ns_(ns),nv_(nv),nl_(nl),nvs_(nv+ns),stride_(stride),stride_row_(stride_row)
	{
		for (int i = 0; i<7; i++){thread_e_term[i] = 0.0;}
	}

	~Energy_Term(){thread_.free();}

    //#define pos(a,b,offset) a + b*offset 

    template<class T, class V>
    double evaluate(const T &t, const double dij, const V &e) {
    // template<class V>
    // double evaluate(const double *t, const double dij, const V &e) {
	double E = 0;
	for (unsigned int b = ns_; b < nvs_; ++b) {
	    const double dijb = dij - e[b];
	    for (unsigned int a = ns_; a < nvs_; ++a) {
		const double ab = t(a,b);
		const double ba = t(b,a);
		// const double ab = t[ pos(a,b,nvs_) ];
		// const double ba = t[ pos(b,a,nvs_) ];
		const double de = dijb - e[a];
		E += ab*(2*ab - ba)/de;
	    }//a
	}//b
	return E;
    };//evaluate


    /**
    *   @brief singly unoccupied term
    *
    *   @details \sum_{\substack_ijsb}\frac{v_i_j^s^b(2v_i_j^s^b-v_j_i^s^b)}{\varepsilon_i+\varepsilon_j-\varepsilon_s^--\varepsilon_b}
    *
    *   @author Luke Roskop
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_m correction vector to singly occupied orbitals
    *   @param E energy component
    */
    template<class T, class V>
    double evaluate2(const T &t, const double dij, const V &e, std::vector<double> &e_m) {
	double E = 0;
	for (size_t s = 0; s < ns_; ++s) {
	    const double dijs = dij - e[s] - e_m[s];
	    for (size_t a = ns_; a < nvs_; ++a) {
		const double as = t(a,s);
		const double sa = t(s,a);
		//		const double de = dij - (e[a] + e[s]+ e_m[s]);
		const double de = dijs - e[a];
		E+= (as*as - as*sa + sa*sa)/de;
	    }//a
	}//s
	return E; // /(double)2;
    } //evaluate2
    /**
    *   @brief singly occupied
    *
    *   @author Luke Roskop
    *
    *   @details \sum_{\substack_sjab}\frac{v_s_j^a^b(2v_s_j^a^b-v_j_s^a^b)}{\varepsilon_s^++\varepsilon_j-\varepsilon_a-\varepsilon_b}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_m correction vector to singly occupied orbitals
    *   @param E energy component
    */
    template<class T, class V>
    double evaluate6(const T &t, const double dij, const V &e, const double e_m) {
	double E = 0;
	for (size_t b = ns_; b < nvs_; ++b) {
	    const double dijb = dij - e[b] - e_m;
	    for (size_t a = ns_; a < nvs_; ++a) {
		const double ab = t(a,b);
		const double ba = t(b,a);
		//		const double de = dij - (e[a] + e[b])-e_m;
		const double de = dijb - e[a];
		E += (ab*(2*ab - ba))/de;
	    }//a
	}//b
	return E/(double)2;
    } //evaluate6
    /**
    *   @brief singly unoccupied/occupied
    *
    *   @author Luke Roskop
    *
    *   @details \sum_{\substack_sjbt}\frac{(v_s_j^b^tv_s_j^b^t)}{\varepsilon_s^++\varepsilon_j-\varepsilon_b-\varepsilon_t^-}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_s correction vector to singly occupied orbitals
    *   @param e_a energy correction to orbital singly occupied orbital
    *   @param E energy component
    */
	template<class T, class V>
	double evaluate4(const T &t, const double dij, const V &e, const double e_a, std::vector<double> &e_s) {
	    double E = 0;
	   	    for (size_t s = 0; s < ns_; ++s) {
			const double dijs = dij - e[s] - e_a - e_s[s];			
	   		for (size_t a = ns_; a < nvs_; ++a) {
//        			double ab = t(s,a);  //this is fucky
			    const double ab = t(a,s);
			    //			    const double de = dij - (e[a] + e[s]) - e_a - e_s[s];
			    const double de = dijs - e[a];
			E += ab*ab/de;
	   		}//a
	    }//s
	    return E/(double)2;
	} //evaluate4


//		//only the master process (rank 0) continues on
//		//(so Fock-like energy is summed correctly)
////	    		BOOST_AUTO(const &ed, wf.e(wf.double_occ()));
	template<class V, class U>
	double evaluate5(const V &ea, const V &ev, U &pfock){
	    double E = 0;
	    for (size_t d = 0; d < nd_; d++){
		double *p = pfock.at(d);
		for (size_t a = 0; a < nv_; a++){
		    E += pow((*p),2)/(ea[d]-ev[a+ns_]);
		    p++;
		}//a
	    }//d
	    return E/(double)2;
	} //evaluate5


	template<class V, class U>
	double evaluate5p5(const V &ea, const V &ev, U &fock){
	    double E = 0;
	    for (size_t d = 0; d < nd_; d++){
		//		double *p = pfock.at(d);
		for (size_t a = 0; a < nv_; a++){
		    //		    double p = fock(d,a);
		    const double p = fock(a,d);
		    E += pow(p,2)/(ea[d]-ev[a+ns_]);
		    //	p++;
		}//a
	    }//d
	    return E/(double)2;
	} //evaluate5p5

	/**
    *   @brief 2 singly occupied
    *
    *   @author Luke Roskop
    *
    *   @details \frac{1}{4} \sum_{\substack_stab}\frac{(v_s_t^a^b-v_t_s^a^b)(v_s_t^a^b-v_t_s^a^b)}{\varepsilon_s^++\varepsilon_t^+-\varepsilon_a-\varepsilon_b}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_s energy correction to orbital singly occupied orbital
    *   @param e_t energy correction to orbital singly occupied orbital
    *   @param E energy component
    */
    template<class T, class V>
    double evaluate3(const T &t, const double dij, const V &e, const double e_s, const double e_t) {
	double E = 0;
	for (size_t b = ns_; b < nvs_; ++b) {
	    const double dijb = dij - e[b] - e_s - e_t;
	    for (size_t a = ns_; a < nvs_; ++a) {
		const double ab = t(a,b);
		const double ba = t(b,a);
		//		const double de =dij -(e[a] + e[b]) - e_s - e_t;
		const double de = dijb - e[a];
		E += ((ab-ba)*(ab-ba))/de;
	    }//a
	}//b
	return E/(double)4;
    } //evaluate3
    
    /**
    *   @brief doubly unoccupied
    *
    *   @author Luke Roskop
    *
    *   @details \frac{1}{4} \sum_{\substack_ijst}\frac{(v_i_j^s^t-v_j_i^s^t)(v_i_j^s^t-v_j_i^s^t)}{\varepsilon_i+\varepsilon_j-\varepsilon_s^--\varepsilon_t^-}
    *
    *   @param t integrals
    *   @param dij half formed energy denominator
    *   @param e singly occupied + virtual orbital energies
    *   @param ns number of singly occupied orbitals
    *   @param e_m correction vector to singly occupied orbitals
    *   @param E energy component
    */
	template<class T, class V>
	double evaluate7(const T &t, const double dij, const V &e, std::vector<double> &e_m) {
	    double E = 0;
	    for (size_t b = 0; b < ns_; ++b) {
		const double dijb =dij - e[b] - e_m[b];
	    	for (size_t a = 0; a < ns_; ++a) {
	    		const double ab = t(a,b);
	    		const double ba = t(b,a);
			//	    		const double de =dij -(e[a] + e[b]) -e_m[a] -e_m[b];
	    		const double de =dijb -e[a] -e_m[a];
	    		E += ((ab-ba)*(ab-ba))/de;
	    	}//a
	    }//b
	    return E/(double)4;
	} //evaluate7

	template<class T, class V>
	double openshell_ii(const T &eri, const double dij, const V &ev,
			std::vector<double> &e_m, const int &iocc){

		if(iocc < nd_ ){
		    thread_e_term[1] += this->evaluate2(eri,dij,ev,e_m);
		    thread_e_term[6] += this->evaluate7(eri, dij, ev, e_m);
		}//)iocc < nd_)

		if(iocc >= nd_){
		    thread_e_term[2] += this->evaluate3(eri, dij, ev, e_m[iocc-nd_], e_m[iocc-nd_]);
		}//(iocc >= nd_)

	}//energy_ii


	template<class T, class V> //first
	double openshell_ij(const T &eri, const double dij, const V &ev,
			std::vector<double> &e_m, const int &iocc, const int &jocc,
			std::vector<double *> &pfock){

		if(iocc < nd_ ){
		    thread_e_term[1] += 2*this->evaluate2(eri,dij,ev,e_m);
		    thread_e_term[6] += 2*this->evaluate7(eri, dij, ev, e_m);
		}//(iocc < nd_)

		//singly occupied E6 / singly occupied/unoccupied E4
		if(iocc >= nd_ && jocc < nd_)	    {
		    thread_e_term[5] += 2*this->evaluate6(eri, dij, ev,  e_m[iocc -nd_]);
		    thread_e_term[3] += 2*this->evaluate4(eri, dij, ev,  e_m[iocc -nd_], e_m);
		    double *p = pfock.at(jocc);
//		    for (int a = ns_; a < eri.cols(); a++){*p += eri(a,iocc-nd_);	p++;}
		    for (int a = ns_,ipos = 0; a < nvs_; a++,ipos++){
			p[ipos] += eri(a,iocc-nd_);}
		} //(iocc >= nd_ && jocc < nd_)

		if(jocc >= nd_){
		    thread_e_term[2] += 2*this->evaluate3(eri, dij, ev, e_m[iocc-nd_], e_m[jocc-nd_]);
		} // (jocc >= nd_)

	}//energy_ij


	template<class T, class U, class V> //second IJ blocks
	void operator()(const T &eri, const int &ipos,
			const int &jblock, const int &jblock_end, const int &iocc,
			const U &ea, const V &ev, std::vector<double> &e_m, std::vector<double *> &pfock){

		for(int jocc = jblock; jocc < jblock_end; jocc++){

			const double dij = ea[iocc] + ea[jocc];
			
			if(iocc < nd_)thread_e_term[0] += 2*this->evaluate(eri.block(ipos*nvs_,(jocc-jblock)*nvs_,nvs_,nvs_),dij,ev);
// const double *ptrTest = &eri.data()[(jocc-jblock)*nvs_];
// if(iocc < nd_)thread_e_term[0] += ( 1+ (iocc!=jocc) )*this->evaluate(ptrTest,dij,ev);

			if(ns_ > 0)this->openshell_ij(eri.block(ipos*nvs_,(jocc-jblock)*nvs_,nvs_,nvs_),dij,ev,e_m,iocc,jocc,pfock);

		} //jocc

	}//operator()


	template<class T, class U, class V> //first  II blocks
	void operator()(const T &eri, const int &jblock, const int &iocc,
			const U &ea, const V &ev, std::vector<double> &e_m, std::vector<double *> &pfock){

		for(int jocc = jblock; jocc <= iocc; jocc++){

			const double dij = ea[iocc] + ea[jocc];

			if(ns_ > 0){
				if(iocc == jocc)this->openshell_ii(eri.block(0,(jocc-jblock)*nvs_,nvs_,nvs_),dij,ev,e_m,iocc);
				if(iocc != jocc)this->openshell_ij(eri.block(0,(jocc-jblock)*nvs_,nvs_,nvs_),dij,ev,e_m,iocc,jocc,pfock);
			}//ns
			if(iocc < nd_ && jocc < nd_)thread_e_term[0] += ( 1+ (iocc!=jocc) )*this->evaluate(eri.block(0,(jocc-jblock)*nvs_,nvs_,nvs_),dij,ev);
						
// const double *ptrTest = &eri.data()[(jocc-jblock)*nvs_];
// if(iocc < nd_ && jocc < nd_)thread_e_term[0] += ( 1+ (iocc!=jocc) )*this->evaluate(ptrTest,dij,ev);
		}//jocc

	}//operator()



//	void operator()(double (&e_term)[7]){
		void operator()(double *e_term){
#pragma omp critical // accumulate ZAPT energy components
		for (int i = 0; i < 7; i++) e_term[i] += thread_e_term[i];
#pragma omp barrier  //make sure threads are done contributing to pfock[i]->vector

	};//operator()

}; // Energy_Term

}//namespace detail
}//namespace rimp2
}//namespace cchem


#endif /* LIBCCHEM_SRC_MP2_RI_ENERGY_TERMS_HPP_ */
