/*
 * libcchem_libint.h
 *
 *  Created on: Oct 29, 2015
 *      Author: luke
 */

#ifndef SRC_LIBCCHEM_LIBINT_H_
#define SRC_LIBCCHEM_LIBINT_H_

//#include <iostream>

//#include "/Users/luke/devel/install/libint_deri/include/libint2.h"
//#include "/Users/luke/devel/install/libint_deri/include/libint2/intrinsic_types.h"
//#include <libint2.hpp>

//#include "core/wavefunction.hpp"



namespace libcchem_libint_interface {


void T_V_1body_deriv_contributions(size_t nmo, double* pmat_data, ::basis::Basis basis, Molecule molecule,
			double * ptr_wmat, Eigen::MatrixXd &TTerm, Eigen::MatrixXd &VTerm);


void twobody_deriv_contributions(size_t nmo, double* pmat_pmp2, double * ptr_pscf, ::basis::Basis basis, Molecule molecule,double *ptr_eg_global,
		Eigen::MatrixXd &schwartz_mat);

void twobody_eri(size_t nmo, ::basis::Basis basis, Molecule molecule,
		Eigen::MatrixXd &schwartz_mat);

void compute_schwartz_ints(::basis::Basis basis, Eigen::MatrixXd &schwartz_mat);

}//namespace libcchem_libint_interface


#endif /* SRC_LIBCCHEM_LIBINT_H_ */
