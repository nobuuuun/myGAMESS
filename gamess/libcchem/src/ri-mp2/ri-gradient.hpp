/*
 * ri-gradient.hpp
 *
 *  Created on: Jun 16, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_MP2_RI_GRADIENT_HPP_
#define LIBCCHEM_SRC_MP2_RI_GRADIENT_HPP_


#include "runtime.hpp"
#include "core/wavefunction.hpp"

namespace cchem {
namespace rimp2_gradient {

    /*
      @brief      setup ri run
      
      @author     LBR
      
      @param      wf Wavefuntion
      @param      rt Runtime
      @param      molecule Molecule
      @param      EG energy gradient
    */
    double gradient(Wavefunction wf, Runtime &rt, Molecule molecule, double *EG);


}//namespace rimp2
}//namespace cchem


#endif /* LIBCCHEM_SRC_MP2_RI_GRADIENT_HPP_ */
