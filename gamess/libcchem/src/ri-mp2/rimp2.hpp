/*
 * rimp2.hpp
 *
 *  Created on: May 7, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_SRC_MP2_RIMP2_HPP_
#define LIBCCHEM_SRC_MP2_RIMP2_HPP_


#include "runtime.hpp"
#include "core/wavefunction.hpp"

namespace cchem {
namespace rimp2 {

    double energy(Wavefunction wf, Runtime &rt);

}
} // namespace cchem


#endif /* LIBCCHEM_SRC_MP2_RIMP2_HPP_ */
