/*
 * twoc-forward.hpp
 *
 *  Created on: Apr 7, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_RI_FORWARD_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_RI_FORWARD_HPP_

namespace rysq {
	namespace kernel {
		namespace twoc_quadrature {

		template<class C>
		struct impl {
			static const bool value = false;
		};

		}
	}
}

namespace rysq {
	namespace kernel {
		namespace twoc_derivative_quadrature {

		template<class C>
		struct impl {
			static const bool value = false;
		};

		}
	}
}


namespace rysq {
	namespace kernel {
		namespace threec_derivative_quadrature {

		template<class C>
		struct impl {
			static const bool value = false;
		};

		}
	}
}



namespace rysq {
	namespace kernel {
		namespace twoc_derivative_overlap_quadrature {

		template<class C>
		struct impl {
			static const bool value = false;
		};

		}
	}
}



#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_RI_FORWARD_HPP_ */
