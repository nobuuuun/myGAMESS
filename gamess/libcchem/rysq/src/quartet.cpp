#include "foreach.hpp"
#include "rysq-core.hpp"
#include "ri-rysq-core.hpp"
#include "transpose.hpp"

using namespace rysq;

void Doublet<Shell>::initialize() {
    shells_[0] = &a_;
    shells_[1] = &b_;

    L_ =  detail::TwoC_L::twoc_add(*this);
    size_ = detail::TwoC_size::twoc_multiply(*this);
    K_ = detail::TwoC_K::twoc_multiply(*this);

    nc_ = 1;
    hybrid_ = 0;
    foreach (const Shell* shell, this->shells_) {
	nc_ *= shell->nc;
	hybrid_ += shell->is_hybrid();
    }
    // std::cout << *this << std::endl;
}

void Triplet<Shell>::initialize() {
    shells_[0] = &a_;
    shells_[1] = &b_;
    shells_[2] = &c_;

    L_ =  detail::ThreeC_L::threec_add(*this);
    size_ = detail::ThreeC_size::threec_multiply(*this);
    K_ = detail::ThreeC_K::threec_multiply(*this);

    nc_ = 1;
    hybrid_ = 0;
    foreach (const Shell* shell, this->shells_) {
	nc_ *= shell->nc;
	hybrid_ += shell->is_hybrid();
    }
    // std::cout << *this << std::endl;
}

void Quartet<Shell>::initialize() {
    shells_[0] = &a_;
    shells_[1] = &b_;
    shells_[2] = &c_;
    shells_[3] = &d_;

    L_ =  detail::L::add(*this);
    size_ = detail::size::multiply(*this);
    K_ = detail::K::multiply(*this);

    nc_ = 1;
    hybrid_ = 0;
    foreach (const Shell* shell, this->shells_) {
	nc_ *= shell->nc;
	hybrid_ += shell->is_hybrid();
    }
    // std::cout << *this << std::endl;
}


// Quartet<Shell> Quartet<Shell>::transposed(int transpose_mask) const {
//     UNPACK((const Shell *A, *B, *C, *D), this->shells_);

//     int shuffle_mask = Transpose::shuffle_mask(transpose_mask);
//     cxx::utility::permute(A, B, C, D, shuffle_mask);
//     // std::cout << shuffle_mask << "\n";
//     return Quartet<Shell>(const_cast<Shell&>(*A),
// 			  const_cast<Shell&>(*B),
// 			  const_cast<Shell&>(*C),
// 			  const_cast<Shell&>(*D));
// }
