///*
// * ri-async.cpp
// *
// *  Created on: Feb 5, 2016
// *      Author: luke
// */
//


#include "god.hpp"

#include "thread.hpp"
#include <Eigen/Dense>
#include <ri-async.hpp>

void ::cchem::rimp2_gradient::detail::DFAsync:: 
set_up_array_info(bool &ReadRow, bool &async, size_t &NWorkingBytes,
		  ArrayData &A1,  size_t NA1RowBStart, size_t NA1ColBStart){

    size_t NBytes = NWorkingBytes;

    A1.NABRows_ = A1.Array_->shape()[0]/(A1.NAB_    + NA1RowBStart); //is this case, A1.NAB_ had NAB1S subtracted out earlier
    A1.NABCols_ = A1.Array_->shape()[1]/(A1.NAColB_ + NA1ColBStart);

    size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);

    //A1.arrayType_ = A1.Array_->str();

    if(ReadRow){

	//we'll need at least one block from A1
	A1.MaxNAB_=1;
	NBytes -= (1+async)*(A1BlockSize);
	//determine A1.MaxNAB_: how many row blocks that can fit into buffer (size NMB)
	while(NBytes > (1+async)*A1BlockSize ){
	    A1.MaxNAB_++;
	    NBytes -= (1+async)*A1BlockSize;
	    if(A1.MaxNAB_ == A1.NAB_)break;
	}

	//figure out how many reads will need to occur
	A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);

    }else{

	//we'll need at least one block from A1
	A1.MaxNAB_=1;
	NBytes -= (1+async)*(A1BlockSize);

	//determine A1.MaxNAB_: how many column blocks that can fit into buffer (size NMB)
	while(NBytes > (1+async)*A1BlockSize ){
	    A1.MaxNAB_++;
	    NBytes -= (1+async)*A1BlockSize;
	    if(A1.MaxNAB_ == A1.NAColB_)break;
	}

	//figure out how many reads will need to occur
	A1.NAReads_ = A1.NAColB_/A1.MaxNAB_ + (A1.NAColB_%A1.MaxNAB_>0);

    }//(ReadRow)

    //set up array access vectors
    this->set_up_block_info(A1, NA1RowBStart, NA1ColBStart, ReadRow);

}//set_up_single_array_info


void ::cchem::rimp2_gradient::detail::DFAsync:: 
set_up_pair_array_info(bool &restrict, bool &ReadRow , bool &async, size_t & NWorkingBytes,
		       ArrayData &A1,  size_t NA1RowBStart,size_t NA1ColBStart,
		       ArrayData &A2,  size_t NA2RowBStart,size_t NA2ColBStart){



    //we need to see how many blocks of A1 and A2 can
    // be loaded into the buffers at the same time
    //		size_t NBytes = NMB_*1024*1024;
    size_t NBytes = NWorkingBytes;
    //		std::cout << NWorkingBytes << std::endl;
    A1.NABRows_ = A1.Array_->shape()[0]/(A1.NAB_    + NA1RowBStart); //is this case, A1.NAB_ had NAB1S subtracted out earlier
    A1.NABCols_ = A1.Array_->shape()[1]/(A1.NAColB_ + NA1ColBStart);
    A2.NABRows_ = A2.Array_->shape()[0]/(A2.NAB_    + NA2RowBStart);
    A2.NABCols_ = A2.Array_->shape()[1]/(A2.NAColB_ + NA2ColBStart);
    //				std::cout << A1.Array_->shape()[0] << " " << A1.NAB_   << " " << NA1RowBStart << std::endl;
    //				std::cout << A1.NABRows_ << " " << A1.NABCols_ << " " << A2.NABRows_ << " " << A2.NABCols_ << std::endl;
    size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);
    size_t A2BlockSize = A2.NABRows_*A2.NABCols_*sizeof(double);

    //		std::cout << A1.NABRows_ << " " << A1.NABRows_ << std::endl;
    //		std::cout << A1.NAB_ << " " << std::endl;
    //		std::cout << "xx " <<  A1.Array_->shape()[0] << " " << A1.NAB_  << " " << NA1RowBStart << std::endl;
    //		std::cout << A2.NABRows_ << " " << A2.NABRows_ << std::endl;
    if(restrict){

	if(ReadRow){
	    //we'll need at least one block from A1 and A2
	    A1.MaxNAB_=1;

	    NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
	    //determine A1.MaxNAB_ (=A2.MaxNAB_):
	    //     how many row blocks that can fit into buffer (size NMB)
	    //				std::cout << "gates "<< NBytes << " " << (1+async)*(A1BlockSize + A2BlockSize) << std::endl;
	    while(NBytes > (1+async)*(A1BlockSize + A2BlockSize) ){
		//					std::cout << "made it" << std::endl;
		A1.MaxNAB_++;
		NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
		if(A1.MaxNAB_ == A1.NAB_)break;
	    }
	    A2.MaxNAB_ = A1.MaxNAB_;

	    //figure out how many reads will need to occur
	    A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
	    A2.NAReads_ = A1.NAReads_;

	    //try to balance out read size
	    A1.MaxNAB_ = A1.NAB_/A1.NAReads_ + (A1.NAB_%A1.NAReads_>0);
	    A2.MaxNAB_ = A2.NAB_/A2.NAReads_ + (A2.NAB_%A2.NAReads_>0);

	}else{

	    //we'll need at least one block from A1 and A2
	    A1.MaxNAB_=1;
	    NBytes -= (1+async)*(A1BlockSize + A2BlockSize);

	    //determine A1.MaxNAB_ (=A2.MaxNAB_):
	    //     how many col blocks that can fit into buffer (size NMB)
	    while(NBytes > (1+async)*(A1BlockSize + A2BlockSize) ){
		A1.MaxNAB_++;
		NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
		if(A1.MaxNAB_ == A1.NAColB_)break;
	    }

	    A2.MaxNAB_ = A1.MaxNAB_;

	    //figure out how many reads will need to occur
	    A1.NAReads_ = A1.NAColB_/A1.MaxNAB_ + (A1.NAColB_%A1.MaxNAB_>0);
	    A2.NAReads_ = A1.NAReads_;

	    //try to balance out read size
	    A1.MaxNAB_ = A1.NAColB_/A1.NAReads_ + (A1.NAColB_%A1.NAReads_>0);
	    A2.MaxNAB_ = A2.NAColB_/A2.NAReads_ + (A2.NAColB_%A2.NAReads_>0);

	}//(ReadRow)

	//			std::cout << " about to setup array access vectors" << std::flush << std::endl;
	//set up array access vectors
	this->set_up_block_info(A1,  NA1RowBStart, NA1ColBStart, ReadRow);
	this->set_up_block_info(A2,  NA2RowBStart, NA2ColBStart, ReadRow);

    }else{ //(restrict)



	if(ReadRow){

	    //we'll need at least one block from A1 and A2
	    A1.MaxNAB_=1;
	    A2.MaxNAB_=1;

	    NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
	    //determine A1.MaxNAB_: how many row blocks that can fit into buffer (size NMB)
	    while(NBytes > (1+async)*A1BlockSize ){
		A1.MaxNAB_++;
		NBytes -= (1+async)*A1BlockSize;
		if(A1.MaxNAB_ == A1.NAB_)break;
	    }

	    //determine A2.MaxNAB_: how many row blocks that can fit into buffer (size NMB)
	    while(NBytes > (1+async)*A2BlockSize ){
		A2.MaxNAB_++;
		NBytes -= (1+async)*A2BlockSize;
		if(A2.MaxNAB_ == A2.NAB_)break;
	    }

	    //figure out how many reads will need to occur
	    A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
	    A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);

	    //try to balance out read size
	    A1.MaxNAB_ = A1.NAB_/A1.NAReads_ + (A1.NAB_%A1.NAReads_>0);
	    A2.MaxNAB_ = A2.NAB_/A2.NAReads_ + (A2.NAB_%A2.NAReads_>0);

	}else{

	    //we'll need at least one block from A1 and A2
	    A1.MaxNAB_=1;
	    A2.MaxNAB_=1;

	    NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
	    //determine A1.MaxNAB_: how many col blocks that can fit into buffer (size NMB)
	    while(NBytes > (1+async)*A1BlockSize ){
		A1.MaxNAB_++;
		NBytes -= (1+async)*A1BlockSize;
		if(A1.MaxNAB_ == A1.NAColB_)break;
	    }

	    //determine A2.MaxNAB_: how many col blocks that can fit into buffer (size NMB)
	    while(NBytes > (1+async)*A2BlockSize ){
		A2.MaxNAB_++;
		NBytes -= (1+async)*A2BlockSize;
		if(A2.MaxNAB_ == A2.NAColB_)break;
	    }

	    //figure out how many reads will need to occur
	    A1.NAReads_ = A1.NAColB_/A1.MaxNAB_ + (A1.NAColB_%A1.MaxNAB_>0);
	    A2.NAReads_ = A2.NAColB_/A2.MaxNAB_ + (A2.NAColB_%A2.MaxNAB_>0);

	    //				std::cout << A1.NAReads_ << std::endl;
	    //try to balance out read size
	    A1.MaxNAB_ = A1.NAColB_/A1.NAReads_ + (A1.NAColB_%A1.NAReads_>0);
	    A2.MaxNAB_ = A2.NAColB_/A2.NAReads_ + (A2.NAColB_%A2.NAReads_>0);

	}//ReadRow


	//set up array access vectors
	this->set_up_block_info(A1,  NA1RowBStart, NA1ColBStart, ReadRow);
	this->set_up_block_info(A2,  NA2RowBStart, NA2ColBStart, ReadRow);

    }//(restrict)

}//set_up_array_info


void ::cchem::rimp2_gradient::detail::DFAsync:: 
set_up_block_info(ArrayData &A,  size_t NARowBStart, size_t NAColBStart, bool ReadRow){

    if(A.NAReads_ == 1){

	std::pair< std::vector<size_t>, std::vector<size_t> > temp;
	temp.first.push_back(NARowBStart*A.NABRows_);
	temp.first.push_back(NAColBStart*A.NABCols_);
	temp.second.push_back( (NARowBStart+A.NAB_)*A.NABRows_);
	temp.second.push_back( (NAColBStart+A.NAColB_)*A.NABCols_);
	//			std::cout << NAColBStart << " " << A.NAColB_ << " " <<A.NABCols_ << std::endl;
	A.AReadPos_.push_back(temp);
	if(ReadRow)A.ABlocks_.push_back( std::pair<size_t,size_t> ( NARowBStart,NARowBStart+A.NAB_) );
	if(!ReadRow)A.ABlocks_.push_back( std::pair<size_t,size_t> ( NAColBStart,NAColBStart+A.NAColB_) );

    }else{ //A.NAReads_ == 1

	if(ReadRow){

	    for(int i =0; i<A.NAReads_-1; i++){

		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
		temp.first.push_back( (NARowBStart+i*A.MaxNAB_)*A.NABRows_);
		temp.first.push_back( 0);
		temp.second.push_back( (NARowBStart+(i+1)*A.MaxNAB_)*A.NABRows_);
		temp.second.push_back(A.NABCols_);
		A.AReadPos_.push_back(temp);
		A.ABlocks_.push_back( std::pair<size_t,size_t> ( NARowBStart+i*A.MaxNAB_ , NARowBStart+(i+1)*A.MaxNAB_) );

	    }//i

	    std::pair< std::vector<size_t>, std::vector<size_t> > temp;
	    temp.first.push_back( (NARowBStart+(A.NAReads_-1)*A.MaxNAB_) *A.NABRows_);
	    temp.first.push_back(0);
	    temp.second.push_back( (NARowBStart+A.NAB_)*A.NABRows_ );
	    temp.second.push_back(A.NABCols_);
	    A.AReadPos_.push_back(temp);
	    A.ABlocks_.push_back( std::pair<size_t,size_t> ( NARowBStart+(A.NAReads_-1)*A.MaxNAB_ , NARowBStart+A.NAB_) );

	}else{  //ReadRow

	    for(int i =0; i<A.NAReads_-1; i++){

		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
		temp.first.push_back( 0 );
		temp.first.push_back( (NAColBStart+i*A.MaxNAB_)*A.NABCols_);
		temp.second.push_back( A.NABRows_ );
		temp.second.push_back( (NAColBStart+(i+1)*A.MaxNAB_)*A.NABCols_ );
		A.AReadPos_.push_back(temp);
		A.ABlocks_.push_back( std::pair<size_t,size_t> ( NAColBStart+i*A.MaxNAB_ , NAColBStart+(i+1)*A.MaxNAB_) );

	    }//i

	    std::pair< std::vector<size_t>, std::vector<size_t> > temp;
	    temp.first.push_back( 0 );
	    temp.first.push_back( (NAColBStart+(A.NAReads_-1)*A.MaxNAB_) *A.NABCols_ );
	    temp.second.push_back( A.NABRows_ );
	    temp.second.push_back( (NAColBStart+A.NAColB_)*A.NABCols_  );
	    A.AReadPos_.push_back(temp);
	    A.ABlocks_.push_back( std::pair<size_t,size_t> ( NAColBStart+(A.NAReads_-1)*A.MaxNAB_ , NAColBStart+A.NAColB_) );

	}//ReadRow

    }//(A.NAReads_ == 1)

}//set_up_block_info







void ::cchem::rimp2_gradient::detail::DFAsync::
PrintArray(int i){

    ArrayData &A = Array_[i];
    std::cout << std::endl << "PRINTING INFORMATION FOR Array_[" << i << "]:" << std::endl ;
    std::cout << "  full array "<< i <<" looks like " << A.Array_->shape()[0] << " x " << A.Array_->shape()[1] << std::endl;
    std::cout << "  smallest array block size: " << A.NABRows_  << " x " <<  A.NABCols_  << std::endl;
    std::cout << "  need to perform:           " << A.NAReads_ << " NAReads_" << std::endl;
    std::cout << "  Largest Block:             " << A.MaxNAB_ << std::endl;
    std::cout << "  number of row blocks:      " << A.NAB_  << std::endl;
    std::cout << "  number of column blocks:   " << A.NAColB_  << std::endl;
    std::cout << "  Rows/Block:                " << A.NABRows_  << std::endl;
    std::cout << "  Cols/Block:                " << A.NABCols_  << std::endl;

    std::cout << "  array access looks like: " << std::endl;
    for(int i = 0; i < A.NAReads_; i++){
	std::cout << "    " <<"["<< i <<"]"<< "   { " << A.AReadPos_[i].first[0] << " " << A.AReadPos_[i].first[1] << " }"
		  <<" --> { "<< A.AReadPos_[i].second[0] << " " << A.AReadPos_[i].second[1] << " }" << std::endl;
    }
    std::cout << "  array blocks range like: " << std::endl;
    for(int i = 0; i < A.NAReads_; i++){
	std::cout << "    " <<"["<< i  <<"]"<< "   " << A.ABlocks_[i].first << " : " << A.ABlocks_[i].second << std::endl;
    }
    std::cout << std::endl;
}//Print







// void ::cchem::rimp2_gradient::detail::RIAsync::
// set_up_array_info(bool &restrict, bool &async, size_t & NWorkingBytes,
// 		  ArrayData &A1,  size_t NAB1S, ArrayData &A2,  size_t NAB2S){

// 	//		ArrayData &A1,  ArrayData &A2){

// 	//we need to see how many blocks of A1 and A2 can
// 	// be loaded into the buffers at the same time
// 	//		size_t NBytes = NMB_*1024*1024;
// 	size_t NBytes = NWorkingBytes;
// //luke
// 	A1.NABRows_ = A1.Array_->shape()[0]/(A1.NAB_+NAB1S); //is this case, A1.NAB_ had NAB1S subtracted out earlier
// 	A1.NABCols_ = A1.Array_->shape()[1];
// 	A2.NABRows_ = A2.Array_->shape()[0]/(A2.NAB_+NAB2S);
// 	A2.NABCols_ = A2.Array_->shape()[1];

// 	size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);
// 	size_t A2BlockSize = A2.NABRows_*A2.NABCols_*sizeof(double);

// 	if(restrict){

// 		//we'll need at least one block from A1 and A2
// 		A1.MaxNAB_=1;

// 		NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
// 		//see how many CijQ/BARE_IA blocks we can fit into the buffer
// 		while(NBytes > (1+async)*(A1BlockSize + A2BlockSize) ){
// 			A1.MaxNAB_++;
// 			NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
// 			if(A1.MaxNAB_ == A1.NAB_)break;
// 		}
// 		A2.MaxNAB_ = A1.MaxNAB_;

// 		//figure out how many reads will need to occur
// 		A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 		A2.NAReads_ = A1.NAReads_;

// 		this->set_up_block_info(A1,  NAB1S);

// 		//restrict: A1 to be the same as A2
// 		A2.AReadPos_ = A1.AReadPos_;
// 		A2.ABlocks_ = A1.ABlocks_;

// 		//			A1.restrict_ = 1;
// 		//			A2.restrict_ = 1;

// 	}else{//(restrict)

// 		//we'll need at least one block from A1 and A2
// 		A1.MaxNAB_=1;
// 		A2.MaxNAB_=1;

// 		NBytes -= (1+async)*(A1BlockSize + A2BlockSize);
// 		//see how many A1.MaxNAB_ blocks we can fit into the buffer
// 		while(NBytes > (1+async)*A1BlockSize ){
// 			A1.MaxNAB_++;
// 			NBytes -= (1+async)*A1BlockSize;
// 			if(A1.MaxNAB_ == A1.NAB_)break;
// 		}

// 		//see how many A2.MaxNAB_ blocks we can fit into the buffer
// 		while(NBytes > (1+async)*A2BlockSize ){
// 			A2.MaxNAB_++;
// 			NBytes -= (1+async)*A2BlockSize;
// 			if(A2.MaxNAB_ == A2.NAB_)break;
// 		}

// 		//figure out how many reads will need to occur
// 		A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 		A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);

// 		//try to balance out read size
// 		A1.MaxNAB_ = A1.NAB_/A1.NAReads_ + (A1.NAB_%A1.NAReads_>0);
// 		A2.MaxNAB_ = A2.NAB_/A2.NAReads_ + (A2.NAB_%A2.NAReads_>0);

// 		this->set_up_block_info(A1, NAB1S);
// 		this->set_up_block_info(A2, NAB2S);


// 	}//(restrict)


// }//set_up_array_info






// void ::cchem::rimp2_gradient::detail::RIAsync::set_up_single_array_info(bool &row,bool &async, size_t & NWorkingBytes,
// 		ArrayData &A1,  size_t NAB1S){

// 	size_t NBytes = NWorkingBytes;

// 	A1.NABRows_ = A1.Array_->shape()[0]; //number of rows
// 	A1.NABCols_ = A1.Array_->shape()[1]/A1.NAB_; // min number of columns per block (whole column)

// 	size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);


// 	//we'll need at least one block from A1
// 	A1.MaxNAB_=1;

// 	NBytes -= (1+async)*(A1BlockSize);
// 	//see how many A1.MaxNAB_ blocks we can fit into the buffer
// 	while(NBytes > (1+async)*A1BlockSize ){

// 		A1.MaxNAB_++;
// 		NBytes -= (1+async)*A1BlockSize;
// 		if(A1.MaxNAB_ == A1.NAB_)break;
// 	}



// 	//figure out how many reads will need to occur
// 	A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);


// 	this->set_up_block_info_cols(A1, NAB1S);




// }//set_up_single_array_info




// void ::cchem::rimp2_gradient::detail::RIAsync::set_up_two_array_info(bool &row,bool &async, size_t & NWorkingBytes,
// 		ArrayData &A1,  size_t NAB1S,
// 		ArrayData &A2,  size_t NAB2S){

// 	size_t NBytes = NWorkingBytes;

// 	A1.NABRows_ = A1.Array_->shape()[0]; //number of rows
// 	A1.NABCols_ = A1.Array_->shape()[1]/A1.NAB_; // min number of columns per block (whole column)

// 	A2.NABRows_ = A2.Array_->shape()[0]; //number of rows
// 	A2.NABCols_ = A2.Array_->shape()[1]/A2.NAB_; // min number of columns per block (whole column)


// 	size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);
// 	size_t A2BlockSize = A2.NABRows_*A2.NABCols_*sizeof(double);


// 	//we'll need at least one block from A1
// 	A1.MaxNAB_=1;
// 	A2.MaxNAB_=1;

// 	NBytes -= (1+async)*(A1BlockSize+A2BlockSize);
// 	//see how many A1.MaxNAB_ blocks we can fit into the buffer
// 	while(NBytes > (1+async)*(A1BlockSize+A2BlockSize) ){

// 		A1.MaxNAB_++;
// 		A2.MaxNAB_++;
// 		NBytes -= (1+async)*(A1BlockSize+A2BlockSize);
// 		if(A1.MaxNAB_ == A1.NAB_)break;
// 	}



// 	//figure out how many reads will need to occur
// 	A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 	A2.NAReads_ = A1.NAReads_;

// 	this->set_up_block_info_cols(A1, NAB1S);
// 	this->set_up_block_info_cols(A2, NAB2S);




// }//set_up_two_array_info




// void ::cchem::rimp2_gradient::detail::RIAsync::set_up_unrestricted_array_info(bool &async, size_t &NWorkingBytes,
// 		ArrayData &A1, ArrayData &A2) {

// 	//we need to see how many blocks of A1 and A2 can
// 	// be loaded into the buffers at the same time
// 	//		size_t NBytes = NMB_*1024*1024;
// 	size_t NBytes = NWorkingBytes;

// 	A2.NABRows_ = A2.Array_->shape()[0]/A2.NAB_;
// 	A2.NABCols_ = A2.Array_->shape()[1];

// 	size_t A1BlockSize = A1.NABRows_*A1.NABCols_*sizeof(double);
// 	size_t A2BlockSize = A2.NABRows_*A2.NABCols_*sizeof(double);

// 	//		if(restrict){

// 	//we'll need at least one block from A2 (A1 is already known)
// 	A2.MaxNAB_=1;

// 	NBytes -= A1BlockSize*A1.MaxNAB_;
// 	NBytes -= (A2BlockSize);
// 	//see how many A1,A2,A3 blocks we can fit into the buffer
// 	while(NBytes > (A2BlockSize) ){
// 		A2.MaxNAB_++;
// 		NBytes -= (A2BlockSize);
// 		if(A2.MaxNAB_ == A2.NAB_)break;
// 	}

// 	//figure out how many reads will need to occur
// 	A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);

// 	this->set_up_block_info(A2, 0);

// 	//this just means that the buffers are restrict to NMB
// 	//			A2.restrict_ = 0;


// 	//		}else{//(restrict)
// 	//
// 	//			//we'll need at least one block from A1 and A2
// 	//			A1.MaxNAB_=1;
// 	//			A2.MaxNAB_=1;
// 	//
// 	//			NBytes -= (A1BlockSize + A2BlockSize);
// 	//			//see how many A1.MaxNAB_ blocks we can fit into the buffer
// 	//			while(NBytes > A1BlockSize ){
// 	//				A1.MaxNAB_++;
// 	//				NBytes -= A1BlockSize;
// 	//				if(A1.MaxNAB_ == A1.NAB_)break;
// 	//			}
// 	//
// 	//			//see how many A2.MaxNAB_ blocks we can fit into the buffer
// 	//			while(NBytes > A2BlockSize ){
// 	//				A2.MaxNAB_++;
// 	//				NBytes -= A2BlockSize;
// 	//				if(A2.MaxNAB_ == A2.NAB_)break;
// 	//			}
// 	//
// 	//			//figure out how many reads will need to occur
// 	//			A1.NAReads_ = A1.NAB_/A1.MaxNAB_ + (A1.NAB_%A1.MaxNAB_>0);
// 	//			A2.NAReads_ = A2.NAB_/A2.MaxNAB_ + (A2.NAB_%A2.MaxNAB_>0);
// 	//
// 	//			this->set_up_block_info(A1);
// 	//			this->set_up_block_info(A2);
// 	//
// 	//			A1.restrict_ = 0;
// 	//			A2.restrict_ = 0;
// 	//
// 	//		}//(restrict)

// }//set_up_3array_info


// void ::cchem::rimp2_gradient::detail::RIAsync::set_up_block_info(ArrayData &A1,  size_t NAB1S){

// 	if(A1.NAReads_ == 1){
// 		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 		temp.first.push_back(NAB1S*A1.NABRows_); temp.first.push_back(0);
// 		temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABRows_); temp.second.push_back(A1.NABCols_);
// 		A1.AReadPos_.push_back(temp);
// 		A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S,NAB1S+A1.NAB_) );
// 	}else{
// 		for(int i =0; i<A1.NAReads_-1; i++){
// 			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 			temp.first.push_back( (NAB1S+i*A1.MaxNAB_)*A1.NABRows_); temp.first.push_back(0);
// 			temp.second.push_back( (NAB1S+(i+1)*A1.MaxNAB_)*A1.NABRows_); temp.second.push_back(A1.NABCols_);
// 			A1.AReadPos_.push_back(temp);
// 			A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+i*A1.MaxNAB_ , NAB1S+(i+1)*A1.MaxNAB_) );
// 		}//i
// 		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 		temp.first.push_back( (NAB1S+(A1.NAReads_-1)*A1.MaxNAB_) *A1.NABRows_); temp.first.push_back(0);
// 		temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABRows_ ); temp.second.push_back(A1.NABCols_);
// 		A1.AReadPos_.push_back(temp);
// 		A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+(A1.NAReads_-1)*A1.MaxNAB_ , NAB1S+A1.NAB_) );
// 	}//(A1.NAReads_ == 1)

// }//set_up_block_info



// void ::cchem::rimp2_gradient::detail::RIAsync::set_up_block_info_cols(ArrayData &A1,  size_t NAB1S){
// 	//almost same as set_up_block_info 1) Rows <--> Cols
// 	//								   2) order or push switched
// 	if(A1.NAReads_ == 1){
// 		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 		temp.first.push_back(0);temp.first.push_back(NAB1S*A1.NABCols_);
// 		temp.second.push_back(A1.NABRows_);temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABCols_);
// 		A1.AReadPos_.push_back(temp);
// 		A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S,NAB1S+A1.NAB_) );
// 	}else{
// 		for(int i =0; i<A1.NAReads_-1; i++){
// 			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 			temp.first.push_back(0);temp.first.push_back( (NAB1S+i*A1.MaxNAB_)*A1.NABCols_);
// 			temp.second.push_back(A1.NABRows_);temp.second.push_back( (NAB1S+(i+1)*A1.MaxNAB_)*A1.NABCols_);
// 			A1.AReadPos_.push_back(temp);
// 			A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+i*A1.MaxNAB_ , NAB1S+(i+1)*A1.MaxNAB_) );
// 		}//i
// 		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// 		temp.first.push_back(0);temp.first.push_back( (NAB1S+(A1.NAReads_-1)*A1.MaxNAB_) *A1.NABCols_);
// 		temp.second.push_back(A1.NABRows_);temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABCols_ );
// 		A1.AReadPos_.push_back(temp);
// 		A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+(A1.NAReads_-1)*A1.MaxNAB_ , NAB1S+A1.NAB_) );
// 	}//(A1.NAReads_ == 1)

// }//set_up_block_info

// //void ::cchem::rimp2_gradient::detail::RIAsync::set_up_block_info_cols(ArrayData &A1,  size_t NAB1S){
// //
// //	if(A1.NAReads_ == 1){
// //		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// //		temp.first.push_back(0); temp.first.push_back(NAB1S*A1.NABCols_);
// //		temp.second.push_back( A1.NABRows_); temp.second.push_back( (NAB1S+A1.NAB_)*A1.NABCols_);
// //		A1.AReadPos_.push_back(temp);
// //		A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+A1.NAB_,NAB1S) );
// //	}else{
// //		for(int i =0; i<A1.NAReads_-1; i++){
// //			std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// //			temp.first.push_back( 0); temp.first.push_back((NAB1S+i*A1.MaxNAB_)*A1.NABCols_);
// //			temp.second.push_back( A1.NABRows_); temp.second.push_back((NAB1S+(i+1)*A1.MaxNAB_)*A1.NABCols_);
// //			A1.AReadPos_.push_back(temp);
// //			A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+(i+1)*A1.MaxNAB_, NAB1S+i*A1.MaxNAB_) );
// //		}//i
// //		std::pair< std::vector<size_t>, std::vector<size_t> > temp;
// //		temp.first.push_back( 0); temp.first.push_back((NAB1S+(A1.NAReads_-1)*A1.MaxNAB_) *A1.NABCols_);
// //		temp.second.push_back( A1.NABRows_ ); temp.second.push_back((NAB1S+A1.NAB_)*A1.NABRows_);
// //		A1.AReadPos_.push_back(temp);
// //		A1.ABlocks_.push_back( std::pair<size_t,size_t> ( NAB1S+A1.NAB_, NAB1S+(A1.NAReads_-1)*A1.MaxNAB_  ) );
// //	}//(A1.NAReads_ == 1)
// //
// //}//set_up_block_info


// void ::cchem::rimp2_gradient::detail::RIAsync::PrintArray(int i){

// 	ArrayData &A = Array_[i];
// 	std::cout << std::endl << "printing information for Array_[" << i << "]:" << std::endl ;
// 	std::cout << "  full array "<< i<<" looks like " << A.Array_->shape()[0] << " x " << A.Array_->shape()[1] << std::endl;
// 	std::cout << "  single array block size: " << A.NABRows_  << " x " <<  A.NABCols_  << std::endl;
// 	std::cout << "  need to perform " << A.NAReads_ << " NAReads_" << std::endl;
// 	std::cout << "  MaxNAB_ " << A.MaxNAB_ << std::endl;
// 	std::cout << "  array reads look like: " << std::endl;
// 	for(int i = 0; i < A.NAReads_; i++){
// 		std::cout << "    " <<"["<< i <<"]"<< "   { " << A.AReadPos_[i].first[0] << " " << A.AReadPos_[i].first[1] << " }"
// 				<<" --> { "<< A.AReadPos_[i].second[0] << " " << A.AReadPos_[i].second[1] << " }" << std::endl;
// 	}
// 	std::cout << "  array blocks range like: " << std::endl;
// 	for(int i = 0; i < A.NAReads_; i++){
// 		std::cout << "    " <<"["<< i  <<"]"<< "   " << A.ABlocks_[i].first << " : " << A.ABlocks_[i].second << std::endl;
// 	}
// 	std::cout << std::endl;
// }//Print




