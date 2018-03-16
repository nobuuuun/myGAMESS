#include "config.h"

#include "god.hpp"

extern "C" {


#if defined(HAVE_ATLAS)
#include <cblas.h>
#include <clapack.h>
#endif


#if defined(HAVE_MKL)
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#endif

#if defined(HAVE_OPENBLAS)
#include <cblas.h>
#include <lapacke.h>
#endif

    //        #include "/opt/local/include/cblas.h"    
    //        #include "/opt/local/include/clapack.h"  

    
    //    #include "/u1/luke/opt/ATLAS/include/cblas.h"
    //    #include "/u1/luke/opt/ATLAS/include/clapack.h"

}
