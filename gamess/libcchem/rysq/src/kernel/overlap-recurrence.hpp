/*
 * overlap-recurrence.hpp
 *
 *  Created on: Sep 3, 2015
 *      Author: luke
 */

#ifndef LIBCCHEM_RYSQ_SRC_KERNEL_OVERLAP_RECURRENCE_HPP_
#define LIBCCHEM_RYSQ_SRC_KERNEL_OVERLAP_RECURRENCE_HPP_

template<int N, int M>
struct overlap_recurrence {

    template<typename T, typename U, typename V>
    static void apply(int m,const int n, double aa12, double x0, double y0, double z0,
		      V *__restrict ri, V *__restrict rj, const U *__restrict H, const U *__restrict W,
		      T *__restrict Gx, T *__restrict Gy, T *__restrict Gz){

#define Gx(i,j) (Gx[(i)+(M*j)])
#define Gy(i,j) (Gy[(i)+(M*j)])
#define Gz(i,j) (Gz[(i)+(M*j)])

    	double xint = 0.0;
    	double yint = 0.0;
    	double zint = 0.0;

    	for (int a = 0; a < N; a++){
	    double dum = W[a];
	    double px = dum;
	    double py = dum;
	    double pz = dum;
	    dum = H[a]*aa12;
	    double ptx = dum + x0;
	    double pty = dum + y0;
	    double ptz = dum + z0;
	    double ax = ptx-ri[0];
	    double ay = pty-ri[1];
	    double az = ptz-ri[2];
	    double bx = ptx-rj[0];
	    double by = pty-rj[1];
	    double bz = ptz-rj[2];

	    for (int i = 0; i < m; i++){
		px = px*ax;
		py = py*ay;
		pz = pz*az;
	    }//i

	    for (int j = 0; j < n; j++){
		px = px*bx;
		py = py*by;
		pz = pz*bz;
	    }//j

	    xint += px;
	    yint += py;
	    zint += pz;

    	}//a

    	Gx(m,n) = xint*aa12;
    	Gy(m,n) = yint*aa12;
    	Gz(m,n) = zint*aa12;

    }//apply

#undef Gx
#undef Gy
#undef Gz

};//namespace overlap_recurrence



#endif /* LIBCCHEM_RYSQ_SRC_KERNEL_OVERLAP_RECURRENCE_HPP_ */
