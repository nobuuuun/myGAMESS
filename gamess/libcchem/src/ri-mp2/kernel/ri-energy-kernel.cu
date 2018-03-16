/* Includes, system */
#include <stdio.h>
#include <stdlib.h>

/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>


__global__ void invokeEnergyAccumIJ(int iblock, int jblock, int jblock_end,const double *d_C,
				     const double *d_ea,
				     const double *d_ev, 
				     int ns, int nvs,
				     double *d_E)
{

#define INDEX(i,j,dim) (i+j*dim)

  double energy = 0.0;

  for(unsigned int jocc = jblock; jocc < jblock_end; jocc++){
   const  double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
    
    const double dij = d_ea[iblock] + d_ea[jocc];
    double energy2 = 0.0;
      
    for(unsigned int a =ns +  blockIdx.x*blockDim.x + threadIdx.x;
  	a < nvs; 
	a += blockDim.x * gridDim.x){

  	const double dija = dij - d_ev[a];

  	for(unsigned int b = ns; b < nvs; b++){
	  
	    const double ab = ptr[ INDEX(a,b,nvs) ];
	    const double ba = ptr[ INDEX(b,a,nvs) ];

	    const double de = dija - d_ev[b];

	 energy2 += ab*(2*ab - ba)/de;

  	}//b
      
    }//a
    
    //energy += (double)2*energy2;
    energy += (1 + (iblock!=jocc) )*energy2;
    
 }//jocc

    d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy;

}//invokeEnergyAccumIJ


__global__ void invokeEnergyAccumOpenIJ(int iblock, int jblock, int jblock_end,const double *d_C,
				     const double *d_ea,
				     const double *d_ev, 
				     const double *d_EM, 
				     int ns, int nvs,
				     double *d_E)
{

#define INDEX(i,j,dim) (i+j*dim)

    double energy = 0.0;



    for(unsigned int jocc = jblock; jocc < jblock_end; jocc++){
	const double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
    
	const double dij = d_ea[iblock] + d_ea[jocc];
	double energy2 = 0.0;

	//evaluate 1 (standard)
	for(unsigned int a =ns +  blockIdx.x*blockDim.x + threadIdx.x;
	    a < nvs; 
	    a += blockDim.x * gridDim.x){
	    const double dija = dij - d_ev[a];
	    for(unsigned int b = ns; b < nvs; b++){
		const double ab = ptr[ INDEX(a,b,nvs) ];
		const double ba = ptr[ INDEX(b,a,nvs) ];
		const double de = dija - d_ev[b];
		energy2 += ab*(2*ab - ba)/de;
	    }//b
	}//a
	energy += (1 + (iblock!=jocc) )*energy2;

	// evaluate 2
	energy2 = 0.0;
	for(unsigned int a = blockIdx.x*blockDim.x + threadIdx.x;
	    a < ns; 
	    a += blockDim.x * gridDim.x ){
	    const double dija = dij - d_ev[a] - d_EM[a];
	    for(unsigned int b = ns; b < nvs; b++){
		const double ab = ptr[ INDEX(a,b,nvs) ];
		const double ba = ptr[ INDEX(b,a,nvs) ];
		const double de = dija - d_ev[b];
		energy2 += (ab*ab - ab*ba + ba*ba)/de;
	    }//b
	}//a
	//	energy += energy2;
	energy += ( 1+ (iblock!=jocc) )*energy2;


	// evaluate 7
	energy2 = 0.0;
	for(unsigned int a = blockIdx.x*blockDim.x + threadIdx.x;
	    a < ns; 
	    a += blockDim.x * gridDim.x ){
	    const double dija = dij - d_ev[a] - d_EM[a];
	    for(unsigned int b = 0; b < ns; b++){
		const double ab = ptr[ INDEX(a,b,nvs) ];
		const double ba = ptr[ INDEX(b,a,nvs) ];
		const double de = dija - d_ev[b] - d_EM[b];
		energy2 += ( (ab-ba)*(ab-ba) )/de;
	    }//b
	}//a
	//	energy += energy2/4;
	energy += ( 1+ (iblock!=jocc) )*energy2/4;
	
    }//jocc
    
    //    d_E[blockIdx.x*blockDim.x + threadIdx.x] += 2*energy;
    d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy;


}//invokeEnergyAccumOpenIJ




__global__ void invokeEnergyAccumOpenIJ_GE_D(int iocc,int jst, int jblock, int jblock_end, const double *d_C,
					     const double *d_ea,
					     const double *d_ev,
					     const double *d_EM,
					     int nd, int ns, int nvs,
					     double *d_E)


{

#define INDEX(i,j,dim) (i+j*dim)

    //evaluate 3
  double energy = 0.0;

  for(unsigned int jocc = jst; jocc < jblock_end; jocc++){
    const double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
    
    const double dij = d_ea[iocc] + d_ea[jocc] - d_EM[iocc-nd] - d_EM[jocc-nd];
    double energy2 = 0.0;

    for(unsigned int a = ns + blockIdx.x*blockDim.x + threadIdx.x;
    	a < nvs; a += blockDim.x * gridDim.x ){
	const double dija = dij - d_ev[a];

    	for(unsigned int b = ns; b < nvs; b++){
    	    const double ab = ptr[ INDEX(a,b,nvs) ];
    	    const double ba = ptr[ INDEX(b,a,nvs) ];
    	    const double de = dija - d_ev[b];

	    energy2 += ((ab-ba)*(ab-ba))/de;

      }//b
    }//a
    energy += ( 1+ (iocc!=jocc) )*energy2;

  }//jocc

  d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy/4;

}//invokeEnergyAccumOpenIJ_GE_D










__global__ void invokeEnergyAccumOpenIJ_LT_D(int iocc, int jst, int jblock, int jblock_end, 
					     const double *d_C,
					     const double *d_ea,
					     const double *d_ev,
					     const double *d_EM,
					     int nd, int ns, int nvs,
					     double *d_E,
					     double *d_FOCK)

{

#define INDEX(i,j,dim) (i+j*dim)

    //evaluate 6
  double energy = 0.0;



  for(unsigned int jocc = jst; jocc < jblock_end; jocc++){
      const double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
    
      const double dij = d_ea[iocc] + d_ea[jocc] - d_EM[iocc-nd];
      double energy2 = 0.0;

      for(unsigned int a = ns + blockIdx.x*blockDim.x + threadIdx.x;
      	  a < nvs; a += blockDim.x * gridDim.x ){
      	  const double dija = dij - d_ev[a];
      	  for(unsigned int b = ns; b < nvs; b++){
      	      const double ab = ptr[ INDEX(a,b,nvs) ];
      	      const double ba = ptr[ INDEX(b,a,nvs) ];
      	      const double de = dija - d_ev[b];
      	      energy2 += (ab*(2*ab-ba))/de;
      	  }//b
      }//a
      energy += energy2;

      //evaluate 4
      energy2 = 0.0;
      for(unsigned int a = blockIdx.x*blockDim.x + threadIdx.x;
      	  a < ns; a += blockDim.x * gridDim.x ){
      	  const double dija = dij - d_ev[a] - d_EM[a];
      	  for(unsigned int b = ns; b < nvs; b++){
      	 //     const double ab = ptr[ INDEX(a,b,nvs) ];
      	      const double ba = ptr[ INDEX(b,a,nvs) ];
      	      const double de = dija - d_ev[b];
      	      energy2 += (ba*ba)/de;
      	  }//b
      }//a
      energy += energy2;


      //evaluate 5 (this is intermediate data)
      int nv = nvs - ns;
      int single = iocc - nd;
      for(unsigned int a = ns + blockIdx.x*blockDim.x + threadIdx.x;
       	  a < nvs;
	  a += blockDim.x * gridDim.x ){
	  const  double ab = ptr[ INDEX(a, single ,nvs) ];
	  d_FOCK[ INDEX( (a-ns) , jocc, nv ) ] += ab;
      }//a
      

      
  }//jocc

  d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy;

  }//invokeEnergyAccumOpenIJ_GE_D




extern "C" void deviceCudaEnergyIJ(int iblock, int jblock, int jblock_end,
				   const double *d_C,
				   const double *d_ea,
				   const double *d_ev,
				   const double* d_EM,
				   int nd, int ns, int nvs,
				   double *d_E,
				   double *d_FOCK,
				   int cudaBlocks,
				   int cudaThreadsPB,
				   cudaStream_t stream)

{
  //if you change block/threads, make sure to adjust d_E


    if(!ns){

	//standard evaluate (1)
	invokeEnergyAccumIJ<<<8, 64,0, stream>>>
	    (iblock, jblock, jblock_end, d_C, d_ea, d_ev, ns, nvs, d_E);

	return;
    }


    // if(iblock < nd ){

    // 	//standard evaluate (1)
    // 	//	invokeEnergyAccumIJ<<<8, 64,0, stream>>>
    // 	    //	    (iblock, jblock, jblock_end, d_C, d_ea, d_ev, ns, nvs, d_E);

    // }//(iocc < nd )


    // if(!ns)return;


    if(iblock < nd ){
	//evaluate 1 & 2 & 7

	 	invokeEnergyAccumOpenIJ<<<8, 64,0, stream>>>
	 (iblock, jblock, jblock_end, d_C, d_ea, d_ev, d_EM, ns, nvs, d_E);


    }else if( jblock < nd) {

	int jst = jblock;
	int jend = jblock_end;
	if(jblock_end > nd) jend = nd; //-1;

	//evaluate 4 & 5 & 6
	invokeEnergyAccumOpenIJ_LT_D<<<8, 64,0,stream>>>
	    (iblock, jst, jblock, jend, 
	     d_C, 
	     d_ea, 
	     d_ev, 
	     d_EM, 
	     nd, ns, nvs, 
	     d_E, 
	     d_FOCK);

    }else{

	//	if(jblock_end > nd){
	    int jst = nd;
	    if(jblock > nd)jst = jblock;
	    //evaluate 3

	    invokeEnergyAccumOpenIJ_GE_D<<<8, 64,0,stream>>>
	    	(iblock, jst, jblock, jblock_end, d_C, d_ea, d_ev, d_EM, nd, ns, nvs, d_E);

	    //	}//(jblock_end > nd)

    }//(iocc < nd)
    


}//deviceCudaEnergyIJ























__global__ void invokeEnergyAccumII(int iocc,int  jblock,const double *d_C,
				     const double *d_ea,
				     const double *d_ev,
				     int ns, int nvs,
				     double *d_E)
{

#define INDEX(i,j,dim) (i+j*dim)

    //evaluate 1 (standard term)
  double energy = 0.0;

  for(unsigned int jocc = jblock; jocc <= iocc; jocc++){
    const double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
    
    const double dij = d_ea[iocc] + d_ea[jocc];
    double energy2 = 0.0;

    for(unsigned int a = ns + blockIdx.x*blockDim.x + threadIdx.x;
	a < nvs; a += blockDim.x * gridDim.x ){
      
      const double dija = dij - d_ev[a];

      for(unsigned int b = ns; b < nvs; b++){
	  
	const double ab = ptr[ INDEX(a,b,nvs) ];
	const double ba = ptr[ INDEX(b,a,nvs) ];
	const double de = dija - d_ev[b];
	
	energy2 += ab*(2*ab - ba)/de;

      }//b
	
    }//a

    energy += ( 1+ (iocc!=jocc) )*energy2;
    
  }//jocc

  d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy;

}//invokeEnergyAccumII




__global__ void invokeEnergyAccumOpenII(int iocc,int  jblock,const double *d_C,
					const double *d_ea,
					const double *d_ev,
					const double *d_EM,
					int ns, int nvs,
					double *d_E)
{

#define INDEX(i,j,dim) (i+j*dim)

    //evaluate 2
  double energy = 0.0;

  for(unsigned int jocc = jblock; jocc <= iocc; jocc++){
    const double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
    
    const double dij = d_ea[iocc] + d_ea[jocc];
    double energy2 = 0.0;

    //term2
    for(unsigned int a = blockIdx.x*blockDim.x + threadIdx.x;
    	a < ns; a += blockDim.x * gridDim.x ){
    	const double dija = dij - d_ev[a] - d_EM[a];
    	for(unsigned int b = ns; b < nvs; b++){
    	    const double ab = ptr[ INDEX(a,b,nvs) ];
    	    const double ba = ptr[ INDEX(b,a,nvs) ];
    	    const double de = dija - d_ev[b];
    	    //	    energy2 += ab*(2*ab - ba)/de;
    	    energy2 += (ab*ab - ab*ba + ba*ba)/de;
      }//b
    }//a
       energy += ( 1+ (iocc!=jocc) )*energy2;
    
       //evaluate 7
    energy2 = 0.0;
    for(unsigned int a = blockIdx.x*blockDim.x + threadIdx.x;
    	a < ns; a += blockDim.x * gridDim.x ){
    	const double dija = dij - d_ev[a] - d_EM[a];
    	for(unsigned int b = 0; b < ns; b++){
    	    const double ab = ptr[ INDEX(a,b,nvs) ];
    	    const double ba = ptr[ INDEX(b,a,nvs) ];
    	    const double de = dija - d_ev[b] - d_EM[b];
    	    //	    energy2 += ab*(2*ab - ba)/de;
    	    energy2 += ( (ab-ba)*(ab-ba) )/de;
      }//b
    }//a
    energy += ( 1+ (iocc!=jocc) )*energy2/4;
    
  }//jocc

  d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy;

}//invokeEnergyAccumOpenII


__global__ void invokeEnergyAccumOpenII_GE_D(int iocc,int jst,int  jblock,const double *d_C,
					     const double *d_ea,
					     const double *d_ev,
					     const double *d_EM,
					     int nd, int ns, int nvs,
					     double *d_E)
{

#define INDEX(i,j,dim) (i+j*dim)

    //evaluate 3
  double energy = 0.0;
  
  //  for(unsigned int jocc = jblock; jocc < iocc; jocc++){
  for(unsigned int jocc = jst; jocc < iocc; jocc++){
      const double *ptr = &d_C[ nvs*nvs*(jocc-jblock) ];
  
      const double dij = d_ea[iocc] + d_ea[jocc] - d_EM[iocc-nd] - d_EM[jocc-nd];
      double energy2 = 0.0;

      for(unsigned int a = ns + blockIdx.x*blockDim.x + threadIdx.x;
	  a < nvs; a += blockDim.x * gridDim.x ){
	  const double dija = dij - d_ev[a];
	  
	  for(unsigned int b = ns; b < nvs; b++){
	      const double ab = ptr[ INDEX(a,b,nvs) ];
	      const double ba = ptr[ INDEX(b,a,nvs) ];
	      const double de = dija - d_ev[b];
	      
	      energy2 += ((ab-ba)*(ab-ba))/de;
	      
	  }//b
      }//a
      energy += ( 1+ (iocc!=jocc) )*energy2;
      
  }//jocc

  d_E[blockIdx.x*blockDim.x + threadIdx.x] += energy/4;

}//invokeEnergyAccumOpenII_GE_D






extern "C" void deviceCudaEnergyII(int iocc, int jblock, const double *d_C,
				   const double *d_ea,
				   const double *d_ev, 
				   const double * d_EM,
				   int nd, int ns, int nvs,
				   double *d_E,
				   double *d_FOCK,
				   int cudaBlocks,
				   int cudaThreadsPB,
				   cudaStream_t stream)
{

  //if you change block/threads, make sure to adjust d_E



    if(!ns){

	int jblock_end = iocc+1;
	invokeEnergyAccumIJ<<<8, 64,0, stream>>>
	    (iocc, jblock, jblock_end, d_C, d_ea, d_ev, ns, nvs, d_E);
	return;
    }


    // if(iocc < nd ){
    // 	//evaluate 1 (standard)

    // 	//	invokeEnergyAccumII<<<8, 64,0,stream>>>
    // 	//	    (iocc, jblock, d_C, d_ea, d_ev, ns, nvs, d_E);

    // 		// int jblock_end = iocc+1;
    // 		// invokeEnergyAccumIJ<<<8, 64,0, stream>>>
    // 		//     (iocc, jblock, jblock_end, d_C, d_ea, d_ev, ns, nvs, d_E);


    // }//(iocc < nd)

    //    if(!ns)return;

    if(iocc < nd ){
	//evaluate 1 & 2 & 7

	//	invokeEnergyAccumOpenII<<<8, 64,0,stream>>>
	//	    (iocc, jblock, d_C, d_ea, d_ev, d_EM, ns, nvs, d_E);

	int jblock_end = iocc +1;
	invokeEnergyAccumOpenIJ<<<8, 64,0, stream>>>
	    (iocc, jblock, jblock_end, 
	     d_C, d_ea, d_ev, d_EM, 
	     ns, nvs, 
	     d_E);


    }else if( jblock < nd) {

    	int jst = jblock;
    	int jend = iocc;
    	if(iocc > nd) jend = nd; //-1;

    	invokeEnergyAccumOpenIJ_LT_D<<<8, 64,0,stream>>>
    	    (iocc, jst, jblock, jend, 
    	     d_C, d_ea, d_ev, d_EM, 
   	     nd, ns, nvs, 
    	     d_E, d_FOCK);
	

	//    }else{


	//	jst = nd;
	//	if(jblock > nd)jst = jblock;
	//	//evaluate 3
	//	invokeEnergyAccumOpenII_GE_D<<<8, 64,0,stream>>>
	//	    (iocc, jst, jblock, d_C, d_ea, d_ev, d_EM, nd, ns, nvs, d_E);

	int jblock_end = iocc;
	jst = nd;
	if(jblock > nd)jst = jblock;
	//evaluate 3

	invokeEnergyAccumOpenIJ_GE_D<<<8, 64,0,stream>>>
	    (iocc, jst, jblock, jblock_end, d_C, d_ea, d_ev, d_EM, nd, ns, nvs, d_E);

    }//logic
 
}//deviceCudaEnergyII












//#include <stdio.h>

__global__ void invokeDeviceAccum(double *d_E){

    int tid = threadIdx.x;

    for  (unsigned int s=1; s < blockDim.x; s *= 2) {
    
	int index = 2 * s * tid;
	if (index < blockDim.x) {
	    d_E[index] += d_E[index + s];
	}
	__syncthreads();
    }//s

    //    if(tid == 0)  printf("energy:%f  %f \n",  d_E[0]  );
}

extern "C" void deviceAccum(double *d_E, 
			    int cudaBlocks,
			    int cudaThreadsPB,
			    cudaStream_t stream){

    invokeDeviceAccum<<<1, 512,0, stream>>>(d_E);

}//deviceAccum





































