/*
 *   Model with slow formation and fluid-filled borehole
 *   last update 11.12.07, O. Hellwig
 */

#include "fd.h"

void model_elastic(float  **  rho, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DH;
	extern int   NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float rhov, muv, vs, y0, t, de, zplat1, zplat2, rplat;
	float *pts, ts, sumu, ws, *ri, *d, *ti, *dl, *vsl, z, r;
	float **checks, **checkrho; 
	int   i, j, l, ii, jj, nk, k, nl;
	char filename_mu[STRING_SIZE];
	char filename_rho[STRING_SIZE]; 
				
	sprintf(filename_mu,"%s.mu",MFILE);
	sprintf(filename_rho,"%s.rho",MFILE);
	
	/* parameters for borehole */  
        const float vs1=3500.0, rho1=7850.0; /* steel */
        	
	/* parameters for layer 1 */
	const float vs2=250.0, rho2=1800.0;
	
	/* parameters for air */
        const float vsa=1e-6, rhoa=1.25;
	
	/* parameters for water */    
        const float vsw=1e-6, rhow=1000.0;         
	
	/* borehole radius */
	float R, a , b, dw; 
	const float rshift=10.0;
	
	/* parameters for the checkquerboard CTS test */
	int nyc, nxc;
	int idcx, idcy;
	float cxsign, cysign;
	
	nyc = 50;
	nxc = 50;
	
	cxsign = 1.0;
	cysign = 1.0;
	
	checks = matrix(1,nyc,1,nxc);
	checkrho = matrix(1,nyc,1,nxc);
	
	/* define parameters for CTS checkerboard test */
        for (i=1;i<=nxc;i++){
	    for (j=1;j<=nyc;j++){
	       checks[j][i]=100.0 ;
	       checkrho[j][i]=0.0;
	    }
	}	
	
	/*-----------------------------------------------------------------------*/

	/* loop over global grid */
	idcx=1;
	for (i=1;i<=NXG;i++){
           idcy=1;
                cysign=-cysign;
		for (j=1;j<=NYG;j++){
	                
			  vs = 164.0;
			  rhov = 1900.0;
	                
	                  /*if((i>nxc)&&(i<NXG-nxc)&&(j<NYG-nyc)){
	                  vs = 1333.0 + cxsign * cysign * checks[idcy][idcx];
	                  vp = 2310.0 + cxsign * cysign * checkp[idcy][idcx];
			  rhov = 1900.0 + cxsign * cysign * checkrho[idcy][idcx];
			  }
			  idcy++;*/

			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii = i-POS[1]*NX;
				jj = j-POS[2]*NY;
				
				u[jj][ii]    = vs;
				rho[jj][ii]  = rhov;
			}
			
		        if(idcy>nyc){
		           idcy = 1;
		           cysign = - cysign;
		        }	
		}
	   idcx++;	
	   if(idcx>nxc){
	      idcx = 1;
	      cxsign = - cxsign;
	   }	
	}	

	/* each PE writes his model to disk */
	writemod(filename_rho,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_rho,3);
	                        
	/* each PE writes his model to disk */
	writemod(filename_mu,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename_mu,3);
	                        
}


