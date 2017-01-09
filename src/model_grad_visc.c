/*
 *   1D gradient model from RAJZEL
 * 
 *
 *   Daniel Koehn
 *   Kiel, 04.01.2017
 *
 */

#include "fd.h"

void model(float  **  rho, float **  u, float **  taus, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
	extern float DH, *FL, TAU, DT;
		/* local variables */
	float vs, rhov, ts, tp, muv, piv, *pts;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE]; 
	
	/* parameters for gradient model */
	const float vs0=170.0, rho0=1900.0, Qs0=10.0, grad0=2.80;
	
	/*-----------------------------------------------------------------------*/
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	        eta[l]=DT/pts[l];
	}
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
			    /* 1D linear gradient model */
                            vs = grad0 * ((j-1)*DH) + vs0;
			    ts = 2.0 / Qs0;
                            rhov = rho0;   
				
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					u[jj][ii]=vs;
					rho[jj][ii]=rhov;
					taus[jj][ii]=ts;
					
				}
			}
		}	

		
sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
writemod(modfile,rho,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_0.bin",INV_MODELFILE);
writemod(modfile,u,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

free_vector(pts,1,L);
}



