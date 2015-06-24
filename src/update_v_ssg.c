/* $Id: update_v_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen 
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"



void update_v_hc(int nx1, int nx2, int ny1, int ny2, int nt, float ** vy, float **  vyp1, float **  vym1, float **  utty, float ** syy,
	float ** sxy, float **rjp, float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff,
	float *hc, int infoout,int sw){

	int i, j,l,fdoh,m;
	float amp, dtdh, drad, angle;
	float vxtmp, vytmp;
	
	extern float DT, DH, FW;
	double time1, time2;
	extern int MYID, QUELLTYP, CHECKPTREAD, FDORDER;
	extern FILE *FP;

	
	fdoh = FDORDER/2;
	dtdh = DT*DT/DH;
        
        drad = PI/180.0; 
        angle = 135.0;
         
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_v (printed by PE %d):\n",MYID);
		fprintf(FP," Updating particle velocities ...");
	}


	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */

	switch (FDORDER){
	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

				/* updating the y-component of the velocity (vy) */
				vyp1[j][i] = 2.0 * vy[j][i] - vym1[j][i] + (  hc[1]*(syy[j+1][i]-syy[j][i])
					    + hc[1]*(sxy[j][i]-sxy[j][i-1]))*dtdh*rjp[j][i];
			}
		}
		
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			     
			     vym1[j][i] = vy[j][i];
			     vy[j][i]   = vyp1[j][i];
			}
		}
		
		break;
		
	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
						
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					   )*dtdh*rjp[j][i];
			}
		}
		break;
		
	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
						
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					   )*dtdh*rjp[j][i];
			}
		}
		break;
		
	case 8:
	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
					    
				utty[j][i] = rjp[j][i]*(hc[1]*(syy[j+1][i]-syy[j][i])
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[4]*(syy[j+4][i]-syy[j-3][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4]))/DH; 	    
						
				vy[j][i] += DT*utty[j][i];
			}
		}
		
		
		
		break;

	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
						
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[4]*(syy[j+4][i]-syy[j-3][i])
					    + hc[5]*(syy[j+5][i]-syy[j-4][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4])
					    + hc[5]*(sxy[j][i+4]-sxy[j][i-5])
					   )*dtdh*rjp[j][i];
			}
		}
		break;

	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
	
						
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
					    + hc[2]*(syy[j+2][i]-syy[j-1][i])
					    + hc[3]*(syy[j+3][i]-syy[j-2][i])
					    + hc[4]*(syy[j+4][i]-syy[j-3][i])
					    + hc[5]*(syy[j+5][i]-syy[j-4][i])
					    + hc[6]*(syy[j+6][i]-syy[j-5][i])
					    + hc[1]*(sxy[j][i]  -sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4])
					    + hc[5]*(sxy[j][i+4]-sxy[j][i-5])
					    + hc[6]*(sxy[j][i+5]-sxy[j][i-6])
					   )*dtdh*rjp[j][i];
			}
		}
		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxtmp = 0;
				vytmp = 0;
				for (m=1; m<=fdoh; m++) {
							
					vytmp +=   hc[m]*( syy[j+m][i]   - syy[j-m+1][i] )
						 + hc[m]*( sxy[j][i+m-1] - sxy[j][i-m]   );
				}
					
				vy[j][i] += vytmp*dtdh*rjp[j][i];
			}
		}
		break;

	} /* end of switch(FDORDER) */
			
	
	/*if (FW>0.0)
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
		
			vy[j][i]*=absorb_coeff[j][i];
			sxy[j][i]*=absorb_coeff[j][i];
			syy[j][i]*=absorb_coeff[j][i];
		
	}}*/

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
