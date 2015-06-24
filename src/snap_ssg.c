/*------------------------------------------------------------------------
 *   Write 2D snapshot for current timestep  to file                                   
 *   last update 24/05/2002
 *
 *  T. Bohlen
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void snap(FILE *fp, int nt, int nsnap, float **vy, float **syy, float **u, float *hc){

	/* 
		different data formats of output:
		SNAP_FORMAT=1  :  SU (IEEE)
		SNAP_FORMAT=2  :  ASCII
		SNAP_FORMAT=3  :  BINARY (IEEE)
		
		different types:
		SNAP=1 : values in vx and vy
		SNAP=2 : -(vx+vy) (pressure field)
		SNAP=3 : divergence of vx and vy (energy of compressional waves)
		         and curl of vx and vy (energy of shear waves)
		SNAP=4 : both particle velocities (type=1) and energy (type=3)
		*/


	int i,j, m, fdoh, nd;
	char snapfile_y[STRING_SIZE];
	char ext[8], wm[2];
	FILE *fpy1;

	extern float DH, DT;
	extern char SNAP_FILE[STRING_SIZE];
	extern int NX, NY,  SNAP_FORMAT, SNAP, FDORDER;
	extern int MYID, POS[3], IDX, IDY;
	
	fdoh = FDORDER/2;

	switch(SNAP_FORMAT){
	case 1:
		sprintf(ext,".su");
		break;
	case 2:
		sprintf(ext,".asc");
		break;
	case 3:
		sprintf(ext,".bin");
		break;
	}
	
	sprintf(snapfile_y,"%s%s.y.%i%i",SNAP_FILE,ext,POS[1],POS[2]);
	fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);
	
	
	if (nsnap==1)
		sprintf(wm,"w");
	else
		sprintf(wm,"a");
		

	switch(SNAP){
	case 1 :
		fprintf(fp,"%s\n\n",snapfile_y);
		
		fpy1=fopen(snapfile_y,wm);

		for (i=1;i<=NX;i+=IDX){
			for (j=1;j<=NY;j+=IDY){
				writedsk(fpy1,vy[j][i],SNAP_FORMAT);
			}
	        }
	        
		fclose(fpy1);
	break;

	}


}


