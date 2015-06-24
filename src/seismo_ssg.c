/*------------------------------------------------------------------------
 *   store amplitudes (particle velocities or pressure or curl and div) 
 *    at receiver positions in arrays
 *   last update 27/12/01, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvy, float **vy, float **syy, float **u, float *hc){ 
		
	extern int NDT, SEISMO, FDORDER;	
	extern float DH;
	int i,j, itr, ins, nxrec, nyrec, m;

	/*ins=lsamp/NDT;*/
	ins=lsamp;
	for (itr=1;itr<=ntr;itr++){
		nxrec=recpos[1][itr];
		nyrec=recpos[2][itr];
		switch (SEISMO){
		case 1 : 
			sectionvy[itr][ins]=vy[nyrec][nxrec];
			break;
		
		}

	}
}
