/*------------------------------------------------------------------------
 *   Initialize FD-FWI variables
 *  
 *   Daniel Koehn
 *   Kiel, 11th of October 2013
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_FWI(float ** grad_vy,float ** gradi_vy,float ** grad_syz,float ** gradi_syz,float ** grad_sxy,float ** gradi_sxy){

     register int i, j;
     extern int NX, NY;
		
     for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
		 
            grad_vy[j][i] = 0.0;
            gradi_vy[j][i] = 0.0;

            grad_syz[j][i] = 0.0;
            gradi_syz[j][i] = 0.0;

            grad_sxy[j][i] = 0.0;
            gradi_sxy[j][i] = 0.0;

		 
        }
      }
		       
}
