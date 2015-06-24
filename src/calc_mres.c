/*------------------------------------------------------------------------
 *   Calculate Data Residuals                                  
 *   last update 29/03/08, D.Koehn
 *  ----------------------------------------------------------------------*/
#include "fd.h"

float calc_mres(int nx1, int nx2, int ny1, int ny2, float ** u, float ** u0){

/* declaration of variables */
int i,j;
float m2, m20, M2tmp;  

M2tmp=0.0;
m20=0.0;

for (j=ny1;j<=ny2;j++){
   for (i=nx1;i<=nx2;i++){

	   M2tmp += (u[j][i] - u0[j][i]) * (u[j][i] - u0[j][i]);
	   m20 += u0[j][i]*u0[j][i];

   }
}

m2=M2tmp/m20;

return m2;
}
