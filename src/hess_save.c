#include "fd.h"
void hess_save(int nt, float ** uy, float ** uyx, float ** vyp1, float *** hess_eyy, float *** hess_eyx, float *** hess_utty){
  
	extern float DT,TIME;
        extern float FC_HESSIAN;        	
	extern int NX, NY, IDX, IDY, DTINV;
	int i, j, k, l;
		
	        for (i=1;i<=NX;i=i+IDX){
		        for (j=1;j<=NY;j=j+IDY){
			
			     hess_eyy[j][i][nt]=uy[j][i];
			     hess_eyx[j][i][nt]=0.5 *uyx[j][i];
			     hess_utty[j][i][nt]=vyp1[j][i];
					  		  
			}
		}				
}