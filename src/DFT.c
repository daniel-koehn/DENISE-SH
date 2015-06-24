#include "fd.h"
void DFT(int ishot, int nt, float ** vy, float ** syy, float ** sxy, float *** green_vy, float *** greeni_vy, float *** green_syy, float *** greeni_syy,
         float *** green_sxy, float *** greeni_sxy){
  
	extern float DT,TIME;
        extern float FC_HESSIAN;        	
	extern int NX, NY, IDX, IDY, DTINV;
	int i, j, k, l;
	double trig1,trig2;
	double t=0.0;
	const double pi=4.0*atan(1.0);
  
	t=nt*DT;

		trig1=cos(2.0*t*FC_HESSIAN*M_PI)*DT*DTINV;
		trig2=sin(2.0*t*FC_HESSIAN*M_PI)*DT*DTINV;
		
		/*printf("M_PI = %e \n",M_PI);*/
		
	        for (i=0;i<=NX+1;i=i+IDX){
		        for (j=0;j<=NY+1;j=j+IDY){
			                  
			    green_vy[j][i][ishot]+=vy[j][i]*trig1;
			    green_syy[j][i][ishot]+=syy[j][i]*trig1;
			    green_sxy[j][i][ishot]+=sxy[j][i]*trig1;
					  
			    greeni_vy[j][i][ishot]+=vy[j][i]*trig2;
			    greeni_syy[j][i][ishot]+=syy[j][i]*trig2;
			    greeni_sxy[j][i][ishot]+=sxy[j][i]*trig2;
					  
		
			}
		}				
}