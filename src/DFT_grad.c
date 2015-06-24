#include "fd.h"
void DFT_grad(int nt, float ** vy, float ** syy, float ** sxy, float ** green_vy, float ** greeni_vy, float ** green_syy, float ** greeni_syy,
         float ** green_sxy, float ** greeni_sxy, int back){
  
	extern float DT,TIME;
        extern float FC_START;        	
	extern int NX, NY, IDX, IDY, DTINV;
	int i, j, k, l;
	double trig1,trig2;
	double t=0.0;
	/*const double pi=4.0*atan(1.0);*/
  
	if(back==0){t=nt*DT;} /* forward wavefield */
	if(back==1){t=TIME-nt*DT;} /* backward wavefield */
        

		trig1=cos(2.0*t*FC_START*M_PI)*DT*DTINV;
		trig2=sin(2.0*t*FC_START*M_PI)*DT*DTINV;
		
	        for (i=1;i<=NX;i=i+IDX){
		        for (j=1;j<=NY;j=j+IDY){
			                  
			    green_vy[j][i]+=vy[j][i]*trig1;
			    green_syy[j][i]+=syy[j][i]*trig1;
			    green_sxy[j][i]+=sxy[j][i]*trig1;
					  
			    greeni_vy[j][i]+=vy[j][i]*trig2;
			    greeni_syy[j][i]+=syy[j][i]*trig2;
			    greeni_sxy[j][i]+=sxy[j][i]*trig2;
					  
		
			}
		}				
}
