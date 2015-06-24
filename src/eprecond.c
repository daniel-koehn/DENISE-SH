#include "fd.h"
void eprecond(float ** W, float ** vy){
  
	extern int NX, NY, IDX, IDY;
	int i, j;
	
		
	        for (i=1;i<=NX;i=i+IDX){
		        for (j=1;j<=NY;j=j+IDY){
					     
			       W[j][i]+=vy[j][i]*vy[j][i];
		
			    }
		    }				
}
