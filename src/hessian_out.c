/* Output of the Hessian */
/*                                                 */
/* Daniel Koehn */
/* Kiel, the 16th of May 2013 */

#include "fd.h"

void hessian_out(float ** hessian_u, float ** hessian_rho, float ** u, float ** rho){
         	
extern int NX, NY, IDX, IDY, MYID, POS[4], INVMAT1;
extern char JACOBIAN[STRING_SIZE];

/* local variables */
int i, j, k, l;
char jac[STRING_SIZE];

FILE *FP4;  

/* calculate Hessian for different model parametrizations */
if(INVMAT1==1){
  
  for (i=1;i<=NX;i=i+IDX){    
     for (j=1;j<=NY;j=j+IDY){
     
         hessian_rho[j][i] += u[j][i] * u[j][i] * u[j][i] * u[j][i] * hessian_u[j][i]; 
         hessian_u[j][i] = 4.0 * rho[j][i] * rho[j][i] * u[j][i] * u[j][i] * hessian_u[j][i];
     
     }
  }
  
}


/* save HESSIAN for mu */
/* ----------------------- */
sprintf(jac,"%s_hessian_u.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){ 
        fwrite(&hessian_u[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);
    
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */   
sprintf(jac,"%s_hessian_u",JACOBIAN);
if (MYID==0) mergemod(jac,3);

/* save HESSIAN for rho */   
/* ----------------------- */
sprintf(jac,"%s_hessian_rho.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){ 
       fwrite(&hessian_rho[j][i],sizeof(float),1,FP4);
   }
}
    
fclose(FP4);
    
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */
sprintf(jac,"%s_hessian_rho",JACOBIAN);
if (MYID==0) mergemod(jac,3);               

}
