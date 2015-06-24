/* Calculation of the gradients in the frequency domain */
/*                                                      */
/* Daniel Koehn */
/* Kiel, the 11th of October 2013 */

#include "fd.h"
#include <complex.h>

void grad_FD(float ** grad_vy, float ** gradi_vy, float ** grad_syz, float ** gradi_syz, float ** grad_sxy, float ** gradi_sxy, 
             float ** grad_vyb, float ** gradi_vyb, float ** grad_syzb, float ** gradi_syzb, float ** grad_sxyb, float ** gradi_sxyb, float ** pu, float ** prho,
             float ** waveconv_u_shot, float ** waveconv_rho_shot){
  
extern float FC_START;        	
extern int NX, NY, IDX, IDY, INVMAT1;

/* local variables */
int i, j;

float complex cv_y, cv_yb, csyz, csyzb, csxy, csxyb;
float muss, omega;

/* calculate angular frequency */
omega = 2.0 * M_PI * FC_START;

/* assemble gradients for mu and rho */
/* ----------------------------------------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
            
   /* assemble complex wavefields */
   cv_y  = (grad_vy[j][i] + gradi_vy[j][i] * I);                 
   cv_yb = (grad_vyb[j][i] + gradi_vyb[j][i] * I);                 

   csyz  = (grad_syz[j][i] + gradi_syz[j][i] * I);                 
   csyzb = (grad_syzb[j][i] + gradi_syzb[j][i] * I);                 
   
   csxy  = (grad_sxy[j][i] + gradi_sxy[j][i] * I);                 
   csxyb = (grad_sxyb[j][i] + gradi_sxyb[j][i] * I);                 
                    
   if(INVMAT1==1){
     muss = prho[j][i] * pu[j][i] * pu[j][i];
   }
                    
   if(INVMAT1==3){
     muss = pu[j][i];
   }
                      
   waveconv_u_shot[j][i]   = creal(omega*I*((csxy*csxyb)+(csyz*csyzb))/(muss*muss));  
   waveconv_rho_shot[j][i] = creal(-omega*I*(cv_y*cv_yb));
                                      
   }
}

}
