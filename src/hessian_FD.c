/* Calculation, assembly and output of the Hessian */
/*                                                 */
/* Daniel Koehn */
/* Kiel, the 3rd of september 2012 */

#include "fd.h"
#include <complex.h>

void hessian(int nshots, int SHOTINC, float *** green_vy, float *** greeni_vy, float *** green_syy, float *** greeni_syy,
             float *** green_sxy, float *** greeni_sxy, float ** prho, float ** pu, int iter){
  
extern float DT,DH,TIME;
extern float FC_HESSIAN;        	
extern int NX, NY, IDX, IDY, DTINV, INVMAT1, MYID, POS[4], FDORDER;
extern char JACOBIAN[STRING_SIZE];

/* local variables */
int i, j, k, l, ns_hess, ishot, irec, nd, NSRC_HESSIAN, NREC_HESSIAN, RECINC;
double trig1,trig2;
double t=0.0;
const double pi=4.0*atan(1.0);
char jac[STRING_SIZE];

float complex uttx, utty, exx, eyy, eyx, Gxx, Gyx, Gxy, Gyy, Gxxx, Gyxx, Gxyy, Gyyy, Gxyx, Gyyx;
float complex tmp_jac_lam, tmp_jac_mu, tmp_jac_rho, tmp_jac_vp, tmp_jac_vs, tmp_fft_fsignal;

float ** abs_green, omega, mulamratio, lamss, muss, HESS_SCALE;
float ** hessian, ** hessian_u, ** hessian_rho, ** hessian_lam, ** hessian_mu;

float hvxx, hvxxi, hvyy, hvyyi, hvxy, hvxyi,  hvyx, hvyxi;
float **exx_shot, **exxi_shot, **eyy_shot, **eyyi_shot, **eyx_shot, **eyxi_shot, **uttx_shot, **uttxi_shot, **utty_shot, **uttyi_shot;
float *psource_hess=NULL, *Hess_for_real=NULL, *Hess_for_complex=NULL;

FILE *FP4;

HESS_SCALE = 1.0;
RECINC = 1;
NSRC_HESSIAN=1;
NREC_HESSIAN=24;

nd = FDORDER/2 + 1;
abs_green = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

/* Diagonal elements of the Hessian*/
hessian_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd); 
hessian_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
eyy_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
eyyi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
eyx_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
eyxi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

utty_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd); 
uttyi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);


Hess_for_real = vector(1,1);
Hess_for_complex = vector(1,1);
        
for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
           hessian_u[j][i]=0.0;  
           hessian_mu[j][i]=0.0; 
           hessian_rho[j][i]=0.0;
    }
}


/* assemble Hessian */
/* ----------------------------------------------------------------- */

/* Circular frequency of the Hessian */
omega = 2.0*M_PI*FC_HESSIAN;

/* calculate absolute values of impulse responses */
for (ishot=1;ishot<=nshots;ishot=ishot+SHOTINC){

    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
        
           /*green_x = green_vx[j][i][ishot] + greeni_vx[j][i][ishot] * I;
           green_y = green_vy[j][i][ishot] + greeni_vy[j][i][ishot] * I;*/
                                                                        
           abs_green[j][i] = 1.0;
           
        }
    }

}    

/*printf("omega = %e \n",omega);
printf("NSRC = %d \n",NSRC_HESSIAN);*/

/*psource_hess=rd_sour(&ns_hess,fopen(SIGNAL_FILE,"r"));
FFT_data(psource_hess,Hess_for_real,Hess_for_complex,NT);
                                                         
MPI_Barrier(MPI_COMM_WORLD);                             

tmp_fft_fsignal = Hess_for_real[1] + Hess_for_complex[1] * I;*/

tmp_fft_fsignal = 1.0;

for (ishot=1;ishot<=nshots;ishot=ishot+SHOTINC){
       
       /* calculate spatial and temporal derivatives of the forward wavefield */ 
       for (i=1;i<=NX;i=i+IDX){
          for (j=1;j<=NY;j=j+IDY){
                                       
                 hvyy = (green_vy[j+1][i][ishot]-green_vy[j][i][ishot])/DH;
                hvyyi = (greeni_vy[j+1][i][ishot]-greeni_vy[j][i][ishot])/DH;
                           
                 hvyx = (green_vy[j][i+1][ishot]-green_vy[j][i][ishot])/DH; 
                hvyxi = (greeni_vy[j][i+1][ishot]-greeni_vy[j][i][ishot])/DH;
                               
                /* calculate strain tensors and integrate FD-wavefield */
       
               eyy_shot[j][i] = hvyy;
               eyyi_shot[j][i] = hvyyi;
               
               eyxi_shot[j][i] = hvyxi;
               eyx_shot[j][i] =  hvyx;
                              
               utty_shot[j][i] = -hvyyi*omega;             
               uttyi_shot[j][i] = hvyy*omega;
               
          }
       }

 for (irec=1;irec<=1;irec=irec+RECINC){

        /* construct Hessian for different material parameters */
            for (i=1;i<=NX;i=i+IDX){
                for (j=1;j<=NY;j=j+IDY){
            
                    /* assemble complex wavefields */
                    utty = (utty_shot[j][i] + uttyi_shot[j][i] * I);
                    
                    eyy = (eyy_shot[j][i] + eyyi_shot[j][i] * I);
                    eyx = (eyx_shot[j][i] + eyxi_shot[j][i] * I);
                    
                    
                    if(INVMAT1==1){
                       muss = prho[j][i] * pu[j][i] * pu[j][i];
                    }
                    
                    if(INVMAT1=3){
                       muss = pu[j][i];
                    }
                      
                    /* Hessian */  
                    
                    tmp_jac_rho = (conj(utty)*utty);                   
                    tmp_jac_mu = (conj(eyx)*abs_green[j][i]*eyx) + (conj(eyy)*abs_green[j][i]*eyy);
                    
                    /* calculate Hessian for lambda, mu and rho by autocorrelation of Frechet derivatives */
                    /*if(INVMAT1==3){*/
                       hessian_u[j][i] += HESS_SCALE * creal(tmp_jac_mu);  
                     hessian_rho[j][i] += HESS_SCALE * creal(tmp_jac_rho);
                    /*}*/

                    /* Assemble Hessian for Vp, Vs and rho by autocorrelation of Frechet derivatives*/
                    /*if(INVMAT1==1){
                    
                         tmp_jac_vp = 2.0 * ppi[j][i] * prho[j][i] * tmp_jac_lam;          
                         tmp_jac_vs = (- 4.0 * prho[j][i] * pu[j][i] * tmp_jac_lam) + (2.0 * prho[j][i] * pu[j][i] * tmp_jac_mu);                  
                         tmp_jac_rho += (((ppi[j][i] * ppi[j][i])-(2.0 * pu[j][i] * pu[j][i])) * tmp_jac_lam) + (pu[j][i] * pu[j][i] * tmp_jac_mu);
                    
                         hessian[j][i] += HESS_SCALE * creal(tmp_jac_vp*conj(tmp_jac_vp));  
                       hessian_u[j][i] += HESS_SCALE * creal(tmp_jac_vs*conj(tmp_jac_vs));  
                     hessian_rho[j][i] += HESS_SCALE * creal(tmp_jac_rho*conj(tmp_jac_rho));
                     
                    }*/
                  
                 }
             }
 }
}


/* apply wavenumber damping for Vp-, Vs- and density Hessian */
/*if(SPATFILTER==1){
    wavenumber(hessian_u);
    wavenumber(hessian_rho);
  }*/

/* save HESSIAN for mu */
/* ----------------------- */
sprintf(jac,"%s_hessian_u_%d.%i%i",JACOBIAN,iter,POS[1],POS[2]);
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
sprintf(jac,"%s_hessian_u_%d",JACOBIAN,iter);
if (MYID==0) mergemod(jac,3);

/* save HESSIAN for rho */   
/* ----------------------- */
sprintf(jac,"%s_hessian_rho_%d.%i%i",JACOBIAN,iter,POS[1],POS[2]);
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
sprintf(jac,"%s_hessian_rho_%d",JACOBIAN,iter);
if (MYID==0) mergemod(jac,3);

free_matrix(hessian_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_mu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_rho,-nd+1,NY+nd,-nd+1,NX+nd);

free_matrix(eyy_shot,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(eyx_shot,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(utty_shot,-nd+1,NY+nd,-nd+1,NX+nd);

free_matrix(eyyi_shot,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(eyxi_shot,-nd+1,NY+nd,-nd+1,NX+nd);  
free_matrix(uttyi_shot,-nd+1,NY+nd,-nd+1,NX+nd);                

}