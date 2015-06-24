/* Calculation, assembly and output of the Hessian */
/*                                                 */
/* Daniel Koehn */
/* Kiel, the 3rd of september 2012 */

#include "fd.h"
#include <complex.h>

void hessian(int nshots, int SHOTINC, float *** green_vx, float *** greeni_vx, float *** green_vy, float *** greeni_vy, float *** green_sxx, float *** greeni_sxx, float *** green_syy, float *** greeni_syy,
             float *** green_sxy, float *** greeni_sxy, float ** prho, float ** pu, float ** ppi){
  
extern float DT,TIME;
extern float FC_HESSIAN;        	
extern int NX, NY, IDX, IDY, DTINV, INVMAT1, MYID, POS[4], FDORDER;
extern char JACOBIAN[STRING_SIZE];

/* local variables */
int i, j, k, l, ns_hess, ishot, irec, nd, NSRC_HESSIAN, RECINC;
double trig1,trig2;
double t=0.0;
const double pi=4.0*atan(1.0);
char jac[STRING_SIZE];

float complex green_x, green_y, tmp_mu1_shot, tmp_mu2_shot, tmp_mu3_shot, tmp_mu5_shot, tmp_mu6_shot;
float tmp_mu1, tmp_mu2, tmp_mu3;
float complex tmp_jac_lam, tmp_jac_mu, tmp_jac_rho, tmp_jac_vp, tmp_jac_vs, tmp_fft_fsignal;

float ** abs_green, omega, mulamratio, lamss, muss, HESS_SCALE;
float ** hessian, ** hessian_u, ** hessian_rho, ** hessian_lam, ** hessian_mu;

float ** hvxx_shot, ** hvxxi_shot, ** hvyy_shot, ** hvyyi_shot, ** hvxy_shot, ** hvxyi_shot,  ** hvyx_shot, ** hvyxi_shot;
float *psource_hess=NULL, *Hess_for_real=NULL, *Hess_for_complex=NULL;

FILE *FP4;

HESS_SCALE = 1e5;
RECINC = 1;
NSRC_HESSIAN=nshots;

nd = FDORDER/2 + 1;
abs_green = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

/* Diagonal elements of the Hessian*/
hessian = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hessian_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd); 
hessian_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

hvxx_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hvxxi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
hvyy_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hvyyi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
hvyx_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hvyxi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        
hvxy_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
hvxyi_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);

Hess_for_real = vector(1,1);
Hess_for_complex = vector(1,1);
        
for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
           hessian[j][i]=0.0;
           hessian_lam[j][i]=0.0;
           hessian_u[j][i]=0.0;  
           hessian_mu[j][i]=0.0; 
           hessian_rho[j][i]=0.0;
    }
}

/* assemble Hessian */
/* ----------------------------------------------------------------- */
/* calculate absolute values of impulse responses */
for (ishot=1;ishot<=nshots;ishot=ishot+SHOTINC){
    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
        
           green_x = green_vx[j][i][ishot] + greeni_vx[j][i][ishot] * I;
           green_y = green_vy[j][i][ishot] + greeni_vy[j][i][ishot] * I;
                                                                        
           /*abs_green[j][i] += creal((green_x*conj(green_x))+(green_y*conj(green_y)));*/
           abs_green[j][i] = 1.0;
        
           /*printf("green_x = %e \n",green_sxx[j][i][ishot]);*/
           
        }
    }
}    
     
omega = 2.0*M_PI*FC_HESSIAN;

/*printf("omega = %e \n",omega);
printf("NSRC = %d \n",NSRC_HESSIAN);*/

/*psource_hess=rd_sour(&ns_hess,fopen(SIGNAL_FILE,"r"));
FFT_data(psource_hess,Hess_for_real,Hess_for_complex,NT);
                                                         
MPI_Barrier(MPI_COMM_WORLD);                             

tmp_fft_fsignal = Hess_for_real[1] + Hess_for_complex[1] * I;*/

tmp_fft_fsignal = 1.0;

for (ishot=1;ishot<=nshots;ishot=ishot+SHOTINC){
        
        /*for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){*/
            
                /* calculate spatial derivatives of the forward wavefields */
                 /*hvxx_shot[j][i] = (green_vx[j][i][ishot]-green_vx[j][i-1][ishot])/DH;
                hvxxi_shot[j][i] = (greeni_vx[j][i][ishot]-greeni_vx[j][i-1][ishot])/DH;
                           
                 hvyy_shot[j][i] = (green_vy[j][i][ishot]-green_vy[j-1][i][ishot])/DH;
                hvyyi_shot[j][i] = (greeni_vy[j][i][ishot]-greeni_vy[j-1][i][ishot])/DH;
                           
                 hvyx_shot[j][i] = (green_vy[j][i+1][ishot]-green_vy[j][i][ishot])/DH; 
                hvyxi_shot[j][i] = (greeni_vy[j][i+1][ishot]-greeni_vy[j][i][ishot])/DH;
                           
                 hvxy_shot[j][i] = (green_vx[j+1][i][ishot]-green_vx[j][i][ishot])/DH;
                hvxyi_shot[j][i] = (greeni_vx[j+1][i][ishot]-greeni_vx[j][i][ishot])/DH;*/

        /*    }
        }*/

 for (irec=nshots;irec<=nshots;irec=irec+RECINC){

        /* construct Hessian for different material parameters */
            for (i=1;i<=NX;i=i+IDX){
                for (j=1;j<=NY;j=j+IDY){
            
                    /* Hessian for Lame parameter lambda, mu and rho */
                    tmp_mu1_shot = (green_sxx[j][i][ishot] + greeni_sxx[j][i][ishot] * I);
                    tmp_mu2_shot = (green_syy[j][i][ishot] + greeni_syy[j][i][ishot] * I);
                    tmp_mu3_shot = (green_sxy[j][i][ishot] + greeni_sxy[j][i][ishot] * I);
                    
                    tmp_mu5_shot = omega*((green_vx[j+1][i][ishot] * I) - greeni_vx[j+1][i][ishot]);
                    tmp_mu6_shot = omega*((green_vy[j+1][i][ishot] * I) - greeni_vy[j+1][i][ishot]);

                    tmp_mu1 = green_sxx[j][i][irec] + greeni_sxx[j][i][irec] * I;
                    tmp_mu2 = green_syy[j][i][irec] + greeni_syy[j][i][irec] * I;
                    tmp_mu3 = green_sxy[j][i][irec] + greeni_sxy[j][i][irec] * I;

                    if(INVMAT1==1){
                       muss = prho[j][i] * pu[j][i] * pu[j][i];
                      lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;}
                    
                    if(INVMAT1=3){
                       muss = pu[j][i];
                      lamss = ppi[j][i];}
                    
                    /*mulamratio = (muss * muss)/((lamss+muss)*(lamss+muss));*/
                                                                    
                    /*tmp_jac_lam = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * ((tmp_mu1_shot + tmp_mu2_shot) * (tmp_mu1_shot + tmp_mu2_shot));*/
                    
                    tmp_jac_lam =  (tmp_mu1_shot + tmp_mu2_shot) * (tmp_mu1 + tmp_mu2);
                    
                    tmp_jac_mu = ((1.0/(muss*muss))*(tmp_mu3_shot * tmp_mu3_shot)) + ((1.0/4.0) * ((tmp_mu1_shot + tmp_mu2_shot) * (tmp_mu1_shot + tmp_mu2_shot)) / ((lamss+muss)*(lamss+muss)))
                                                                                  + ((1.0/4.0) * ((tmp_mu1_shot - tmp_mu2_shot) * (tmp_mu1_shot - tmp_mu2_shot)) / (muss*muss));
                    
                    /*tmp_jac_mu = (tmp_mu3_shot * tmp_mu3) + (((mulamratio*((tmp_mu1_shot + tmp_mu2_shot) * (tmp_mu1 + tmp_mu2))) + ((tmp_mu1_shot - tmp_mu2_shot) * (tmp_mu1 - tmp_mu2)))/4.0);*/
                   
                    tmp_jac_rho = (tmp_mu5_shot*tmp_mu5_shot) + (tmp_mu6_shot*tmp_mu6_shot);
                    
                    /* Assemble Hessian for lambda, mu and rho */
                    if(INVMAT1==3){
                         hessian[j][i] += HESS_SCALE * creal(tmp_jac_lam*abs_green[j][i]*conj(tmp_jac_lam));
                       hessian_u[j][i] += HESS_SCALE * creal(tmp_jac_mu*abs_green[j][i]*conj(tmp_jac_mu));  
                     hessian_rho[j][i] += HESS_SCALE * creal(tmp_jac_rho*abs_green[j][i]*conj(tmp_jac_rho));
                    }

                    /* Assemble Hessian for Vp, Vs and rho*/
                    if(INVMAT1==1){
                    
                         tmp_jac_vp = 2.0 * ppi[j][i] * prho[j][i] * tmp_jac_lam;          
                         tmp_jac_vs = (- 4.0 * prho[j][i] * pu[j][i] * tmp_jac_lam) + (2.0 * prho[j][i] * pu[j][i] * tmp_jac_mu);                  
                         tmp_jac_rho += (((ppi[j][i] * ppi[j][i])-(2.0 * pu[j][i] * pu[j][i])) * tmp_jac_lam) + (pu[j][i] * pu[j][i] * tmp_jac_mu);
                    
                         hessian[j][i] += HESS_SCALE * creal(tmp_jac_vp*abs_green[j][i]*conj(tmp_jac_vp));  
                       hessian_u[j][i] += HESS_SCALE * creal(tmp_jac_vs*abs_green[j][i]*conj(tmp_jac_vs));  
                     hessian_rho[j][i] += HESS_SCALE * creal(tmp_jac_rho*abs_green[j][i]*conj(tmp_jac_rho));
                     
                    }
                  
                 }
             }
 }
}

/* save Hessian for Vp */
/* ----------------------- */
sprintf(jac,"%s_hessian.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){ 
        fwrite(&hessian[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);
    
MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */   
sprintf(jac,"%s_hessian",JACOBIAN);
if (MYID==0) mergemod(jac,3);

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

free_matrix(hessian,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_lam,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_mu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_rho,-nd+1,NY+nd,-nd+1,NX+nd);

}