/*------------------------------------------------------------------------
 * Module for the Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
 * method for the elastic multiparameter inversion of
 * vp, vs, rho and lambda, mu, rho respectively
 * 
 * Daniel Koehn
 * ----------------------------------------------------------------------*/

#include "fd.h"

void LBFGS1(float ** waveconv, float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float C_vp, float ** gradp, int nfstart_jac,
	     float ** waveconv_u, float C_vs, float ** gradp_u, float ** waveconv_rho, float C_rho, float ** gradp_rho, float *** y_LBFGS_vp, float *** s_LBFGS_vp, float * rho_LBFGS_vp, float * rho_LBFGS_vs, 
		 float * alpha_LBFGS_vp, float * alpha_LBFGS_vs, float *** y_LBFGS_vs, float *** s_LBFGS_vs, float *** y_LBFGS_rho, float *** s_LBFGS_rho, float ** ppi, float ** pu, float ** prho, int nxnyi, float ** qvs, 
		 float ** rvs, float * beta_LBFGS_vs, float eps_true){

	extern int NX, NY, IDX, IDY, SPATFILTER;
	extern int HESSIAN, INVMAT, SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, k;
	float betaz, betan, gradplastiter, gradclastiter, betar, beta;
	float gamma_LBFGS_vp, gamma_LBFGS_vs, gamma_LBFGS_rho, Vp_sum, Vs_sum;
        float LBFGSTMP, LBFGSTMP1, LBFGSTMP2, LBFGSTMP3, modellastiter, norm_fac, norm_fac_u, norm_fac_rho;
        int NLBFGS, ki, itershift, iter1;
	extern FILE *FP;
	FILE *FP3, *FP4, *FP6, *FP5, *FP7;
	
	NLBFGS = 200;
        itershift = 1;

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT ZP ================================================================================== */
/* ===================================================================================================================================================== */


if((INVMAT==1)||(INVMAT==0)){

/* Normalize gradient to maximum value */
/*norm_fac=norm(waveconv,iter,1);*/

/* Normalization of the gradient */
/* ------------------------------- */
        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
                 
                 waveconv[j][i] = C_vp*waveconv[j][i];
                  
           }
        }

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
   taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,4);
}


for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
	    gradp[j][i] = waveconv[j][i];
    }
}

/* apply spatial wavelength filter */
if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
        spat_filt(waveconv,iter,1);
}


/* save gradient for output as inversion result */
    if(iter==nfstart_jac){
	sprintf(jac,"%s_p_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
	    fwrite(&waveconv[j][i],sizeof(float),1,FP3);
        }
    }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Zs ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==3)||(INVMAT==0)){
	
/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_u[j][i] = C_vs * waveconv_u[j][i];
   }
}

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,5);
}

/* Normalize gradient to maximum value */
/*norm_fac_u=norm(waveconv_u,iter,2);
if(MYID==0){printf("norm_fac_u=%e \n",norm_fac_u);}*/

for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp_u[j][i] = waveconv_u[j][i];
   }
}

/* apply spatial wavelength filter */
if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
        spat_filt(waveconv_u,iter,2);}

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_u_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_u_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}

}

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT rho ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==2)||(INVMAT==0)){

/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_rho[j][i] = C_rho * waveconv_rho[j][i];
   }
}

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,6);
}

/* Normalize gradient to maximum value */
/*norm_fac_rho=norm(waveconv_rho,iter,3);
if(MYID==0){printf("norm_fac_rho=%e \n",norm_fac_rho);}*/

for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp_rho[j][i] = waveconv_rho[j][i];
   }
} 

/* apply spatial wavelength filter */
if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv_rho,iter,3);}

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_rho_it%d.old.%i%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}
}

/* calculate H^-1 * waveconv, using the L-BFGS method, if iter > 1 */
/* --------------------------------------------------------------------- */

if(iter>1){

   
   /* load old models and gradients - Vp */
   /* ---------------------------------- */

   sprintf(jac,"%s_p.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_vp.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");

   iter1 = iter-itershift; /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */

   LBFGSTMP = 0.0;
   
   /*printf("k = %d \t iter1 = %d \t MYID = %d \n",k,iter1, MYID);*/
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS_vp[j][i][iter1] = waveconv[j][i]-gradplastiter;

	  fread(&modellastiter,sizeof(float),1,FP7);
          /*s_LBFGS_vp[j][i][iter1] = -eps_true*gradplastiter;*/
          s_LBFGS_vp[j][i][iter1] = ppi[j][i]-modellastiter; 

       }
     }
     
     fclose(FP6);
     fclose(FP7);

   /* load old models and gradients - Rho */
   /* ----------------------------------- */

   sprintf(jac,"%s_p_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");
   
   sprintf(jac1,"%s_p_mrho.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");

   iter1 = iter-itershift; /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */
   
   /*printf("k = %d \t iter1 = %d \t MYID = %d \n",k,iter1, MYID);*/
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS_rho[j][i][iter1] = waveconv_rho[j][i]-gradplastiter;

	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS_rho[j][i][iter1] = -eps_true*gradplastiter;

       }
     }
     
     fclose(FP6);
     fclose(FP7);
   
   /* load old models and gradients - Vs */
   sprintf(jac,"%s_p_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_vs.old.%i%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");
   
   /*printf("eps_true = %f \n",eps_true);*/
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS_vs[j][i][iter1] = waveconv_u[j][i]-gradplastiter;

    	  fread(&modellastiter,sizeof(float),1,FP7);
          /*s_LBFGS_vs[j][i][iter1] = -eps_true*gradplastiter;*/
          s_LBFGS_vs[j][i][iter1] = pu[j][i]-modellastiter;  

       }
     }
     
     fclose(FP6);
     fclose(FP7);
     
     /* calculate improved first guess Hessian */
     LBFGSTMP = 0.0; 
     LBFGSTMP1 = 0.0;

     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate nominator and denominator of gamma */
		  LBFGSTMP += y_LBFGS_vs[j][i][iter1] * s_LBFGS_vs[j][i][iter1];
		  LBFGSTMP1 += y_LBFGS_vs[j][i][iter1] * y_LBFGS_vs[j][i][iter1];
	 
       }
     }

	 /* Sum nominator and denominator of gamma of all CPUs */
         Vp_sum = 0.0;
         MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     
	 Vs_sum = 0.0;
         MPI_Allreduce(&LBFGSTMP1,&Vs_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	 gamma_LBFGS_vs = Vp_sum/Vs_sum;
     
     /*printf("gamma_LBFGS_vp = %e \n",gamma_LBFGS_vp);*/
         
	 /* calculate rho */
	 for(k=iter1;k>=1;k--){
       
           LBFGSTMP = 0.0;
	   /* calculate rho at iteration step iter for all parameter classes*/
	   for (i=1;i<=NX;i=i+IDX){
             for (j=1;j<=NY;j=j+IDY){
          
	          LBFGSTMP += y_LBFGS_vs[j][i][k] * s_LBFGS_vs[j][i][k];

              }
	   }
	   
           /* Sum over Rho of all CPUs */
	   Vp_sum = 0.0;
	   MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	   rho_LBFGS_vs[k] = 1.0/Vp_sum;
	   
	   /*if(MYID==0){                                                
	   printf("rho_LBFGS_vp = %e of k = %d \n",rho_LBFGS_vp[k],k);}*/
	                                                       
	   
	 }
         
         /* save q for all material parameters */
	 for (i=1;i<=NX;i=i+IDX){
          for (j=1;j<=NY;j=j+IDY){
          
		  /*qvs[j][i][iter1] = waveconv[j][i];*/
		  qvs[j][i] = waveconv_u[j][i];
		   
          }
	 }

     /* update alpha and q*/
     for(k=iter1;k>=1;k--){
		
		LBFGSTMP = 0.0;
	    
          for (i=1;i<=NX;i=i+IDX){
             for (j=1;j<=NY;j=j+IDY){
          
		  /* calculate s'q */
		  LBFGSTMP += s_LBFGS_vs[j][i][k]*qvs[j][i];
		     
             }
	  }

     /* Sum over alpha of all CPUs */
     Vp_sum = 0.0;
     MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     alpha_LBFGS_vs[k] = rho_LBFGS_vs[k] * Vp_sum;

	 /* update q for all material parameters */
      for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
          
		  qvs[j][i] = qvs[j][i] - alpha_LBFGS_vs[k] * y_LBFGS_vs[j][i][k]; 

         }
      }

     }
	 
	 /* Multiply gradient with approximated Hessian */
	 for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	           
	           rvs[j][i] = gamma_LBFGS_vs * qvs[j][i];
	           /*rvs[j][i][1] = taper_coeff[j][i] * qvs[j][i][1];*/
		     
		   }
	 }

	 /* calculate H^-1 * waveconv[j][i] */
	 for(k=1;k<=iter1;k++){

	   LBFGSTMP = 0.0;
	   for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               
		  /* calculate beta */
		  LBFGSTMP += y_LBFGS_vs[j][i][k]*rvs[j][i];

           }
	 }

	 /* Sum over beta of all CPUs */
         Vp_sum = 0.0;
         MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
         beta_LBFGS_vs[k] = rho_LBFGS_vs[k] * Vp_sum;

         for (i=1;i<=NX;i=i+IDX){
	    for (j=1;j<=NY;j=j+IDY){
                               
	        rvs[j][i] = rvs[j][i] + s_LBFGS_vs[j][i][k]*(alpha_LBFGS_vs[k]-beta_LBFGS_vs[k]);
            }
         }
         
	 }

	 /* update gradient */
	  for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               
		      waveconv_u[j][i] = rvs[j][i];
		  
           }
	  }

	 /* Denormalize Gradients */
	 for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
               
		   waveconv[j][i] = waveconv[j][i] * C_vp;  
		   waveconv_u[j][i] = waveconv_u[j][i];
		  /*waveconv_rho[j][i] = waveconv_rho[j][i] * norm_fac_rho;*/  /*XX no density update !!! */

            }
	 }

}

}

    /* save old model - Vp */
	sprintf(jac,"%s_p_vp.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&ppi[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_vp.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
                         



/* save old models Vs */
/* ------------------ */

    /* save old model */
	sprintf(jac,"%s_p_vs.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&pu[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_vs.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_u[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c_u.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);


/* save old models Rho */
/* ------------------ */

	sprintf(jac,"%s_p_mrho.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&prho[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_mrho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_rho[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g_rho */
        sprintf(jac,"%s_c_rho.old.%i%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
}
