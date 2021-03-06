/**/
/*------------------------------------------------------------------------
 *   Read FWI-workflow-parameters from Stdin                          
 *
 *  D. Koehn
 *  Kiel, the 8th of August 2013
 *  ----------------------------------------------------------------------*/

/* reading FWI-workflow parameters for DENISE */

#include "fd.h"

void read_par_inv(FILE *fp,int nstage,int stagemax){

/* declaration of extern variables */
extern int MYID;
extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER, NORMALIZE;
extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
extern int TIME_FILT, ORDER, EPRECOND;
extern int LNORM, OFFSET_MUTE;
extern int INV_STF;
extern float PRO, FC_START, FC_END, EPS_STF, OFFSETC, OFFSETC_STF; 
extern int TIMEWIN;
extern float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
extern float WD_DAMP, WD_DAMP1, SCALERHO;

/* definition of local variables */
int i;
char str [80];

fscanf(fp,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str,str);
for (i=1;i<=nstage;i++){
    fscanf(fp,"%f%i%f%f%i%i%f%f%f%i%i%i%i%f%f%i%i%i%f%f%i%i%f%f",&PRO,&TIME_FILT,&FC_START,&FC_END,&ORDER,&TIMEWIN,&GAMMA,&TWLENGTH_MINUS,&TWLENGTH_PLUS,&INV_VP_ITER,&INV_VS_ITER,&INV_RHO_ITER,&SPATFILTER,&WD_DAMP1,&WD_DAMP,&EPRECOND,&LNORM,&INV_STF,&OFFSETC_STF,&EPS_STF,&NORMALIZE,&OFFSET_MUTE,&OFFSETC,&SCALERHO);
}

fclose(fp);

if(MYID==0){
   printf("=========================================== \n");
   printf("       FWI-stage %d of %d \n",nstage,stagemax);
   printf("=========================================== \n");
   printf(" Density is inverted from iteration step %d.\n",INV_RHO_ITER);
   printf("\n");
   printf(" Vp is inverted from iteration step %d.\n",INV_VP_ITER);
   printf("\n");
   printf(" Vs is inverted from iteration step %d.\n",INV_VS_ITER);
   
   printf("\n\n");
   printf(" Smoothing (spatial filtering) of the gradients: \n ");
   if(SPATFILTER==1){
     printf(" \tSPATFILTER=%d: Gradients are smoothed.\n",SPATFILTER);
     printf(" \t(SPAT_FILT_SIZE=%d, SPAT_FILT_1=%d, SPAT_FILT_ITER=%d, WD_DAMP = %e, WD_DAMP1 = %e)\n",SPAT_FILT_SIZE,SPAT_FILT_1,SPAT_FILT_ITER,WD_DAMP,WD_DAMP1);}
     else 	printf(" \tSPATFILTER=%d: Gradients are not smoothed.\n",SPATFILTER);
   
   printf("\n\n");
   printf(" --------------- Frequency filtering -------------------\n");
	if (TIME_FILT){
		printf(" TIME_FILT=%d: Time domain filtering is applied \n",TIME_FILT);
		
		if(TIME_FILT==1){
		  printf(" Applying FWI with lowpass filter for corner frequency %f Hz\n",FC_END);
		  printf(" Order of lowpass filter is:\t%d\n\n",ORDER);
		}
		
                if(TIME_FILT==2){
                  printf(" Applying FWI with bandpass filter for corner frequencies %f Hz and %f Hz \n",FC_START,FC_END);
                  printf(" Order of bandpass filter at lower and upper frequencies are:\t%d\n\n",ORDER);
                }
	}	  
	else printf(" TIME_FILT=%d: No time domain filtering is applied.\n",TIME_FILT);

   printf("\n\n");
   printf(" --------------- Time windowing and damping -------------------\n");
	if (TIMEWIN){
		printf(" TIMEWIN=%d: Time windowing and damping is applied \n",TIMEWIN);
		/*printf(fp," Reading picked times from files: %s \n",PICKS_FILE);*/
		printf(" length of window after pick in s is: %f \n",TWLENGTH_PLUS);
		printf(" length of window befor pick in s is: %f \n",TWLENGTH_MINUS);
		printf(" gamma is : %f \n\n",GAMMA);}
	else printf(" TIMEWIN=%d: No time windowing and damping is applied \n",TIMEWIN);

   printf("\n\n");
   printf(" --------------- termination of the program -------------------\n");
   printf(" Misfit change during the last two iterations is smaller than %f percent.\n\n",(PRO*100.0));

   printf("\n\n");
   printf(" --------------- Energy preconditioning -------------------\n");
   if(EPRECOND==1){
      printf("EPRECOND = %d - Hessian approximation by Shin et al. (2001) \n\n",EPRECOND);
   }
   if(EPRECOND==3){
      printf("EPRECOND = %d - Hessian approximation by Plessix & Mulder (2004) \n\n",EPRECOND); 
   }
      else printf("EPRECOND = %d - energy preconditioning deactivated \n\n",EPRECOND);

   printf("\n\n");
   printf(" --------------- Gradient calculation -------------------\n");
   printf(" Misfit function:\n");
   printf("   LNORM==1 corresponds to L1 Norm\n");
   printf("   LNORM==2 corresponds to L2 Norm\n");
   printf("   LNORM==3 corresponds to Cauchy\n");
   printf("   LNORM==4 corresponds to SECH\n");
   printf("   LNORM==5 corresponds to global correlation\n\n");
   printf(" Used LNORM=%d\n\n",LNORM);
   
  printf("\n\n");
  if(INV_STF==1){
  printf(" --------------- STF inversion  -------------------\n");
  printf(" EPS_STF=%f\n\n",EPS_STF);
  }
  if(INV_STF==2){
  printf(" --------------- STF inversion using only near offset traces -------------------\n");
  printf(" EPS_STF=%f\n\n",EPS_STF);
  printf(" critical OFFSETC=%f\n\n",OFFSETC_STF);
  }
  if(INV_STF==0){printf(" --------------- No STF inversion  -------------------\n\n");}
  
  printf("\n\n");
  if(NORMALIZE){printf(" --------------- TRACE normalization  -------------------\n\n");
  printf(" NORMALIZE=%d\n\n",NORMALIZE);
  }
  
  printf("\n\n");
  if(OFFSET_MUTE==1){
    printf(" --------------- Apply far-offset mute  -------------------\n");
    printf(" OFFSET_MUTE=%d\n",OFFSET_MUTE);
    printf(" critical OFFSETC=%f\n\n",OFFSETC);
  }
  if(OFFSET_MUTE==2){
    printf(" --------------- Apply near-offset mute  -------------------\n");
    printf(" OFFSET_MUTE=%d\n",OFFSET_MUTE);
    printf(" critical OFFSETC=%f\n\n",OFFSETC);
  }
  
  if(OFFSET_MUTE==0){printf(" --------------- No Offset mute  -------------------\n\n");}
  printf(" --------------- Scale density gradient by factor SCALERHO = %f  -------------------\n\n",SCALERHO);
        
  
}

}
