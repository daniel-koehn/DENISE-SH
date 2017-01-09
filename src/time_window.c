/*------------------------------------------------------------------------
 *   Apply time damping (after Brossier (2009))                                 
 *   last update 04/01/2017, D.Koehn
 *   modified    02/02/2012, S.Heider
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void time_window(float **sectiondata, int iter, int ntr_glob, int ns, int ishot){

/* declaration of variables */
extern float DT;
extern int REC1, REC2, MYID, TIMEWIN;
extern int POS[3];
extern float GAMMA, TWLENGTH_PLUS, TWLENGTH_MINUS, DT;
extern char PICKS_FILE[STRING_SIZE];
int READ_PICKED_TIMES;
char pickfile_char[STRING_SIZE];
float time, dump, taper, taper1;
float *pick_tmp;
int i, j, h;

FILE *fptime;

pick_tmp = vector(1,ntr_glob);

/* read picked first arrival times from external pick file */
/* ------------------------------------------------------- */
if(TIMEWIN==1){   

   sprintf(pickfile_char,"%s%i.dat",PICKS_FILE,ishot);

   fptime=fopen(pickfile_char,"r");
   if (fptime == NULL) {
      err(" picks_?.dat could not be opened !");
   }

   for(i=1;i<=ntr_glob;i++){
      fscanf(fptime,"%f",&dump);
      pick_tmp[i] = dump + TWLENGTH_PLUS;    
   }

   fclose(fptime);
    
} /* end of if(TIMEWIN==1) */


/* Define constant time windows */
/* ---------------------------- */
if(TIMEWIN==2){
  h=1;
    for(i=1;i<=ntr_glob;i++){
    
        pick_tmp[h] = DT;
                
        h++;
    }                    
}

/* apply time damping */
for(i=1;i<=ntr_glob;i++){
      for(j=1;j<=ns;j++){
      
         if(TIMEWIN){time = (float)(j * DT);}

         dump = (time-pick_tmp[i]);
         taper = exp(-GAMMA*dump*dump);                  
	 
	 if(TIMEWIN){

           if(time>=pick_tmp[i]){
             sectiondata[i][j] = sectiondata[i][j] * taper;
           }                         
                          
	 }	 	 
                         	   
      }     
}

free_vector(pick_tmp,1,ntr_glob);

} /* end of function time_window.c */
