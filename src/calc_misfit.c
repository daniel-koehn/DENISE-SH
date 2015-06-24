/*------------------------------------------------------------------------
 *   Calculate Misfit                                  
 *   last update 18/04/11, L. Rehor
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_misfit(float **sectiondata, float **section, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot,int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot){

/* declaration of variables */
extern float DT;
extern int MYID;
extern int TRKILL;
extern char TRKILL_FILE[STRING_SIZE];
float intseis, intseis_data, intseis_synthetics, abs_synthetics, abs_data;
int i,j;
float l2;
float L2_dummy;
int umax=0, h;
	

/* declaration of variables for trace killing */
int ** kill_tmp, *kill_vector;
char trace_kill_file[STRING_SIZE];	
FILE *ftracekill;

if(TRKILL){
	kill_tmp = imatrix(1,ntr_glob,1,nsrc_glob);
	kill_vector = ivector(1,ntr);

	ftracekill=fopen(TRKILL_FILE,"r");

	if (ftracekill==NULL) err(" Trace kill file could not be opened!");

	for(i=1;i<=ntr_glob;i++){
		for(j=1;j<=nsrc_glob;j++){
			fscanf(ftracekill,"%d",&kill_tmp[i][j]);
		}
	}

	fclose(ftracekill);

	h=1;
	for(i=1;i<=ntr;i++){
	   kill_vector[h] = kill_tmp[recpos_loc[3][i]][ishot];
	   h++;
	}
} /* end if(TRKILL)*/



/* calculate misfit */

for(i=1;i<=ntr;i++){
	if((TRKILL==1)&&(kill_vector[i]==1))
    	continue;	


    intseis = 0.0;
    intseis_data = 0.0;
    intseis_synthetics = 0.0;
    abs_data = 0.0;
    abs_synthetics = 0.0;
    L2_dummy=0.0;
      for(j=1;j<=ns;j++){
                        /*printf("%d \t %d \t %e \t %e \n",i,j,sectionpdata[i][j],sectionp[i][j]);*/
                        
			
			/* calculate L2 residuals */
			if((LNORM==2) ||(LNORM==6)){
			intseis += section[i][j]-sectiondata[i][j];}
			
			if(LNORM==5){
			intseis_data += sectiondata[i][j]*DT;
			intseis_synthetics += section[i][j]*DT;
			abs_data += intseis_data*intseis_data;
			abs_synthetics += intseis_synthetics*intseis_synthetics;
			L2_dummy+=(intseis_data*intseis_synthetics);}
			                        
			/* calculate norm */
			/*if((sws==1)&&(swstestshot==1)){*/
			if(((LNORM==2) ||(LNORM==6)) &&(swstestshot==1)){
			/*L2+=sectiondiff[i][invtime]*sectiondiff[i][invtime];*/
			L2+=intseis*intseis*DT*DT; 
			}
			
			}
			if(LNORM==5){
			L2 -= L2_dummy/(sqrt(abs_data)*sqrt(abs_synthetics));}

}
		    l2=L2;
return l2;

/* free memory for trace killing */
if(TRKILL){
free_imatrix(kill_tmp,1,ntr_glob,1,nsrc_glob);
free_ivector(kill_vector,1,ntr);
}

}
