/*------------------------------------------------------------------------
 *   write seismograms to files 
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis(FILE *fp, float **sectionvy,  int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot, int ns, int iter){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID, RUN_MULTIPLE_SHOTS, INVMAT;	
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE];
	extern char  SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];

        char vyf[STRING_SIZE];
        int nsrc=1;		
	
	 /*if (RUN_MULTIPLE_SHOTS){*/
         
		sprintf(vyf,"%s.shot%d.%d",SEIS_FILE_VY,ishot,MYID);
        
		/*}
        else{
                sprintf(vxf,"%s.%d",SEIS_FILE_VX,MYID);
                sprintf(vyf,"%s.%d",SEIS_FILE_VY,MYID);
                sprintf(curlf,"%s.%d",SEIS_FILE_CURL,MYID);
                sprintf(divf,"%s.%d",SEIS_FILE_DIV,MYID);
                sprintf(pf,"%s.%d",SEIS_FILE_P,MYID);
                 
        }*/

	
	if(INVMAT==0){
	
	    if(iter>0){
	
	        sprintf(vyf,"%s.shot%d_it%d.%d",SEIS_FILE_VY,ishot,iter,MYID);
	
	    }
	
	    if(iter==-1){
	
	        sprintf(vyf,"%s.shot%d_it%d.%d",SEIS_FILE_VY,ishot,iter,MYID);
	
	    }
	}
	
	if(INVMAT==10){
	
	        sprintf(vyf,"%s.shot%d.%d",SEIS_FILE_VY,ishot,MYID);

	}
	
	
	switch (SEISMO){
	case 1 : /* particle velocity vy only */
		fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",MYID,ntr,vyf);
		outseis(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
		break;	
      }     
}
