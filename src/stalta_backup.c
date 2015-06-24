/*
 *------------------------------------------------------------------------
 *
 *   STA/LTA first arrival time picker
 *
 *------------------------------------------------------------------------
*/


#include "fd.h"

/* subroutine prototypes */
float mean(float *vec, int istart, int iend);
float maximumpos(float *a, int n);

void stalta(float **sectionp, int ntr, int nst, float *picked_times, int ishot)
{
	
	extern float DT;
	extern int POS[3];

	int	stawin;		/* sta time window in samples			*/
	int	ltawin;		/* lta time window in samples			*/
	float	*sta;		/* average of samples over window stawin	*/
	float	*lta;		/* average of samples over window ltawin	*/
	float	*stalta;	/* sta/lta					*/
	int	staltawin;	/* search max(sta/lta) in (1,...,staltawin)	*/
	int	maxpos=1;	/* position of maximum in stalta		*/
	float	*stmp, *ltmp;	/* temporary variables			*/
	int	j,k,i,ntshift;
	float	*data, *data1, maxdata, maxstalta, tshift, tshift1;
	char    jac[225];
        /*FILE    *fp1;*/

	/* define time windows */
        stawin    = ceil(0.007/DT); /*10*/
	ltawin    = ceil(0.08/DT);  /*30*/
	staltawin = nst; /*500*/
	tshift = 0.0;  /* 0.2 */  /* timeshift [s] to avoid picking errors at small offsets */
	tshift1 = 0.0; /* 0.2 */  /* arbitarily timeshift [s] for picked arrival times */
	ntshift = floor(tshift/DT);  /* timeshift [timesamples] */
	
	/* allocate memory for vector variables */
	lta = vector(1,nst);
	sta = vector(1,nst);
	ltmp = vector(1,ltawin);
	stmp = vector(1,stawin);
	stalta = vector(1,nst);
	data = vector(1,nst);
        data1 = vector(1,nst);

        /*sprintf(jac,"%s_shot%d.%d%d",PICKS_FILE,ishot,POS[1],POS[2]);
        fp1=fopen(jac,"w");*/

	/* Main loop over traces */
	for(i=1;i<=ntr;i++) { 

		/*****************/
		/* extract trace */
		/*****************/
		for (j=1;j<=nst;j++) {
			if (j==1) {
				data[1] = 0.0;
				continue;
			}
			data[j] = sectionp[i][j];
			if (j==2)  maxdata = data[j];
			else { if (maxdata<data[j]) maxdata = data[j]; }
		}
		
		/****************************/
		/* Apply time shift to trace*/
		/****************************/
		for (j=1;j<=nst;j++) {
		    
		    if((j+ntshift)<=nst){data1[j+ntshift]=data[j];}
		    
		}
				
		/***********/
		/* STA/LTA */
		/***********/
		/* square of trace data */
		for (j=1;j<=nst;j++)  data[j] = data1[j]*data1[j]/(maxdata*maxdata);
		
		/* apply STA/LTA-algorithm to trace */
		for (j=1;j<=nst-ltawin;j++) {
			/* copy trace piece to temporary variable */
			for (k=0;k<=ltawin-1;k++) ltmp[k] = data[j+k];
			/* average */			
			lta[j+ltawin-1] = mean(ltmp, 1, ltawin);
		}
		for (j=1;j<=nst-stawin;j++) {
			/* copy trace piece to temporary variable */
			for (k=0;k<=stawin-1;k++) stmp[k] = data[j+k];
			/* average */			
			sta[j+stawin-1] = mean(stmp, 1, stawin);
		}
		for (j=1;j<=nst;j++) {
			if (lta[j] == 0 || sta[j] == 0)
				stalta[j] = -1; /*setnan(stalta[j]);*/
			else
				stalta[j] = sta[j]/lta[j];	

		}
		maxpos = maximumpos(stalta, staltawin);
                maxpos = maxpos - (int)(((float)(ltawin-stawin))/2);
		picked_times[i] = (float)maxpos*DT-tshift+tshift1;
	}	
		/* output of picked times to file */
		/*fprintf(fp1,"%d \t %e \n",i,picked_times[i]);*/
	
	/*fclose(fp1);*/
	
	/*free_vector(lta,1,nst);
	free_vector(sta,1,nst);
	free_vector(ltmp,1,ltawin);
	free_vector(stmp,1,stawin);
	free_vector(stalta,1,nst);
	free_vector(data,1,nst);
        free_vector(data1,1,nst);*/

}


/*========================================================================*/
/*                              SUBROUTINES                               */
/*========================================================================*/

float mean(float *vec, int istart, int iend)
{
	int i;
	float m=0;
	
	for (i=istart;i<=iend;i++) m += vec[i];
	m /= (iend-istart+1);
	
	return m;
}


float maximumpos(float *a, int n){
	float maxi=0.0;
	int j;
	int maxpos = 1;

	for (j=1;j<=n;j++)
	{
		/* if (isnan(a[j])) continue; */
		/* if (!finite(a[j])) continue; */
		if (a[j]<0) continue;
		if (fabs(a[j])>maxi)
		{
			maxi = fabs(a[j]);
			maxpos = j;
		}
	}
		
	return maxpos;
}
