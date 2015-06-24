
/* Normalize data 
*
*  Daniel Koehn
*  Kiel, the 5th of October 2013 
*/


#include "fd.h"

void normalize_data_rms_MPI(float **data, float **data_obs, int ntr, int ns){

	float rms, rms_obs, rms1, rms_obs1;
	int i,j;
 
        rms=0.0;
      	rms_obs=0.0;

	for(i=1;i<=ntr;i++){        
     	    for(j=2;j<=ns;j++){

      	       /* Estimate rms-values of model and field data */
               rms += fabs(data[i][j]);
               rms_obs += fabs(data_obs[i][j]);
               
            }

	}

        rms=rms/(ntr*(ns-1));
        rms_obs=rms_obs/(ntr*(ns-1));

        /* set rms-values=0.0 to 1.0 */
	if(rms==0.0){rms=1.0;}
        if(rms_obs==0.0){rms_obs=1.0;}

        /*MPI_Allreduce(&rms,&rms1,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&rms_obs,&rms_obs1,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);*/

	for(i=1;i<=ntr;i++){        

	    for(j=2;j<=ns;j++){
		data_obs[i][j] = data_obs[i][j]*(rms1/rms_obs1);
	    }

	}

}
