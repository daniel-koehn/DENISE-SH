/*------------------------------------------------------------------------
 *   Inversion of source time function 
 *
 *   Daniel Koehn
 *   Kiel, the 12th of September 2013
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <fftw3.h>

void stf(float **sectionvy_obs, float **sectionvy, int ntr_glob, int ishot, int ns, int iter, int nshots, float **signals, int **recpos, float **srcpos){ 

	/* declaration of global variables */
	extern float DT, DH, OFFSETC, EPS_STF;
	extern int SEIS_FORMAT, MYID, NT, NORMALIZE, TIMEWIN;
	extern char  SEIS_FILE_VY[STRING_SIZE], PARA[STRING_SIZE], DATA_DIR[STRING_SIZE];
	extern int TRKILL, OFFSET_MUTE;
	extern char TRKILL_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
	
	/* declaration of variables for trace killing */
        float *STF_vector;
	int ** kill_tmp, *kill_vector, h, j, i, Npad;
	char trace_kill_file[STRING_SIZE];	
        double npadd;
	FILE *ftracekill, *STF;

        /* declaration of variables for offset-muting */
        float offset, xr, yr, xs, ys;

        /* complex variables for source wavelet estimation */
        fftw_complex *sumn, *sumd, *D_s, *D_ss, *D_ss_fd, *D_s_td;
        float maxdr, maxdi;
        char signal_wave[STRING_SIZE];

        sprintf(signal_wave,"%s_shot_%i.dat",SIGNAL_FILE,ishot);
	
	printf("\n================================================================================================\n\n");
	printf("\n ***** Inversion of Source Time Function - shot: %d - it: %d ***** \n\n",ishot,iter);

        /*Npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );*/
        Npad = ns;
    
        /* Allocate memory */
        sumn  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        sumd  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_s  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_ss  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_ss_fd  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
        D_s_td  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

        STF_vector = vector(1,ns);

        /*printf("Npad=%d \n",Npad);   */
        
        /* apply time window */                                            
        if(TIMEWIN==1){
          /*time_window_stf(sectionvy, iter, ntr_glob, ns, ishot);*/
          time_window_stf(sectionvy_obs, iter, ntr_glob, ns, ishot);
        }
                               
	/* TRKILL==1 - trace killing is applied */
	if(TRKILL){
	  kill_tmp = imatrix(1,nshots,1,ntr_glob);
	  kill_vector = ivector(1,ntr_glob);

	  ftracekill=fopen(TRKILL_FILE,"r");

	  if (ftracekill==NULL) err(" Trace kill file could not be opened!");

		for(i=1;i<=nshots;i++){
			for(j=1;j<=ntr_glob;j++){
				fscanf(ftracekill,"%d",&kill_tmp[i][j]);
			}
		}

		fclose(ftracekill);

		for(i=1;i<=ntr_glob;i++){
	   	   kill_vector[i] = kill_tmp[ishot][i];
		}

	}
	
	if(TRKILL){
	  for(i=1;i<=ntr_glob;i++){

	     if(kill_vector[i]==1){
	       
	        for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }
	     }	
    	     
	  }
	
	}	
	/* trace killing ends here */

        /* apply offset mute */
        /*if(OFFSET_MUTE){

         printf("OFFSETC = %f \n",OFFSETC);
         printf("OFFSET_MUTE = %d \n",OFFSET_MUTE);      

         for (i=1;i<=ntr_glob;i++){
             
             /* calculate source and receiver positions */
      	    /* xr = recpos[1][i]*DH;
      	     xs = srcpos[1][ishot];
      	     yr = recpos[2][i]*DH;
      	     ys = srcpos[2][ishot];

             /* calculate absolute offset */
             /*offset = sqrt(((xs-xr)*(xs-xr))+((ys-yr)*(ys-yr)));

             printf("offset = %f \n",offset);

             /* mute far-offset data*/
             /*if((OFFSET_MUTE==1)&&(offset>=OFFSETC)){
                
                for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }
                    
             } */

             /* mute near-offset data*/
             /*if((OFFSET_MUTE==2)&&(offset<=OFFSETC)){

	        for(j=1;j<=ns;j++){
		   sectionvy[i][j]=0.0;
		   sectionvy_obs[i][j]=0.0;
	        }

             }*/
        /* } 

        } /* end of OFFSET_MUTE */	

        /* Trace normalization to maximum amplitude of each trace */
        if(NORMALIZE==1){
          normalize_data(sectionvy,ntr_glob,ns);
          normalize_data(sectionvy_obs,ntr_glob,ns);
        }
        
        /* Trace normalization of field data with respect to maximum amplitude of model data */
        if(NORMALIZE==2){
          normalize_data_rel(sectionvy,sectionvy_obs,ntr_glob,ns);
        }

        /* initialize nominator and denominator for Wiener deconvolution */	
        for(j=0;j<Npad;j++){

           sumn[j][0]=0.0;
           sumd[j][0]=0.0;

           sumn[j][1]=0.0;
           sumd[j][1]=0.0;

        }

        /* FFT spike wavelet */
        for(j=0;j<Npad;j++){

           if(j<ns){
	     D_ss[j][0]=signals[1][j+1];
           }
           else{D_ss[j][0]=0.0;}

           D_ss[j][1]=0.0;           

        }

        fftw_plan p_s;

        p_s = fftw_plan_dft_1d(Npad, D_ss, D_ss_fd, 1, FFTW_ESTIMATE);         
        fftw_execute(p_s);

        /* FFT of each data and model trace and calculation of nominator and denominator of Wiener deconvolution */
        for(i=1;i<=ntr_glob;i++){

           /* allocate memory for complex variables */           
           fftw_complex *in_data, *out_data, *in_model, *out_model;
           fftw_plan p_data,p_model;
         
           in_data  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
           out_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

           in_model  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);
           out_model = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Npad);

           /* define real and imaginary parts of data vectors and apply zero-padding */
           for(j=0;j<Npad;j++){


              if(j<ns){
	      	 in_model[j][0]=sectionvy[i][j+1];
                 in_data[j][0]=sectionvy_obs[i][j+1];	
              }
              else{
	      	in_model[j][0]=0.0;
                 in_data[j][0]=0.0;
              }

              in_model[j][1]=0.0;
              in_data[j][1]=0.0; 
              
	   }
           
           /* apply FFTW */           
           p_data  = fftw_plan_dft_1d(Npad, in_data, out_data, 1, FFTW_ESTIMATE);
           p_model = fftw_plan_dft_1d(Npad, in_model, out_model, 1, FFTW_ESTIMATE);         

           fftw_execute(p_data);
           fftw_execute(p_model);

           /* estimate nominator and denominator of the Wiener deconvolution */
           for(j=0;j<Npad;j++){

              /* real parts of the nominator and denominator */
              sumn[j][0] += (out_data[j][0]*out_model[j][0]+out_data[j][1]*out_model[j][1]);
              sumd[j][0] += (out_model[j][0]*out_model[j][0]+out_model[j][1]*out_model[j][1]);

              /* imaginary parts of the nominator and denominator */
              sumn[j][1] += (out_data[j][0]*out_model[j][1]-out_data[j][1]*out_model[j][0]);
              sumd[j][1] += (out_model[j][0]*out_model[j][0]+out_model[j][1]*out_model[j][1]);

           }           
 
           fftw_destroy_plan(p_data);
           fftw_free(in_data); 
           fftw_free(out_data);          

           fftw_destroy_plan(p_model);
           fftw_free(in_model); 
           fftw_free(out_model);

        }

        /* estimate maximum of real and imaginary part of the denominator*/
        maxdr = 0.0;
        maxdi = 0.0;
        for(j=0;j<Npad;j++){
        
           if(fabs(sumd[j][0])>maxdr){
             maxdr = fabs(sumd[j][0]);
           } 
           
           if(fabs(sumd[j][1])>maxdi){
             maxdi = fabs(sumd[j][1]);
           }
                                              
        }

        /* construct source wavelet in frequency domain by Wiener deconvolution */        		
        for(j=0;j<Npad;j++){
           
           /*D_s[j][0] = D_ss_fd[j][0]*(sumn[j][0]/(sumd[j][0]+EPS_STF*maxdr));
           D_s[j][1] = D_ss_fd[j][1]*(sumn[j][1]/(sumd[j][1]+EPS_STF*maxdi));*/

           D_s[j][0] = sumn[j][0]/(sumd[j][0]+EPS_STF*maxdr);
           D_s[j][1] = sumn[j][1]/(sumd[j][1]+EPS_STF*maxdi);

           /*D_s[j][0] = sumn[j][0]/sumd[j][0];
           D_s[j][1] = -sumn[j][1]/sumd[j][1];*/

                      
        }

        /* inverse FFTW of the estimated STF */
        fftw_plan p_stf;
        p_stf  = fftw_plan_dft_1d(Npad, D_s, D_s_td, -1, FFTW_ESTIMATE);
        fftw_execute(p_stf);

        /* extract real part */
        for(j=0;j<ns;j++){
           STF_vector[j+1]=D_s_td[j+1][0];
        }	

        /* normalization of source wavelet to maximum amplitude for better numerical performance */
        if(NORMALIZE){
          normalize_STF(STF_vector,ns);
        }

        /* output of the STF */
        STF=fopen(signal_wave,"w");                                	
        for(j=1;j<=ns;j++){
           fprintf(STF,"%e\n",STF_vector[j]);
        }	
        fclose(STF);
						
	printf("\n\n================================================================================================\n");
	
	
	/* free memory for trace killing and FFTW */
	if(TRKILL){
	free_imatrix(kill_tmp,1,nshots,1,ntr_glob);
	free_ivector(kill_vector,1,ntr_glob);
	}

        fftw_free(sumn); 
        fftw_free(sumd); 
	fftw_free(D_s);
 	fftw_free(D_ss);
	fftw_free(D_ss_fd);
	fftw_free(D_s_td);
        fftw_destroy_plan(p_stf);
        fftw_destroy_plan(p_s);

        free_vector(STF_vector,1,ns);
}

