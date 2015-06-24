/* Calculation and assembly of the Hessian in time domain*/
/*                                                               */
/* Daniel Koehn */
/* Kiel, the 16th of May 2013 */

#include "fd.h"
#include <fftw3.h>

void hessian(int NTDTINV, float *** hess_eyy, float *** hess_eyx, float *** hess_utty, float ** hessian_u, float ** hessian_rho){
  
extern float DT,DH,TIME;      	
extern int NX, NY, IDX, IDY, INVMAT1, MYID, POS[4];
extern char JACOBIAN[STRING_SIZE];

/* local variables */
int i, j, k, l, nt, ishot, irec, nd, NSRC_HESSIAN, NREC_HESSIAN, RECINC, npad;
char jac[STRING_SIZE];

/*if(MYID==0){printf("Bis hier gehts !");}*/

/* assemble Hessian */
/* ---------------- */

/*if(MYID==0){printf("Bis hier gehts 1!");}*/
        
        npad = (int)(pow(2.0, ceil(log((double)(NTDTINV))/log(2.0))+2.0) );  /* NTDTINV -> npad for usage in FFT */
        /*npad = NTDTINV;*/
        
       /* calculate spatial and temporal derivatives of the forward wavefield */ 
       for (i=1;i<=NX;i=i+IDX){
          for (j=1;j<=NY;j=j+IDY){

              /*fftw_complex *eyy, *eyx, *utty, *eyyfft, *eyxfft, *uttyfft;
              fftw_plan  p1, p2, p3;

	      eyy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	      eyx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
              utty = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);

              eyyfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
	      eyxfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);
              uttyfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npad);*/
              

              /* write data for each time column to vector */
              /*for (nt=1;nt<NTDTINV;nt++){

		  eyy[nt][0] = hess_eyy[j][i][nt];
	          eyx[nt][0] = hess_eyx[j][i][nt];
                  utty[nt][0] = hess_utty[j][i][nt];	

	          eyy[nt][1] = 0.0;
	          eyx[nt][1] = 0.0;
                  utty[nt][1] = 0.0;		

	      }
	      
              for (nt=NTDTINV;nt<npad;nt++){   
	                    
                  eyy[nt][0] = 0.0;
                  eyx[nt][0] = 0.0;
                  utty[nt][0] = 0.0;
	                                                                          
                  eyy[nt][1] = 0.0;
                  eyx[nt][1] = 0.0;
                  utty[nt][1] = 0.0;
	                                                                                                                                
              }*/
	        
	      /*if(MYID==0){printf("Bis hier gehts 21!");}*/
	                                                                                                                                              
                
              /* calculate FFTW */
	      /*p1  = fftw_plan_dft_1d(npad, eyy, eyyfft, FFTW_FORWARD, FFTW_ESTIMATE );
	      fftw_execute(p1);
	      
              p2  = fftw_plan_dft_1d(npad, eyx, eyxfft, FFTW_FORWARD, FFTW_ESTIMATE );
              fftw_execute(p2);
              
              p3  = fftw_plan_dft_1d(npad, utty, uttyfft, FFTW_FORWARD, FFTW_ESTIMATE );
              fftw_execute(p3);*/

              /*if(MYID==0){printf("Bis hier gehts 22!");} */

              /* assemble diagonal elements of Pseudo-Hessian */
              /*for (nt=1;nt<npad;nt++){
                hessian_u[j][i] += (eyyfft[nt][0]*eyyfft[nt][0]+eyyfft[nt][1]*eyyfft[nt][1]) + (eyxfft[nt][0]*eyxfft[nt][0]+eyxfft[nt][1]*eyxfft[nt][1]); 
                hessian_rho[j][i] += (uttyfft[nt][0]*uttyfft[nt][0]+uttyfft[nt][1]*uttyfft[nt][1]); 
              }*/
              
                for (nt=1;nt<NTDTINV;nt++){
                   hessian_u[j][i] += (hess_eyy[j][i][nt]*hess_eyy[j][i][nt]) + (hess_eyx[j][i][nt]*hess_eyx[j][i][nt]);
                   hessian_rho[j][i] += (hess_utty[j][i][nt]*hess_utty[j][i][nt]); 
                }
              
              /* clean up memory */
     
              /*fftw_free(eyy);
     	      fftw_free(eyx);
     	      fftw_free(utty);
     
     	      fftw_free(eyyfft);
     	      fftw_free(eyxfft);
     	      fftw_free(uttyfft);
     
              fftw_destroy_plan(p1);
              fftw_destroy_plan(p2);  
              fftw_destroy_plan(p3);*/

               
          }
       }
     
     /*if(MYID==0){printf("Bis hier gehts 3!");}*/

/* apply wavenumber damping for Vp-, Vs- and density Hessian */
/*if(SPATFILTER==1){
    wavenumber(hessian_u);
    wavenumber(hessian_rho);
  }*/                

}
