/* $Id: update_s_visc_PML.c,v 1.1.1.1 2011/10/06 22:44:52 groos Exp $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_visc_PML(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   uy, float **   uyx, float **   syz,
      float **   sxy, float ** ujp, float ** uip, float **rho, float *hc, int infoout,
      float ***r, float ***p, float ***q, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, 
      float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
      float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
      float ** psi_vyy, float ** psi_vyx){



	int i,j, m, fdoh, h, h1, l;
	float  vxx, vyy, vxy, vyx;
	float  dhi, dthalbe;	
	extern float DT, DH;
	extern int MYID, FDORDER, FW, L;
        extern int FREE_SURF, BOUNDARY, GRAD_FORM;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;
	
	float sumr=0.0, sump=0.0, sumq=0.0;
	
	

	/*dhi = DT/DH;*/
	dhi=1.0/DH;
	fdoh = FDORDER/2;
	dthalbe = DT/2.0;

	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_s (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}
	


	switch (FDORDER){

	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			
			   vyx = (  hc[1]*(vy[j][i+1]-vy[j][i]))*dhi;
                           vyy = (  hc[1]*(vy[j+1][i]  -vy[j][i]))*dhi; 
        
	                  if(FW>0){
        
        	             /* left boundary */                                         
        	             if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                               psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                               vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
		             }

		             /* right boundary */                                         
		             if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                                h1 = (i-nx2+2*FW);
                                h = i;
                        
                                psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			        vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
		              }

	                      /* top boundary */                                         
	                      if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                                psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                                vyy = vyy / K_y[j] + psi_vyy[j][i];

	                      }
	
	                      /* bottom boundary */                                         
	                      if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                                 h1 = (j-ny2+2*FW);                                        
                                 h = j;
                                                
                                 psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                                 vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		              }      
                          }
        
	
	                  /* computing sums of the old memory variables */
	                  sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     sumr+=r[j][i][l];
  			     sumq+=q[j][i][l];
		          }
			
                          /* updating components of the stress tensor, partially */
		          sxy[j][i] += (fipjp[j][i]*vyx)+(dthalbe*sumr);
		          syz[j][i] += (f[j][i]*vyy)+(dthalbe*sumq);
				
			
		          /* now updating the memory-variables and sum them up*/
		          sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*vyx));
			     q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(d[j][i][l]*vyy));
			     sumr += r[j][i][l];
		 	     sumq += q[j][i][l];
		          }
			
			  if(GRAD_FORM==2){
			     uy[j][i] = vyy;
		             uyx[j][i] = vyx;
			  }
			
		          /* and now the components of the stress tensor are completely updated */
		          sxy[j][i]+=(dthalbe*sumr);
		          syz[j][i]+=(dthalbe*sumq);        
			
		    }
		}
		break;

	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
                          
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
			       + hc[2]*(vy[j][i+2]-vy[j][i-1]))*dhi;

                        
		        vyy = (  hc[1]*(vy[j+1][i]  -vy[j][i])
			       + hc[2]*(vy[j+2][i]-vy[j-1][i]))*dhi; 
			  
	                  if(FW>0){
        
        	             /* left boundary */                                         
        	             if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                               psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                               vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
		             }

		             /* right boundary */                                         
		             if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                                h1 = (i-nx2+2*FW);
                                h = i;
                        
                                psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			        vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
		              }

	                      /* top boundary */                                         
	                      if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                                psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                                vyy = vyy / K_y[j] + psi_vyy[j][i];

	                      }
	
	                      /* bottom boundary */                                         
	                      if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                                 h1 = (j-ny2+2*FW);                                        
                                 h = j;
                                                
                                 psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                                 vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		              }      
                          }
        
	
	                  /* computing sums of the old memory variables */
	                  sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     sumr+=r[j][i][l];
  			     sumq+=q[j][i][l];
		          }
			
                          /* updating components of the stress tensor, partially */
		          sxy[j][i] += (fipjp[j][i]*vyx)+(dthalbe*sumr);
		          syz[j][i] += (f[j][i]*vyy)+(dthalbe*sumq);
				
			
		          /* now updating the memory-variables and sum them up*/
		          sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*vyx));
			     q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(d[j][i][l]*vyy));
			     sumr += r[j][i][l];
		 	     sumq += q[j][i][l];
		          }
			
			  if(GRAD_FORM==2){
			     uy[j][i] = vyy;
		             uyx[j][i] = vyx;
			  }
			
		          /* and now the components of the stress tensor are completely updated */
		          sxy[j][i]+=(dthalbe*sumr);
		          syz[j][i]+=(dthalbe*sumq);			

			}
		}
		break;

	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
			       + hc[2]*(vy[j][i+2]-vy[j][i-1])
			       + hc[3]*(vy[j][i+3]-vy[j][i-2]))*dhi;

                
                        vyy = (  hc[1]*(vy[j+1][i]  -vy[j][i])
			       + hc[2]*(vy[j+2][i]-vy[j-1][i])
			       + hc[3]*(vy[j+3][i]-vy[j-2][i]))*dhi; 

	                  if(FW>0){
        
        	             /* left boundary */                                         
        	             if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                               psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                               vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
		             }

		             /* right boundary */                                         
		             if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                                h1 = (i-nx2+2*FW);
                                h = i;
                        
                                psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			        vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
		              }

	                      /* top boundary */                                         
	                      if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                                psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                                vyy = vyy / K_y[j] + psi_vyy[j][i];

	                      }
	
	                      /* bottom boundary */                                         
	                      if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                                 h1 = (j-ny2+2*FW);                                        
                                 h = j;
                                                
                                 psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                                 vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		              }      
                          }
        
	
	                  /* computing sums of the old memory variables */
	                  sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     sumr+=r[j][i][l];
  			     sumq+=q[j][i][l];
		          }
			
                          /* updating components of the stress tensor, partially */
		          sxy[j][i] += (fipjp[j][i]*vyx)+(dthalbe*sumr);
		          syz[j][i] += (f[j][i]*vyy)+(dthalbe*sumq);
				
			
		          /* now updating the memory-variables and sum them up*/
		          sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*vyx));
			     q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(d[j][i][l]*vyy));
			     sumr += r[j][i][l];
		 	     sumq += q[j][i][l];
		          }
			
			  if(GRAD_FORM==2){
			     uy[j][i] = vyy;
		             uyx[j][i] = vyx;
			  }
			
		          /* and now the components of the stress tensor are completely updated */
		          sxy[j][i]+=(dthalbe*sumr);
		          syz[j][i]+=(dthalbe*sumq);
			

			}
		}
		break;

	case 8:

               for (j=ny1;j<=ny2;j++){
        	       for (i=nx1;i<=nx2;i++){
		
			vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
			       + hc[2]*(vy[j][i+2]-vy[j][i-1])
			       + hc[3]*(vy[j][i+3]-vy[j][i-2])
			       + hc[4]*(vy[j][i+4]-vy[j][i-3]))*dhi;

        
                        vyy = (  hc[1]*(vy[j+1][i]  -vy[j][i])
		               + hc[2]*(vy[j+2][i]-vy[j-1][i])
		               + hc[3]*(vy[j+3][i]-vy[j-2][i])
		               + hc[4]*(vy[j+4][i]-vy[j-3][i]))*dhi;		
		
	                  if(FW>0){
        
        	             /* left boundary */                                         
        	             if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                               psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * vyx;
                               vyx = vyx / K_x_half[i] + psi_vyx[j][i];                 
		             }

		             /* right boundary */                                         
		             if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                                h1 = (i-nx2+2*FW);
                                h = i;
                        
                                psi_vyx[j][h1] = b_x_half[h1] * psi_vyx[j][h1] + a_x_half[h1] * vyx;
			        vyx = vyx / K_x_half[h1] + psi_vyx[j][h1];
                                           
		              }

	                      /* top boundary */                                         
	                      if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                                                
                                psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
                                vyy = vyy / K_y[j] + psi_vyy[j][i];

	                      }
	
	                      /* bottom boundary */                                         
	                      if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                                 h1 = (j-ny2+2*FW);                                        
                                 h = j;
                                                
                                 psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                                 vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		              }      
                          }
        
	
	                  /* computing sums of the old memory variables */
	                  sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     sumr+=r[j][i][l];
  			     sumq+=q[j][i][l];
		          }
			
                          /* updating components of the stress tensor, partially */
		          sxy[j][i] += (fipjp[j][i]*vyx)+(dthalbe*sumr);
		          syz[j][i] += (f[j][i]*vyy)+(dthalbe*sumq);
				
			
		          /* now updating the memory-variables and sum them up*/
		          sumr=sumq=0.0;
		          for (l=1;l<=L;l++){
			     r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*vyx));
			     q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(d[j][i][l]*vyy));
			     sumr += r[j][i][l];
		 	     sumq += q[j][i][l];
		          }
			
			  if(GRAD_FORM==2){
			     uy[j][i] = vyy;
		             uyx[j][i] = vyx;
			  }
			
		          /* and now the components of the stress tensor are completely updated */
		          sxy[j][i]+=(dthalbe*sumr);
		          syz[j][i]+=(dthalbe*sumq);

                    }
                }
		break;

/*	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				
				
/*			}
		}
		break;
		
	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				
			}
		}
		break;
*/
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vyy = 0.0;
				vyx = 0.0;
				for (m=1; m<=fdoh; m++) {
					vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
					vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
				}
		
				vyy *= dhi;
				vyx *= dhi;

				sumr=sump=sumq=0.0;
				for (l=1;l<=L;l++){
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				sxy[j][i] += (fipjp[j][i]*vyx)+(dthalbe*sumr);
				syz[j][i] += (f[j][i]*vyy)+(dthalbe*sumq);

				sumr=sump=sumq=0.0;
				for (l=1;l<=L;l++){
					r[j][i][l] = bip[l]*(r[j][i][l]*cip[l]-(dip[j][i][l]*vyx));
					q[j][i][l] = bjm[l]*(q[j][i][l]*cjm[l]-(d[j][i][l]*vyy));
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i]+=(dthalbe*sumr);
				syz[j][i]+=(dthalbe*sumq);
			}
		}
		break;
		
	} /* end of switch(FDORDER) */


	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
