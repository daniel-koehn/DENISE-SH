/* $Id: update_v_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen 
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_v_PML(int nx1, int nx2, int ny1, int ny2, int nt, float ** vy, float **  vyp1, float **  vym1, float **  utty, float ** syz,
	float ** sxy, float ** rho, float **  srcpos_loc, float ** signals1, int nsrc, float ** absorb_coeff,
	float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_syy_y, float ** psi_sxy_x){

	int i, j,l,fdoh,m, h, h1;
	float amp, dtdh, azi_rad;
	float vxtmp, vytmp, b;
        float syy_y, sxy_x;
	
	extern float DT, DH, ANGLE;
	double time1, time2;
	extern int MYID, QUELLTYP, QUELLTYPB, CHECKPTREAD, FDORDER;
        extern int FDORDER, INVMAT1, GRAD_FORM;
        extern int FREE_SURF, BOUNDARY, FW;
        extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;

	
	fdoh = FDORDER/2;
	dtdh = DT*DT/DH;
        
     /* drad = PI/180.0; 
        angle = 135.0; */
         
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_v (printed by PE %d):\n",MYID);
		fprintf(FP," Updating particle velocities ...");
	}


	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */

	switch (FDORDER){
	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1]);
                           syy_y = hc[1]*(syz[j][i]-syz[j-1][i]);
                          

        /* apply PML */
        if(FW>0){
        
		/* left boundary */                                         
		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                	           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                        	   sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
		}

		/* right boundary */                                         
		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                     
                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
		}

         
		/* top boundary */                                         
		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
			   psi_syy_y[j][i] = b_y[j] * psi_syy_y[j][i] + a_y[j] * syy_y;                                                
                           syy_y = syy_y / K_y[j] + psi_syy_y[j][i]; 
                        
                    
		}
		

		/* bottom boundary */                                         
		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
			    h = j;
		            		
                           
                           psi_syy_y[h1][i] = b_y[h1] * psi_syy_y[h1][i] + a_y[h1] * syy_y;                                                
			   syy_y = syy_y / K_y[h1] + psi_syy_y[h1][i]; 

		}
                              
        } /* end of PML */ 
                         
			   b = 1.0 / rho[j][i]; 
			   if(rho[j][i]<1e-20){b = 0.0;}			    
			 
                           if(sw==0){vyp1[j][i] = b * (sxy_x+syy_y)/DH;}


			   if(sw==1){
			   
			      if(GRAD_FORM==2){   /* FWI without data integration */
			        vyp1[j][i] = vy[j][i];
			      }
			      
			      if(GRAD_FORM==1){   /* FWI with data integration */
			        vyp1[j][i] += vy[j][i]*DT;
			      }
			   
			   }                                           
 
		           vy[j][i] += b * DT * (sxy_x+syy_y) / DH; 		         
      
     }
}
		break;
		
	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2]);

                           syy_y = hc[1]*(syz[j][i]-syz[j-1][i])
					    + hc[2]*(syz[j+1][i]-syz[j-2][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                           
                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		                   psi_syy_y[j][i] = b_y[j] * psi_syy_y[j][i] + a_y[j] * syy_y;                                                
                           syy_y = syy_y / K_y[j] + psi_syy_y[j][i]; 
                        
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
							h = j;
		            			   
       		               
                           psi_syy_y[h1][i] = b_y[h1] * psi_syy_y[h1][i] + a_y[h1] * syy_y;                                                
						   syy_y = syy_y / K_y[h1] + psi_syy_y[h1][i];
                           
        }                       
                           
			   b = 1.0 / rho[j][i]; 
			   if(rho[j][i]<1e-20){b = 0.0;}			    
			 
                           if(sw==0){vyp1[j][i] = b * (sxy_x+syy_y)/DH;}			                              

                           if(sw==1){
			   
			      if(GRAD_FORM==2){   /* FWI without data integration */
			        vyp1[j][i] = vy[j][i];
			      }
			      
			      if(GRAD_FORM==1){   /* FWI with data integration */
			        vyp1[j][i] += vy[j][i]*DT;
			      }
			   
			   }
                           
		           vy[j][i] += b * DT * (sxy_x+syy_y) / DH; 		         
      
     }
}
	     
		break;
		
	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3]);

                           syy_y = hc[1]*(syz[j][i]-syz[j-1][i])
					    + hc[2]*(syz[j+1][i]-syz[j-2][i])
					    + hc[3]*(syz[j+2][i]-syz[j-3][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                           h=i; 

                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
		                   psi_syy_y[j][i] = b_y[j] * psi_syy_y[j][i] + a_y[j] * syy_y;                                                
                           syy_y = syy_y / K_y[j] + psi_syy_y[j][i]; 
         }
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
							h = j;
		            			                              
                           psi_syy_y[h1][i] = b_y[h1] * psi_syy_y[h1][i] + a_y[h1] * syy_y;                                                
						   syy_y = syy_y / K_y[h1] + psi_syy_y[h1][i];

        }                       
                           b = 1.0 / rho[j][i]; 
			   if(rho[j][i]<1e-20){b = 0.0;}			    
			 
                           if(sw==0){vyp1[j][i] = b * (sxy_x+syy_y)/DH;}                           
                           
                           if(sw==1){
			   
			      if(GRAD_FORM==2){   /* FWI without data integration */
			        vyp1[j][i] = vy[j][i];
			      }
			      
			      if(GRAD_FORM==1){   /* FWI with data integration */
			        vyp1[j][i] += vy[j][i]*DT;
			      }
			   
			   }
                           
                           vy[j][i] += b * DT * (sxy_x+syy_y) / DH; 		         
      
     }
}

	     
		break;
		
	case 8:

for (j=ny1;j<=ny2;j++){
	for (i=nx1;i<=nx2;i++){

                           
                           sxy_x = hc[1]*(sxy[j][i]-sxy[j][i-1])
					    + hc[2]*(sxy[j][i+1]-sxy[j][i-2])
					    + hc[3]*(sxy[j][i+2]-sxy[j][i-3])
					    + hc[4]*(sxy[j][i+3]-sxy[j][i-4]);

		          
                           syy_y = hc[1]*(syz[j][i]-syz[j-1][i])
					    + hc[2]*(syz[j+1][i]-syz[j-2][i])
					    + hc[3]*(syz[j+2][i]-syz[j-3][i])
					    + hc[4]*(syz[j+3][i]-syz[j-4][i]);
                          

	/* left boundary */                                         
        if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){

                           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
         }

        /* right boundary */                                         
        if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
         }

         
	  /* top boundary */                                         
        if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
			   psi_syy_y[j][i] = b_y[j] * psi_syy_y[j][i] + a_y[j] * syy_y;                                                
                           syy_y = syy_y / K_y[j] + psi_syy_y[j][i]; 
                                 
		}
		

	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
		            h = j;
		            			   
       	                   
                           psi_syy_y[h1][i] = b_y[h1] * psi_syy_y[h1][i] + a_y[h1] * syy_y;                                                
						   syy_y = syy_y / K_y[h1] + psi_syy_y[h1][i];

        }                       
                           
			   b = 1.0 / rho[j][i]; 
			   if(rho[j][i]<1e-20){b = 0.0;}			    
			 
                           if(sw==0){vyp1[j][i] = b * (sxy_x+syy_y)/DH;}			                              
                           
                           if(sw==1){
			   
			      if(GRAD_FORM==2){   /* FWI without data integration */
			        vyp1[j][i] = vy[j][i];
			      }
			      
			      if(GRAD_FORM==1){   /* FWI with data integration */
			        vyp1[j][i] += vy[j][i]*DT;
			      }
			   
			   }
                           
                           vy[j][i] += b * DT * (sxy_x+syy_y) / DH; 		         
      
     }
}

	      
		break;
		
	case 10:
	
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			   
			   sxy_x =  hc[1]*(sxy[j][i]-sxy[j][i-1])
				 +  hc[2]*(sxy[j][i+1]-sxy[j][i-2])
				 +  hc[3]*(sxy[j][i+2]-sxy[j][i-3])
				 +  hc[4]*(sxy[j][i+3]-sxy[j][i-4])
				 +  hc[5]*(sxy[j][i+4]-sxy[j][i-5]);			   
                          
			   syy_y =  hc[1]*(syz[j][i]-syz[j-1][i])
				 +  hc[2]*(syz[j+1][i]-syz[j-2][i])
				 +  hc[3]*(syz[j+2][i]-syz[j-3][i])
				 +  hc[4]*(syz[j+3][i]-syz[j-4][i])
				 +  hc[5]*(syz[j+4][i]-syz[j-5][i]);
			  

        /* apply PML */
        if(FW>0){
        
		/* left boundary */                                         
		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                	           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                        	   sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
		}

		/* right boundary */                                         
		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                     
                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
		}

         
		/* top boundary */                                         
		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
			   psi_syy_y[j][i] = b_y[j] * psi_syy_y[j][i] + a_y[j] * syy_y;                                                
                           syy_y = syy_y / K_y[j] + psi_syy_y[j][i]; 
                        
                    
		}
		

		/* bottom boundary */                                         
		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
			    h = j;
		            		
                           
                           psi_syy_y[h1][i] = b_y[h1] * psi_syy_y[h1][i] + a_y[h1] * syy_y;                                                
			   syy_y = syy_y / K_y[h1] + psi_syy_y[h1][i]; 

		}
                              
        } /* end of PML */ 
                         
			   b = 1.0 / rho[j][i]; 
			   if(rho[j][i]<1e-20){b = 0.0;}			    
			 
                           if(sw==0){vyp1[j][i] = b * (sxy_x+syy_y)/DH;}


			   if(sw==1){
			   
			      if(GRAD_FORM==2){   /* FWI without data integration */
			        vyp1[j][i] = vy[j][i];
			      }
			      
			      if(GRAD_FORM==1){   /* FWI with data integration */
			        vyp1[j][i] += vy[j][i]*DT;
			      }
			   
			   }                                           
 
		           vy[j][i] += b * DT * (sxy_x+syy_y) / DH; 		         
      
     }
}

		break;

	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			   
			   sxy_x =  hc[1]*(sxy[j][i]-sxy[j][i-1])
				 +  hc[2]*(sxy[j][i+1]-sxy[j][i-2])
				 +  hc[3]*(sxy[j][i+2]-sxy[j][i-3])
				 +  hc[4]*(sxy[j][i+3]-sxy[j][i-4])
				 +  hc[5]*(sxy[j][i+4]-sxy[j][i-5])
				 +  hc[6]*(sxy[j][i+5]-sxy[j][i-6]);			   
                          
			   syy_y =  hc[1]*(syz[j][i]-syz[j-1][i])
				 +  hc[2]*(syz[j+1][i]-syz[j-2][i])
				 +  hc[3]*(syz[j+2][i]-syz[j-3][i])
				 +  hc[4]*(syz[j+3][i]-syz[j-4][i])
				 +  hc[5]*(syz[j+4][i]-syz[j-5][i])
				 +  hc[5]*(syz[j+5][i]-syz[j-6][i]);
			  

        /* apply PML */
        if(FW>0){
        
		/* left boundary */                                         
		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
			         
                	           psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * sxy_x;                                                
                        	   sxy_x = sxy_x / K_x[i] + psi_sxy_x[j][i];
              
		}

		/* right boundary */                                         
		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

                           h1 = (i-nx2+2*FW);
                            h=i; 
                                
                     
                           psi_sxy_x[j][h1] = b_x[h1] * psi_sxy_x[j][h1] + a_x[h1] * sxy_x;                                                
                           sxy_x = sxy_x / K_x[h1] + psi_sxy_x[j][h1];
		}

         
		/* top boundary */                                         
		if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
					   
			   psi_syy_y[j][i] = b_y[j] * psi_syy_y[j][i] + a_y[j] * syy_y;                                                
                           syy_y = syy_y / K_y[j] + psi_syy_y[j][i]; 
                        
                    
		}
		

		/* bottom boundary */                                         
		if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        
                           h1 = (j-ny2+2*FW);
			    h = j;
		            		
                           
                           psi_syy_y[h1][i] = b_y[h1] * psi_syy_y[h1][i] + a_y[h1] * syy_y;                                                
			   syy_y = syy_y / K_y[h1] + psi_syy_y[h1][i]; 

		}
                              
        } /* end of PML */ 
                         
			   b = 1.0 / rho[j][i]; 
			   if(rho[j][i]<1e-20){b = 0.0;}			    
			 
                           if(sw==0){vyp1[j][i] = b * (sxy_x+syy_y)/DH;}


			   if(sw==1){
			   
			      if(GRAD_FORM==2){   /* FWI without data integration */
			        vyp1[j][i] = vy[j][i];
			      }
			      
			      if(GRAD_FORM==1){   /* FWI with data integration */
			        vyp1[j][i] += vy[j][i]*DT;
			      }
			   
			   }                                           
 
		           vy[j][i] += b * DT * (sxy_x+syy_y) / DH; 		         
      
     }
}

		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxtmp = 0;
				vytmp = 0;
				for (m=1; m<=fdoh; m++) {
							
					vytmp +=   hc[m]*( syz[j+m][i]   - syz[j-m+1][i] )
						 + hc[m]*( sxy[j][i+m-1] - sxy[j][i-m]   );
				}

				vy[j][i] += vytmp*dtdh/rho[j][i];
			}
		}
		break;

	} /* end of switch(FDORDER) */


 		/* Forward Modelling (sw==0) */
	        if(sw==0){
	        for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];
		    azi_rad=srcpos_loc[7][l]*PI/180;
		    QUELLTYP=(int)srcpos_loc[8][l];
		
		    if(QUELLTYP==3){vy[j][i] +=  signals1[l][nt];}  /* single force in y */ 
		              
		}}
		
		/* Backpropagation (sw==1) */
		if(sw==1){
		for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];

		    if(QUELLTYPB==2){vy[j][i] += signals1[l][nt];}  /* single force in y */

		}}                         
			

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
