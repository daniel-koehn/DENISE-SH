/* $Id: update_s_elastic_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic_PML(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   uy, float **   uyx, float **   syz,
	float **   sxy, float ** u, float ** uip, float ** ujp, float ** absorb_coeff, float **rho, float *hc, int infoout,
      float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
      float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
      float ** psi_vyy, float ** psi_vyx){


	int i,j, m, fdoh, h, h1;
	float  vyy, vyx;
	float  dhi;	
	extern float DT, DH;
	extern int MYID, FDORDER, INVMAT1, FW;
        extern int FREE_SURF, BOUNDARY, GRAD_FORM;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;
	
	

	dhi = DT/DH;
	fdoh = FDORDER/2;

	
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
                                                
                        psi_vyy[j][i] = b_y_half[j] * psi_vyy[j][i] + a_y_half[j] * vyy;                                            
                     
                        vyy = vyy / K_y_half[j] + psi_vyy[j][i];

	        }
	
	        /* bottom boundary */                                         
	        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                        vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		}
       }
                         	if(GRAD_FORM==2){
				  uy[j][i] = vyy;
				  uyx[j][i] = vyx;
				}
				
				sxy[j][i] += uip[j][i]*vyx;
				syz[j][i] += ujp[j][i]*vyy;
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
                                                
                        psi_vyy[j][i] = b_y_half[j] * psi_vyy[j][i] + a_y_half[j] * vyy;                                            
                        vyy = vyy / K_y_half[j] + psi_vyy[j][i];
                        
        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                        vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
                                
        }

				if(GRAD_FORM==2){
				   uy[j][i] = vyy;
				   uyx[j][i] = vyx;
				}
				
				sxy[j][i] += uip[j][i]*vyx;
				syz[j][i] += ujp[j][i]*vyy;
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
                                                
                        psi_vyy[j][i] = b_y_half[j] * psi_vyy[j][i] + a_y_half[j] * vyy;                                            
                        vyy = vyy / K_y_half[j] + psi_vyy[j][i];
                        
        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                        vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
        }
                                if(GRAD_FORM==2){
				   uy[j][i] = vyy;
				   uyx[j][i] = vyx;
				}
				
				sxy[j][i] += uip[j][i]*vyx;
				syz[j][i] += ujp[j][i]*vyy;
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
                                                
                        psi_vyy[j][i] = b_y_half[j] * psi_vyy[j][i] + a_y_half[j] * vyy;                                            
                        vyy = vyy / K_y_half[j] + psi_vyy[j][i];

        }
	
	  /* bottom boundary */                                         
        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                        vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
                                
        }
				if(GRAD_FORM==2){
				   uy[j][i] = vyy;
				   uyx[j][i] = vyx;
				}
				
				sxy[j][i] += uip[j][i]*vyx;
				syz[j][i] += ujp[j][i]*vyy;

   }}
		break;

	case 10:

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){					

                       vyx = (hc[1]*(vy[j][i+1]-vy[j][i]) 
			     + hc[2]*(vy[j][i+2]-vy[j][i-1])
			     + hc[3]*(vy[j][i+3]-vy[j][i-2])
			     + hc[4]*(vy[j][i+4]-vy[j][i-3])
			     + hc[5]*(vy[j][i+5]-vy[j][i-4]))*dhi;

                        vyy = (hc[1]*(vy[j+1][i]-vy[j][i])
		             + hc[2]*(vy[j+2][i]-vy[j-1][i])
		             + hc[3]*(vy[j+3][i]-vy[j-2][i])
			     + hc[4]*(vy[j+4][i]-vy[j-3][i])
			     + hc[5]*(vy[j+5][i]-vy[j-4][i]))*dhi;
        
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
                                                
                        psi_vyy[j][i] = b_y_half[j] * psi_vyy[j][i] + a_y_half[j] * vyy;                                            
                     
                        vyy = vyy / K_y_half[j] + psi_vyy[j][i];

	        }
	
	        /* bottom boundary */                                         
	        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                        vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		}
       }
                         	if(GRAD_FORM==2){
				  uy[j][i] = vyy;
				  uyx[j][i] = vyx;
				}
				
				sxy[j][i] += uip[j][i]*vyx;
				syz[j][i] += ujp[j][i]*vyy;
			}
		}


		break;
		
	case 12:


		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){					

                       vyx = (hc[1]*(vy[j][i+1]-vy[j][i]) 
			     + hc[2]*(vy[j][i+2]-vy[j][i-1])
			     + hc[3]*(vy[j][i+3]-vy[j][i-2])
			     + hc[4]*(vy[j][i+4]-vy[j][i-3])
			     + hc[5]*(vy[j][i+5]-vy[j][i-4])
			     + hc[6]*(vy[j][i+6]-vy[j][i-5]))*dhi;

                        vyy = (hc[1]*(vy[j+1][i]-vy[j][i])
		             + hc[2]*(vy[j+2][i]-vy[j-1][i])
		             + hc[3]*(vy[j+3][i]-vy[j-2][i])
			     + hc[4]*(vy[j+4][i]-vy[j-3][i])
			     + hc[5]*(vy[j+5][i]-vy[j-4][i])
			     + hc[6]*(vy[j+6][i]-vy[j-5][i]))*dhi;
        
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
                                                
                        psi_vyy[j][i] = b_y_half[j] * psi_vyy[j][i] + a_y_half[j] * vyy;                                            
                     
                        vyy = vyy / K_y_half[j] + psi_vyy[j][i];

	        }
	
	        /* bottom boundary */                                         
	        if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

                        h1 = (j-ny2+2*FW);                                        
                        h = j;
                                                
                        psi_vyy[h1][i] = b_y_half[h1] * psi_vyy[h1][i] + a_y_half[h1] * vyy;                                            
                        vyy = vyy / K_y_half[h1] + psi_vyy[h1][i];
        
		}
       }
                         	if(GRAD_FORM==2){
				  uy[j][i] = vyy;
				  uyx[j][i] = vyx;
				}
				
				sxy[j][i] += uip[j][i]*vyx;
				syz[j][i] += ujp[j][i]*vyy;
			}
		}


		break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				
				vyy = 0.0;
				vyx = 0.0;
				
				for (m=1; m<=fdoh; m++) {
					vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
					vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
					
				}	

				sxy[j][i]+=uip[j][i]*vyx*dhi;
				syz[j][i]+=ujp[j][i]*vyy*dhi;
			}
		}
		break;
		
	} /* end of switch(FDORDER) */


	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
