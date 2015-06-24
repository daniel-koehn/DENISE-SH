/*------------------------------------------------------------------------
 * Pure forward modelling code
 *
 * D. Koehn
 * Kiel, 3.2.2012
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void forward_mod(FILE *fprec, float ** waveconv, float ** waveconv_rho, float ** waveconv_u, float ** prho, float ** prhonp1, float ** ppi, float ** ppinp1, int iter, float eps_scale, int nfstart,
	int nsrc, float ** puipjp, float ** prip, float ** prjp, float L2, int partest, float ** srcpos_loc, float ** srcpos, float ** srcpos1, float ** signals, int ns,
	int nd, float ** pvx, float ** pvy, float ** psxx, float ** psyy, float ** psxy, float ** ux, float ** uy, float ** pvxp1, float ** pvyp1, float ** psi_sxx_x, float ** psi_sxy_x,
	float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ** pvxm1, float ** pvym1, float ** uttx, 
	float ** utty, float ** absorb_coeff, float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y,  
	float * K_y_half, float * a_y_half, float * b_y_half, float ** uxy, float ** uyx, int ntr, int **recpos_loc, float **sectionvx, float **sectionvy, float **sectionp, float **sectioncurl, 
	float **sectiondiv, float **sectionread, int ntr_glob, float ** sectionvxdata, float ** sectionvxdiff, float ** sectionvxdiffold, float ** sectionvydata, float ** sectionvydiff, 
	float ** sectionvydiffold, int LNORM, int TIMEWIN, float * epst1, float * L2t, float L2sum, float energy_sum, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
	float ** buffertop_to_bot, float ** bufferbot_to_top, float **pu, float **punp1, int itest, int nsrc_loc,MPI_Request * req_send, MPI_Request * req_rec){

	extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY, FDORDER, INVMAT, MYID, QUELLART, QUELLTYP, SEISMO, QUELLTYPB;
	extern int RUN_MULTIPLE_SHOTS, TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, NT, NDT, FREE_SURF, TRKILL, REC1, REC2;
	
	extern float TSNAP1, DT;

	int h, i, j, n, nshots, ishot, nt, lsnap, lsamp, nsnap, infoout;
	float alphanom, alphadenom, tmp;
	
/* calculate change in the material parameters */
tmp=calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,INVMAT,eps_scale,1,nfstart);

if(MYID==0){printf("eps_scale = %e \n",eps_scale);}

/* For the calculation of the material parameters beteween gridpoints
   the have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
matcopy_elastic(prhonp1, ppinp1, punp1);
MPI_Barrier(MPI_COMM_WORLD);

av_mue(punp1,puipjp,prhonp1);
av_rho(prhonp1,prip,prjp);


/* initialization of L2 calculation */
L2=0.0;

alphanom = 0.0;
alphadenom = 0.0;

exchange_par();
 
if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
    
        for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){
 
        if(MYID==0){
        printf("\n=================================================================================================\n");
        printf("\n MYID=%d *****  Starting simulation (test-forward model) no. %d for shot %d of %d \n",MYID,itest,ishot,nshots);
        printf("\n=================================================================================================\n\n");}
		
        for (nt=1;nt<=6;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
        if (RUN_MULTIPLE_SHOTS){

	    /* find this single source positions on subdomains */
           if (nsrc_loc>0) free_matrix(srcpos_loc,1,6,1,1);
           srcpos_loc=splitsrc(srcpos1,&nsrc_loc,1);
		}
		           
		else{
		/* Distribute multiple source positions on subdomains */
		   srcpos_loc = splitsrc(srcpos,&nsrc_loc,nsrc);
		}

/* calculate wavelet for each source point */
signals=wavelet(srcpos_loc,nsrc_loc);
if (nsrc_loc){
	if(QUELLART==6){FFT_filt(signals,1.0,1,ns,1);}}
		    
/* initialize wavefield with zero */
zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);
     
/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");

		
   infoout = 0;

   /* update of particle velocities */
   update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);

   /* exchange of particle velocities between PEs */
   exchange_v(pvx, pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);

   update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, absorb_coeff, prhonp1, hc, infoout,
                        K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx); 
 
    /* explosive source */
   if (QUELLTYP==1) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc);

   if ((FREE_SURF) && (POS[2]==0))
         surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppinp1, punp1, prhonp1, hc, K_x, a_x, b_x, psi_vxxs);

   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);

	/* store amplitudes at receivers in section-arrays */
	if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppinp1, punp1, hc);
		lsamp+=NDT;
	} 		


   }/*--------------------  End  of loop over timesteps (test forward model) ----------*/


if (MYID==0){
printf("Calculate residuals between test forward model m - mu * dm and actual model m \n");
printf("----------------------------------------------------------------------------- \n");
}

if (ntr > 0){

/* read seismic data from SU file vx */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==3)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,1);

if(TRKILL==1){
     sectionread[3][1]=20000;
     sectionread[20][1]=20000;
}
          
/* assign input data to each PE */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,1,TIMEWIN,iter,ishot);
} /* end QUELLTYPB*/

/* read seismic data from SU file vy */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */
inseis(fprec,ishot,sectionread,ntr_glob,ns,2);

if(TRKILL==1){
     sectionread[1][1]=20000;
     sectionread[2][1]=20000;
     sectionread[3][1]=20000;
     sectionread[19][1]=20000;
     sectionread[20][1]=20000;
     sectionread[21][1]=20000;
     sectionread[22][1]=20000;
     sectionread[23][1]=20000;
}
          
/* assign input data to each PE */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,1,TIMEWIN,iter,ishot);			   	    

} /* end QUELLTYPB */

}

} /* ===========================================================================================================================*/   
/* ==================================== end of loop over shots (test forward) ==================================================*/
/* =============================================================================================================================*/

epst1[itest]=eps_scale;
epst1[1] = 0.0;

L2sum=0.0;
MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
if(LNORM==2){   
L2t[itest] = L2sum/energy_sum;}
else {L2t[itest] = L2sum;}
      
}

