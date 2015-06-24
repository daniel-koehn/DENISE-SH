/*------------------------------------------------------------------------
 * Step length estimation by parabolic line search
 *
 * D. Koehn
 * Kiel, 3.2.2012
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float line_search(FILE *fprec, float ** waveconv_rho, float ** waveconv_u, float ** prho, float ** prhonp1, int iter, int nfstart,
	int nsrc, float L2, int partest, float ** srcpos_loc, float ** srcpos, float ** srcpos1, float ** signals, int ns,
	int nd, float ** pvy, float ** psyz, float ** psxy, float ** uy, float ** pvyp1, float ** psi_sxy_x,
        float ** psi_vyx, float ** psi_syz_y, float ** psi_vyy, float ** pvym1, 
	float ** utty, float ** absorb_coeff, float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y,  
	float * K_y_half, float * a_y_half, float * b_y_half, float ** uyx, int ntr, int **recpos_loc, float **sectionvy, 
	float **sectionread, int ntr_glob, float ** sectionvydata, float ** sectionvydiff, 
	float ** sectionvydiffold, float * epst1, float * L2t, float L2sum, float energy_sum, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
	float ** buffertop_to_bot, float ** bufferbot_to_top, float **pu, float **punp1, float **puip, float **pujp, float ***pr, float ***pp, float ***pq, float **fipjp, float **f, float **g, 
	float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, float *etajm, float *etaip, float *peta, float **ptaus, float **ptausipjp, int itest, int nsrc_loc, int min_iter_help,
        int ** recpos, int nsrc_glob, int *step1, int *step3, float **gradg_u, float C_vs, float FC, MPI_Request * req_send, MPI_Request * req_rec){

/* global variables */
extern int MYID, TIME_FILT, MIN_ITER, INVMAT, TIMELAPSE, QUELLTYPB, INV_STF, QUELLART, LNORM, NDT,ORDER_HESSIAN, ORDER, SEISMO;
extern int NX, NY, L, NT, RUN_MULTIPLE_SHOTS, TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR;
extern float EPS_SCALE, SCALEFAC, TSNAP1, DT, FC_HESSIAN, FC_START;
extern int STEPMAX;

/* local variables */	
float opteps_vp, tmp, alphanom, alphadenom;
float eps_scale, scalefac ;
int h, i, j, n, nshots, ishot, nt, lsnap, lsamp, nsnap, infoout;

/* Variables for step length calculation */
int step2, itests, iteste, countstep, stepmax;

scalefac=SCALEFAC;
stepmax=STEPMAX;

*step1=0;
step2=0;

/* start with first guess for step length alpha */
eps_scale=EPS_SCALE; /* maximum model change = 1% of the maximum model value */
countstep=0; /* count number of forward calculations */

itests=2;
iteste=2;

/* set min_iter_help to initial global value of MIN_ITER */
if(iter==1){min_iter_help=MIN_ITER;}

while((step2!=1)||(*step1!=1)){

for (itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

/* calculate change in the material parameters */
tmp=calc_mat_change_test(waveconv_rho,waveconv_u,prho,prhonp1,pu,punp1,iter,1,INVMAT,eps_scale,1,nfstart);

/*char modfile[STRING_SIZE];

sprintf(modfile,"%s_vp_it_countstep%d.bin",INV_MODELFILE,countstep);
writemod(modfile,ppinp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_countstep%d.bin",INV_MODELFILE,countstep);

writemod(modfile,punp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_rho_it_countstep%d.bin",INV_MODELFILE,countstep);
writemod(modfile,prhonp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);*/


/*if(MYID==0){printf("EPSILON = %e \n",EPSILON);}*/

/* For the calculation of the material parameters beteween gridpoints
   the have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
matcopy_elastic(prhonp1, punp1);	/* no differentiation of elastic and viscoelastic modelling because the 
						viscoelastic parameters did not change during the forward modelling */
MPI_Barrier(MPI_COMM_WORLD);

av_mue(punp1,puip,pujp,prhonp1);
/*av_rho(prhonp1,prjp);*/


/* Preparing memory variables for update_s (viscoelastic) */
if (L) prepare_update_s(etajm,etaip,peta,fipjp,pujp,puip,prho,ptaus,ptausipjp,f,g,bip,bjm,cip,cjm,dip,d,e);
		
/* initialization of L2 calculation */
L2=0.0;

alphanom = 0.0;
alphadenom = 0.0;

exchange_par();
 
if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
    
for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){
 
        if(MYID==0){
           printf("\n=================================================================================================\n");
           printf("\n *****  Starting simulation (test-forward model) no. %d for shot %d of %d (rel. step length %.8f) \n",itest,ishot,nshots,eps_scale);
	   printf("\n=================================================================================================\n\n");
        }
		
        for (nt=1;nt<=8;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
        if (RUN_MULTIPLE_SHOTS){

	    /* find this single source positions on subdomains */
           if (nsrc_loc>0) free_matrix(srcpos_loc,1,8,1,1);
           srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		}
		           
		else{
		/* Distribute multiple source positions on subdomains */
		   srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
		}

/* calculate wavelet for each source point */
signals=wavelet(srcpos_loc,nsrc_loc,ishot);

if (nsrc_loc){if(QUELLART==6){timedomain_filt(signals,FC_HESSIAN,ORDER_HESSIAN,nsrc_loc,ns,1);}}

if((TIME_FILT)&&(INV_STF==0)){

  /*time domain low-pass filtering of the source signal */
  timedomain_filt(signals,FC,ORDER,nsrc_loc,ns,1);

  if(TIME_FILT==2){ /* band-pass */  
    timedomain_filt(signals,FC_START,ORDER,nsrc_loc,ns,2);
  }

}
		    
/* initialize wavefield with zero */
if (L){
	zero_fdveps_visc(-nd+1,NY+nd,-nd+1,NX+nd,pvy,psyz,psxy,uy,pvyp1,psi_sxy_x,psi_vyx,psi_syz_y,psi_vyy,pr,pp,pq);
}else{	
	zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvy,psyz,psxy,uy,pvyp1,psi_sxy_x,psi_vyx,psi_syz_y,psi_vyy);	
}     

/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");

		
   infoout = !(nt%10000);

   /*if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }*/

   /* update of particle velocities */
   update_v_PML(1, NX, 1, NY, nt, pvy, pvyp1, pvym1, utty, psyz, psxy, prhonp1, srcpos_loc,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_syz_y, psi_sxy_x);


   /*if (MYID==0){
      time4=MPI_Wtime();
      time_av_v_update+=(time4-time3);
     if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
      }*/

   /* exchange of particle velocities between PEs */
   exchange_v(pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);

   /*if (MYID==0){
      time5=MPI_Wtime();
	time_av_v_exchange+=(time5-time4);
      if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
      }*/

   if (L)    /* viscoelastic */
    	update_s_visc_PML(1, NX, 1, NY, pvy, uy, uyx, psyz, psxy, pujp, puip, prhonp1, hc, infoout,
				pr, pp, pq, fipjp, f, g, bip, bjm, cip, cjm, d, e, dip,
				K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vyx);
    else
   	update_s_elastic_PML(1, NX, 1, NY, pvy, uy, uyx, psyz, psxy, punp1, puip, pujp, absorb_coeff, prhonp1, hc, infoout,
         	               K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vyy, psi_vyx);  

   /*if (MYID==0){
      time6=MPI_Wtime();
	time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
      } */

   /*if ((FREE_SURF) && (POS[2]==0)){
   	if (L)   
		surface_PML(1, pvy, psyz, psxy, pp, pq, punp1, prhonp1, ptaus, etajm, peta, hc, K_x, a_x, b_x);
	else   
   		surface_elastic_PML(1, pvy, psyz, psxy, punp1, prhonp1, hc, K_x, a_x, b_x);
   }*/

   /* stress exchange between PEs */
    exchange_s(psyz,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);

   /*if (MYID==0){
      time7=MPI_Wtime();
 	time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }  */

	/* store amplitudes at receivers in section-arrays */
	if (SEISMO){
		seismo_ssg(nt, ntr, recpos_loc, sectionvy, pvy, psyz, punp1, hc);
		/*lsamp+=NDT;*/
	}

      
   /*if (MYID==0){
      time8=MPI_Wtime();
	time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }  */ 		


   }/*--------------------  End  of loop over timesteps (test forward model) ----------*/


if (MYID==0){
printf("Calculate residuals between test forward model m - mu * dm and actual model m \n");
printf("----------------------------------------------------------------------------- \n");
}

if (ntr > 0){

/* read seismic data from SU file vy */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,2,iter);

if (TIME_FILT){
   
   timedomain_filt(sectionread,FC,ORDER,ntr_glob,ns,1);

   if(TIME_FILT==2){ /* band-pass */
     timedomain_filt(sectionread,FC_START,ORDER,ntr_glob,ns,2);
   }
}

/* assign input data to each PE */
h=1;
for(i=1;i<=ntr;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[recpos_loc[3][i]][j];
   }
   h++;
}

/* Calculate v_mod(t1) - v_mod(t0) if TIMELAPSE == 1 */
/* ------------------------------------------------- */
if(TIMELAPSE==1){
  
      /* read synthetic seismic data at time step t0 vy */
      inseis(fprec,ishot,sectionread,ntr_glob,ns,10,iter);
                  
      /* calculate vy_mod(t1) - vy_mod(t0) */
      h=1;
      for(i=1;i<=ntr;i++){
	 for(j=1;j<=ns;j++){
	      sectionvy[h][j]=sectionvy[h][j]-sectionread[recpos_loc[3][i]][j];
	 }
	 h++;
      }
                                                                                                   
} /* end of TIMELAPSE */

L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,1,ntr_glob,recpos,recpos_loc,srcpos,nsrc_glob,ishot,iter);			   	    

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
     
} /* end of L2 test */

/* Did not found a step size which reduces the misfit function */
if((*step1==0)&&(L2t[1]<=L2t[2])){
 eps_scale = eps_scale/scalefac; 
 countstep++;
}

/* Found a step size with L2t[2] < L2t[3]*/
if((*step1==1)&&(L2t[2]<L2t[3])){
 epst1[3]=eps_scale;
 step2=1;
}

/* Could not found a step size with L2t[2] < L2t[3]*/
if((*step1==1)&&(L2t[2]>=L2t[3])){
 epst1[3]=eps_scale;
 /* increase step length to find  a larger misfit function than L2t[2]*/
 eps_scale = eps_scale + (eps_scale/scalefac);
 countstep++;                       
}         

/* found a step size which reduces the misfit function */
if((*step1==0)&&(L2t[1]>L2t[2])){
 epst1[2]=eps_scale; 
 *step1=1;
 iteste=3;
 itests=3;
 countstep=0;
 /* find a second step length with a larger misfit function than L2t[2]*/
 eps_scale = eps_scale + (eps_scale/scalefac);
}

*step3=0;

if((*step1==0)&&(countstep>stepmax)){
  if(MYID==0){
  printf(" Steplength estimation failed!");}
  *step3=1;
  break;
}

if((*step1==1)&&(countstep>stepmax)){
  if(MYID==0){
  printf("Could not found a proper 3rd step length which brackets the minimum\n");}
  *step1=1;
  step2=1;
}

if(MYID==0){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \n",iteste,itests,*step1,step2,eps_scale,countstep,stepmax,scalefac,MYID,L2t[1],L2t[2],L2t[3]);}

} /* end of while loop */

if(*step1==1){ /* only find an optimal step length if step1==1 */
/* calculate optimal step length epsilon for Vs and density */
if(MYID==0){
printf("======================================================== \n");
printf("calculate optimal step length epsilon for Vs and density \n");
printf("======================================================== \n");
}
eps_scale=calc_opt_step(L2t,waveconv_u,gradg_u,epst1,1,C_vs);
}

return eps_scale;
}

