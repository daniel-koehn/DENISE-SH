/*------------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD programs          
 *  last update  03/12/2000 
 *
 *  Copyright (c) 1998 T. Bohlen 
 *  See COPYING file for copying and redistribution conditions.
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 90
#define STRING_SIZE 74
#define STRING_SIZE2 256
#define REQUEST_COUNT 4


/* declaration of functions */

void window_cos(float **win, int npad, int nsrc, float it1, float it2, float it3, float it4);

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, MPI_Comm newcomm_nodentr);

int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch);

void absorb(float ** absorb_coeff);

void smooth_model(float ** unp1, float ** rho, int iter);

void taper_grad(float ** waveconv, float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int iter, int sws);

void taper_grad_shot(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int ishot);

void spat_filt(float ** waveconv, int iter, int sws);

float norm(float ** waveconv, int iter, int sws);

void av_mat(float **  u, float **  puip, float ** pujm);

void av_mue(float ** u, float ** uip, float ** ujp, float ** rho);

void av_rho(float **rho, float **rjp);

void av_tau(float **taus, float **tausipjp);

float median2d(float **mat, int ny, int nx);

float calc_mat_change_test(float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float  **  rhonp1, float **  u, float **  unp1, int iter, int epstest, int INVMAT, float eps_scale, int itest, int nfstart);

double calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot, int ntr_glob, int **recpos, int **recpos_loc, float **srcpos, int nsrc_glob, int ishot, int iter);

float calc_mres(int nx1, int nx2, int ny1, int ny2, float ** u, float ** u0);

double calc_misfit(float **sectiondata, float **section, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot);

float calc_opt_step(float *  L2t, float ** waveconv, float ** gradp, float * epst, int sws, float C_vp);

float calc_opt_step_test(float *  L2t, float ** waveconv, float ** gradp, float * epst, int sws, float C_vp);

double calc_energy(float **sectiondata, int ntr, int ns, float energy, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot);

void checkfd(FILE *fp, float ** prho, float ** pu,
float ** ptaus, float *peta);

void checkfd_hc(FILE *fp, float ** prho, float ** pu,
float ** ptaus, float *peta, float *hc);

void checkfd_ssg_elastic(FILE *fp, float ** prho, float ** pu, float *hc);
void checkfd_ssg_visc(FILE *fp, float ** prho, float ** pu, float ** ptaus, float * peta, float *hc);

void checkfd_rsg(FILE *fp, float ** prho, float ** pu,
float ** ptaus, float *peta);

void comm_ini(float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top, 
MPI_Request *req_send, MPI_Request *req_rec);

void conv_FD(float * temp_TS, float * temp_TS1, float * temp_conv, int ns);

void descent(float ** grad, float ** gradm);

void DFT(int ishot, int nt, float ** vy, float ** syy, float ** sxy, float *** green_vy, float *** greeni_vy, float *** green_syy, float *** greeni_syy, float *** green_sxy, float *** greeni_sxy);

void DFT_grad(int nt, float ** vy, float ** syy, float ** sxy, float ** green_vy, float ** greeni_vy, float ** green_syy, float ** greeni_syy,
         float ** green_sxy, float ** greeni_sxy, int back);

float dotp(float * vec1, float *vec2, int n1, int n2, int sw);

float exchange_L2(float L2, int sw, int bcast_l2);

void eprecond(float ** W, float ** vy);

void eprecond1(float ** We, float ** Ws, float ** Wr);

void exchange_rsg(float ** vy, float ** vz,
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top);

void exchange_rsg_4th(float ** vy, float ** vz,
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top);

void exchange_v(float ** vy,  
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top,
MPI_Request * req_send, MPI_Request * req_rec);

void exchange_s(float ** syy, 
float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top,
MPI_Request * req_send, MPI_Request * req_rec);

void exchange_par(void);

void exchange_mod_es(float ** matmod, int ncptot, int nparameter);

void extract_LBFGS_SH(int iter, float ** waveconv_u, float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float ** pu, float ** prho, float * r_LBFGS);

void extract_PCG_SH(float * PCG_old, float ** waveconv_u, float ** waveconv_rho);

void gauss_filt(float ** waveconv);

void grad_FD(float ** grad_vy, float ** gradi_vy, float ** grad_syz, float ** gradi_syz, 
float ** grad_sxy, float ** gradi_sxy, float ** grad_vyb, float ** gradi_vyb, float ** grad_syzb, 
float ** gradi_syzb, float ** grad_sxyb, float ** gradi_sxyb, float ** pu, float ** prho,
float ** waveconv_u_shot, float ** waveconv_rho_shot);

void hess_save(int nt, float ** uy, float ** uyx, float ** vyp1, float *** hess_eyy, float *** hess_eyx, float *** hess_utty);

void hessian(int NTDTINV,  float *** hess_eyy, float *** hess_eyx, float *** hess_utty, float ** hessian_u, float ** hessian_rho);

void hessian_out(float **hessian_u, float **hessian_rho, float ** u, float ** rho);

float *holbergcoeff(void);

int householder(int m, int n, float **mat, float *b);

void info(FILE *fp);

void initproc(void);

void interpol(int ni1, int ni2, float **  intvar, int cfgt_check);

void LBFGS(int iter, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, float * alpha_LBFGS, float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec);

float line_search(FILE *fprec, float ** waveconv_rho, float ** waveconv_u, float ** prho, float ** prhonp1, int iter, int nfstart,
	int nsrc, float L2, int partest, float ** srcpos_loc, float ** srcpos, float ** srcpos1, float ** signals, int ns,
	int nd, float ** pvy, float ** psyz, float ** psxy, float ** uy, float ** pvyp1, float ** psi_sxy_x,
        float ** psi_vyx, float ** psi_syz_y, float ** psi_vyy, float ** pvym1, 
	float ** utty, float ** absorb_coeff, float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y,  
	float * K_y_half, float * a_y_half, float * b_y_half, float ** uyx, int ntr, int **recpos_loc, float **sectionvy, 
	float **sectionread, int ntr_glob, float ** sectionvydata, float ** sectionvydiff, 
	float ** sectionvydiffold, float * epst1, float * L2t, float L2sum, float energy_sum, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
	float ** buffertop_to_bot, float ** bufferbot_to_top, float **pu, float **punp1, float ** puip, float ** pujp, float ***pr, float ***pp, float ***pq, float **fipjp, float **f, float **g, 
	float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, float *etajm, float *etaip, float *peta, float **ptaus, float **ptausipjp, int itest, int nsrc_loc, 
        int min_iter_help, int ** recpos, int nsrc_glob, int *step1, int *step3, float **gradg_u, float C_vs, float FC, MPI_Request * req_send, MPI_Request * req_rec);

double LU_decomp(double  **A, double *x, double *b,int n);

float minimum_m(float **mat, int nx, int ny);
float maximum_m(float **mat, int nx, int ny);

void mer(float **sectionp, int ntr, int nst, float *picked_times, int ishot);

void model(float  **  rho, float **  u, 
float **  taus, float *  eta);

void model_elastic(float  **  rho, float **  u);

void model_ani(float  **  rho, float **  c11, float **  c15, float **  c13, 
float **  c35, float **  c33, float **  c55, 
float **  taus, float *  eta);

void model_freq_out(float  **  rho, float **  unp1, int iter, float freq);

void matcopy(float ** prho, float ** pu, float ** ptaus);

void matcopy_elastic(float ** prho, float ** pu);

void matcopy_ani(float ** rho, float **  c11, float **  c15, float **  c13, 
float **  c35, float **  c33, float **  c55, float ** taus);

void matinv(float ** matrix, int n);

void matp(float ** mat1, float **mat2, float **matr, int n1, int n2, int m1, int m2, int sws);

void max_grad(float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float **  u);

void median_src(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int iter, int sws);
			  
void merge(int nsnap, int type);

void merge2(int nsnap, int type);

void mergemod(char modfile[STRING_SIZE], int format);

void note(FILE *fp);

void normalize_data_rms_MPI(float **data, float **data_obs, int ntr, int ns);

void  outseis(FILE *fp, FILE *fpdata, int comp, float **section,
int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
int nsrc, int ns, int seis_form, int ishot, int sws);

void  outseis_glob(FILE *fp, FILE *fpdata, int comp, float **section,
int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
int nsrc, int ns, int seis_form, int ishot, int sws);

void  outseis_vector(FILE *fp, FILE *fpdata, int comp, float *section,
int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
int nsrc, int ns, int seis_form, int ishot, int sws);

void  inseis(FILE *fp, int comp, float **section, int ntr, int ns, int sws, int iter);

void  taper(float **sectionpdiff, int ntr, int ns);

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

void PCG(float * PCG_new, float * PCG_old, float * PCG_dir, int PCG_class);

void PML_pro(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
float * d_y, float * K_y, float * alpha_prime_y, float * a_y, float * b_y, 
float * d_y_half, float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half);

float *rd_sour(int *nts,FILE* fp_source);

float readdsk(FILE *fp_in, int format);

void readbufs(float ** syy, 
float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top);

void readbufv(float ** vy,  
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2, float ** vy, float ** syy, float ** sxy);

void read_par(FILE *fp_in);

void read_par_json(FILE *fp, char *fileinp);

void read_par_inv(FILE *fp, int nstage, int stagemax);

void readmod(float  **  rho, float **  u, 
float **  taus, float *  eta);

void readmod_elastic(float  **  rho, float **  u);

void readmod_elastic_es(float  **  rho, float **  u, float ** matmod, int is);

int **receiver(FILE *fp, int *ntr);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2, float ** vy, float ** syy, float ** sxy);

void saveseis(FILE *fp, float **sectionvy, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos_loc, int nsrc,int ns, int iter);

void saveseis_glob(FILE *fp, float **sectionvy, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos_loc, int nsrc,int ns, int iter);

void  smooth2(float ** grad);

void snap(FILE *fp,int nt, int nsnap, float **vy,
	float **syy, float **u, float *hc);


void snap_rsg(FILE *fp,int nt, int nsnap, float **vy, float **syy, float **u);

void snapmerge(int nsnap);

void snapmerge2(int nsnap);

float **sources(int *nsrc);

void solvelin(float  **AA, float *bb, float *x, int e, int method);

void seismo(int lsamp, int ntr, int **recpos, 
float **sectionvy, float **pvy, float **psyy, float **pu); 

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvy, float **pvy, float **psyy, float **pu,
float *hc); 

void seismo_rsg(int lsamp, int ntr, int **recpos, 
float **sectionvy, float **pvy, float **psyy, float **pu); 

float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc);

float **splitsrc_back(int **recpos,int *nsrc_loc, int nsrc);

void stalta(float **sectionp, int ntr, int nst, float *picked_times, int ishot);

void stf(float **sectionvy_obs, float **sectionvy, int ntr_glob, int ishot, int ns, int iter, int nshots, float **signals, int **recpos, float **srcpos);

void store_LBFGS_SH(float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float ** waveconv_u, float ** gradp_u, float ** waveconv_rho, float ** gradp_rho, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float ** pu, float ** prho, int nxnyi, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec);

void store_PCG_SH(float * PCG_old, float ** waveconv_u, float ** waveconv_rho);

void time_window_stf(float **sectiondata, int iter, int ntr_glob, int ns, int ishot);

void surface(int ndepth, float ** pvy, 
float ** psyy,
float ** psxy, float *** pp, float *** pq, float  **  pu,
float ** ptaus, float * etajm, float * peta);

void surface_elastic(int ndepth, float ** vy, float ** syy,
float ** sxy, float  **  u, float ** rho, float * hc);

void surface_elastic_PML(int ndepth, float ** vy, float ** syy,
float ** sxy, float  **  u, float ** rho, float * hc, float * K_x, float * a_x, float * b_x);

void surface_PML(int ndepth, float ** vy, float ** syy,
float ** sxy, float ***p, float ***q, float  **  pu, float **prho, float **ptaus, 
float *etajm, float *peta, float * hc, float * K_x, float * a_x, float * b_x);

void  timedomain_filt(float ** data, float fc, int order, int ntr, int ns, int method);
void  timedomain_filt_vector(float * data, float fc, int order, int ntr, int ns, int method);

void time_window(float **sectiondata, float * picked_times, int iter, int ntr_glob, int **recpos_loc, int ntr, int ns, int ishot);

void prepare_update_s(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
float **puipjp, float **prho, float **ptaus,
float **ptausipjp, float **f, float **g, float *bip, float *bjm,
float *cip, float *cjm, float ***dip, float ***d, float ***e);

void update_s(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   syy,
float **   sxy, float *** r, float *** p, float *** q, float ** pu, float ** taus, 
float *   etaip, float *   etajm, float * peta);

void update_s_visc_hc(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   syy,
float **   sxy, float *** r, float *** p, float *** q, float ** pu, float **uipjp, float ** taus, 
float **tausipjp, float *   etaip, float *   etajm, float * peta, float *hc);


void update_s_rsg(int nx1, int nx2, int ny1, int ny2, float ** pvy, float ** psyy,
float ** psxy, float *** pr, float *** pp, float ***pq, float  **  pu,
float ** ptaus, float * etaip,
float * etajm, float * peta, float ** absorb_coeff);

void update_s_rsg_4th(int nx1, int nx2, int ny1, int ny2, float ** pvy, float ** psyy,
float ** psxy, float *** pr, float *** pp, float ***pq, float  **  pu,
float ** ptaus, float * etaip,
float * etajm, float * peta);

void update_s_ani(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   syy,
float **   sxy, float *** r, float *** p, float *** q,
float  ** c11, float  **  c15, float  ** c13, float  **  c35, 
float  ** c33, float  **  c55, float **   ptaus, 
float *   etaip, float *   etajm, float * peta);

void update_s_elastic(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   syy,
float **   sxy, float ** u, float ** uipjm, float ** absorb_coeff);

void update_s_elastic_rsg(int nx1, int nx2, int ny1, int ny2, float **  vy, float **  syy,
float **  sxy, float  **   u, float ** absorb_coeff);

void update_s_elastic_hc(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   uy, float **   uyx, float **   syy,
	float **   sxy, float ** u, float ** uipjp, float ** absorb_coeff, float ** rho,
	float *hc, int infoout);

void update_s_elastic_PML(int nx1, int nx2, int ny1, int ny2, 
float **   vy, float **   uy, float **   uyx, float **   syy,
float **   sxy, float ** u, float ** uip, float ** ujp, float ** absorb_coeff, float **rho, float *hc, int infoout,
float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
float ** psi_vyy, float ** psi_vyx);

void update_s_elastic_hh(int nx1, int nx2, int ny1, int ny2, float **   vy, float **   syy,
float **   sxy, float ** u );

void update_s_visc_PML(int nx1, int nx2, int ny1, int ny2,
float **   vy, float **   uy, float **   uyx, float **   syy,
float **   sxy, float ** u, float ** uipjp, float **rho, float *hc, int infoout,
float ***r, float ***p, float ***q, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip, 
float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
float ** psi_vyy, float ** psi_vyx);


void update_v(int nx1, int nx2, int ny1, int ny2, int nt, 
float ** pvy, float ** psyy,
float ** psxy, float  ** prho,  
float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

	
void update_v_hc(int nx1, int nx2, int ny1, int ny2, int nt, 
float ** vy, float **  vyp1, float **  vym1, float **  utty, float ** syy,
float ** sxy, float **rjp, float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff,
float *hc, int infoout, int sw);

void update_v_PML(int nx1, int nx2, int ny1, int ny2, int nt,
float ** vy, float **  vyp1, float **  vym1, float **  utty, float ** syy,
float ** sxy, float **rjp, float **  srcpos_loc, float ** signals1, int nsrc, float ** absorb_coeff,
float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
float ** psi_syy_y, float ** psi_syx_x);

void update_v_hh(int nx1, int nx2, int ny1, int ny2, int nt,
float ** pvy, float ** psyy,
float ** psxy, float  ** prho,  
float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

void update_v_rsg(int nx1, int nx2, int ny1, int ny2, int nt,
float ** pvy, float ** psyy,
float ** psxy, float  ** prho,  
float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

void update_v_rsg_4th(int nx1, int nx2, int ny1, int ny2, int nt,
float ** pvy, float ** psyy,
float ** psxy, float  ** prho,  
float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

float ** wavelet(float ** srcpos_loc, int nsrc, int ishot);
float ** wavelet_stf(int nsrc, int ishot, float ** signals_stf);

void  wavenumber(float ** grad);

void writebufs(float ** syy, 
float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top);

void writebufv(float ** vy,
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top);

void write_par(FILE *fp);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float ** array, int format);

void zero_fdveps(int ny1, int ny2, int nx1, int nx2, float ** vy,
                 float ** syy, float ** sxy, float ** vym1, float ** vyp1,
                 float ** psi_sxy_x, float ** psi_vyx, float ** psi_syy_y, float ** psi_vyy);

void zero_fdveps_visc(int ny1, int ny2, int nx1, int nx2, float ** vy,
                 float ** syy, float ** sxy, float ** vym1, float ** vyp1,
                 float ** psi_sxy_x, float ** psi_vyx, float ** psi_syy_y, float ** psi_vyy,
                 float ***pr, float ***pp, float ***pq);

void zero_FWI(float ** grad_vy,float ** gradi_vy,float ** grad_syz,float ** gradi_syz,float ** grad_sxy,float ** gradi_sxy);

void zero_hessian(int ny1, int ny2, int nx1, int nx2, int nshots, float *** green_vy, float *** greeni_vy, float *** green_syy, 
                 float *** greeni_syy, float *** green_sxy, float *** greeni_sxy);

void zero_LBFGS(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float * r_LBFGS, 
                 float * alpha_LBFGS, float * beta_LBFGS, float * rho_LBFGS);

void zero_PCG(float * PCG_old, float * PCG_new, float * PCG_dir, int PCG_vec);
		 
void FLnode(float  **  rho, float **  u, float **  taus, float *  eta);

void smooth_grad(float ** waveconv);

/* declaration of functions for parser*/

/* declaration of functions for json parser in json_parser.c*/
int read_objects_from_intputfile(FILE *fp, char input_file[STRING_SIZE],char ** varname_list,char ** value_list);

void print_objectlist_screen(FILE *fp, int number_readobject,char ** varname_list,char ** value_list);

int count_occure_charinstring(char stringline[STRING_SIZE], char teststring[]);

void copy_str2str_uptochar(char string_in[STRING_SIZE], char string_out[STRING_SIZE], char teststring[]);

int get_int_from_objectlist(char string_in[STRING_SIZE], int number_readobject, int * int_buffer,
		char ** varname_list,char ** value_list);

int get_float_from_objectlist(char string_in[STRING_SIZE], int number_readobject, float * double_buffer,
		char ** varname_list,char ** value_list);

int get_string_from_objectlist(char string_in[STRING_SIZE], int number_readobject, char string_buffer[STRING_SIZE],
		char ** varname_list,char ** value_list);

int is_string_blankspace(char string_in[STRING_SIZE]);

void remove_blankspaces_around_string(char string_in[STRING_SIZE] );

void add_object_tolist(char string_name[STRING_SIZE],char string_value[STRING_SIZE], int * number_read_object,
		char ** varname_list,char ** value_list );

/* utility functions */
void err(char err_text[]);
void warning(char warn_text[]);
double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
double *dvector(int nl, int nh);
float **fmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);

float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh); 
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, 
int ndh);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);
void normalize_data_rel(float ** sectionvy, float ** sectionvy_obs, int ntr_glob, int ns);
void normalize_data_rms(float **data, float **data_obs, int ntr, int ns);
void normalize_STF(float *data, int ns);
void tripd(float *d, float *e, float *b, int n);

