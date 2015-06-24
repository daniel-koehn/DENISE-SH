/**/
/*------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs                         
 *   last update 29/06/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){

	/* declaration of extern variables */
	extern int   NX, NY, FDORDER, MAXRELERROR, QUELLART, QUELLTYP, SNAP, SNAP_FORMAT, L;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, PLANE_WAVE_DEPTH, PHI;
	extern float XREC1, XREC2, YREC1, YREC2, FPML;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST, MUN, EPSILON, EPSILON_u, EPSILON_rho;
	extern int SEISMO, NDT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC;
	extern int BOUNDARY, REC_ARRAY, DRX, LOG, FW;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4], ANGLE;
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
	extern char JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE], INV_MODELFILE[STRING_SIZE];
	extern int RUN_MULTIPLE_SHOTS, TAPERLENGTH, INVTYPE;
	extern int NPROC, NPROCX, NPROCY, MYID, IDX, IDY, CHECKPTREAD, CHECKPTWRITE;
	extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT1, INVMAT, QUELLTYPB;
	extern int HESSIAN, GRAD_METHOD, ORDER_HESSIAN, NFREQ;
	extern float FC_HESSIAN, FC_HESS_START, FC_HESS_INC; 
	extern int MODEL_FILTER, FILT_SIZE;
	extern int FILT_SIZE_GRAD, GRAD_FILTER;
	
	extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR; 
	extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
	extern int SWS_TAPER_FILE;
	extern float SRTRADIUS, WD_DAMP;
	extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
	extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
	extern int MIN_ITER;
	extern int nfstart, nf;
	extern int nfstart_jac, nf_jac;
	extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
	extern float npower, k_max_PML;
	extern int INV_STF, N_STF, N_STF_START;
	extern char PARA[STRING_SIZE];
	extern int TIME_FILT, ORDER;
	extern float FC_START, FC_END, FC_INCR;
	extern int LNORM, DTINV;
	extern int STEPMAX;
	extern float EPS_SCALE, SCALEFAC;
	extern float PRO;
	extern int TRKILL;
	extern char TRKILL_FILE[STRING_SIZE];
	extern int TIMEWIN, NORMALIZE;
	extern float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
	extern char PICKS_FILE[STRING_SIZE];
	extern char MISFIT_LOG_FILE[STRING_SIZE];
        extern int TIMELAPSE;
        extern char DATA_DIR_T0[STRING_SIZE]; 
        
	/* definition of local variables */
	int idum[NPAR];
	float fdum[NPAR];
	
	
	if (MYID == 0){ 

	fdum[1]  = DH;                                                                                 
        fdum[2]  = TIME;                                                                               
        fdum[3]  = DT;                                                                                 
        fdum[4]  = TS;                                                                                 
        fdum[5]  = 0.0;                                                                                 
        fdum[6]  = 0.0;                                                                                 

        fdum[8]  = TAU;                                                                                                                                                                 
        fdum[10]  = TSNAP1;                                                                            
        fdum[11]  = TSNAP2;                                                                            
        fdum[12]  = TSNAPINC;                                                                          
        fdum[13]  = REFREC[1];                                                                         
        fdum[14]  = REFREC[2];                                                                         
        fdum[15]  = PHI;                                                                             

        fdum[16]  = XREC1;                                                                             
        fdum[17]  = YREC1;                                                                             

        fdum[19]  = XREC2;                                                                             
        fdum[20]  = YREC2;                                                                             

        fdum[22]  = DAMPING;                                                                           
        fdum[23]  = REC_ARRAY_DEPTH;                                                                   
        fdum[24]  = REC_ARRAY_DIST;                                                                    
        fdum[25]  = PLANE_WAVE_DEPTH;
	
	fdum[26]  = MUN;
	fdum[27]  = EPSILON;
	fdum[28]  = EPSILON_u;
	fdum[29]  = EPSILON_rho;
	fdum[30]  = ANGLE;
	fdum[31]  = FPML;
	
	fdum[32]  = SRTRADIUS;
	
	fdum[33]  = VPUPPERLIM;	
	fdum[34]  = VPLOWERLIM;	
	fdum[35]  = VSUPPERLIM;	
	fdum[36]  = VSLOWERLIM;	
	fdum[37]  = RHOUPPERLIM;	
	fdum[38]  = RHOLOWERLIM;
	
	fdum[39]  = npower;
	fdum[40]  = k_max_PML;
	
	fdum[41]  = FC_HESSIAN;
	
	fdum[42]  = FC_START;
	fdum[43]  = FC_END;
	fdum[44]  = FC_INCR;
	
	fdum[45]  = EPS_SCALE;
	fdum[46]  = SCALEFAC;
	fdum[47]  = PRO;
	fdum[48]  = WD_DAMP;
	fdum[49]  = FC_HESS_START;
	fdum[50]  = FC_HESS_INC;	
	
	                                                                                                                                                                                                                                                             
        idum[1]  = NPROCX;                                                                             
        idum[2]  = NPROCY;                                                                             
        idum[3]  = LOG;                                                                             

        idum[4]  = NPROCX*NPROCY;                                                               
        idum[5]  = NX;                                                                                 
        idum[6]  = NY;                                                                                 

        idum[8]  = QUELLART;                                                                           
        idum[9]  = QUELLTYP;                                                                           
        idum[10]  = READMOD;                                                                           
        idum[11]  = L;                                                                                 
        idum[12]  = FREE_SURF;                                                                         
        idum[13]  = SNAP;                                                                              
        idum[14]  = DRX;                                                                               

        idum[16]  = BOUNDARY;                                                                          
        idum[17]  = REC_ARRAY;                                                                         
        idum[18]  = SRCREC;                                                                                 
        idum[19]  = IDX;                                                                                 
        idum[20]  = IDY;                                                                                 
        idum[21]  = 0;                                                                                 
        idum[22]  = 0;                                                                                 
        idum[23]  = SNAP_FORMAT;                                                                       
        idum[24]  = SEISMO;                                                                            
        idum[25]  = READREC;                                                                           
        idum[26]  = NGEOPH;                                                                            
        idum[27]  = NDT;                                                                               
        idum[28]  = SEIS_FORMAT;                                                                       
        idum[29]  = CHECKPTREAD;                                                                       
        idum[30]  = CHECKPTWRITE;                                                                       
                                                                                                      
	idum[31]  = FDORDER;
	idum[32]  = MAXRELERROR;
	idum[33]  = RUN_MULTIPLE_SHOTS;	
        idum[34]  = TAPERLENGTH;
	idum[35]  = INVTYPE;
	idum[36]  = GRADT1;
	idum[37]  = GRADT2;
	idum[38]  = GRADT3;
	idum[39]  = GRADT4;
	idum[40]  = ITERMAX;
	idum[41]  = INVMAT1;
	idum[42]  = FW;
	idum[43]  = INVMAT;
	idum[44]  = QUELLTYPB;
	
	idum[45]  = TESTSHOT_START;
	idum[46]  = TESTSHOT_END;
	idum[47]  = TESTSHOT_INCR;
	
	idum[48]  = SWS_TAPER_GRAD_VERT;
	idum[49]  = SWS_TAPER_GRAD_HOR;
	idum[50]  = SWS_TAPER_GRAD_SOURCES;
	idum[51]  = SWS_TAPER_CIRCULAR_PER_SHOT;
	idum[52]  = SRTSHAPE;
	idum[53]  = FILTSIZE;
	
	idum[54]  = SPATFILTER;
	idum[55]  = SPAT_FILT_SIZE;
	idum[56]  = SPAT_FILT_1;
	idum[57]  = SPAT_FILT_ITER;
	
	idum[58]  = INV_RHO_ITER;
	idum[59]  = nfstart;
	idum[60]  = nf;
	
	idum[61]  = nfstart_jac;
	idum[62]  = nf_jac;
	idum[63]  = SWS_TAPER_FILE;
	idum[64]  = HESSIAN;
	idum[65]  = GRAD_METHOD;
	
	idum[66]  = MODEL_FILTER;
	idum[67]  = FILT_SIZE;

	idum[68]  = ORDER_HESSIAN;
	
	idum[69]  = INV_STF;
	idum[70]  = N_STF;
	idum[71]  = N_STF_START;
	
	idum[72]  = TIME_FILT;
	idum[73]  = ORDER;
	
	idum[74]  = LNORM;
	idum[75]  = DTINV;
	
	idum[76]  = STEPMAX;
	
	idum[77]  = TRKILL;

	idum[78]  = TIMEWIN;
	fdum[79]  = TWLENGTH_PLUS;
	fdum[80]  = TWLENGTH_MINUS;
	fdum[81]  = GAMMA;

	idum[82]  = NORMALIZE;
	
	idum[83]  = INV_VP_ITER;
	idum[84]  = INV_VS_ITER;
	
	idum[85]  = MIN_ITER;
	
	idum[86]  = GRAD_FILTER;
	idum[87]  = FILT_SIZE_GRAD;
	idum[88]  = TIMELAPSE;
	idum[89]  = NFREQ;

	
	} /** if (MYID == 0) **/
	
	if (MYID != 0) FL=vector(1,L);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);

	MPI_Bcast(&SOURCE_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SNAP_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&REC_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_VX,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_VY,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_CURL,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_DIV,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE_P,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&CHECKPTFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&JACOBIAN,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&DATA_DIR,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&INV_MODELFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&PARA,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MISFIT_LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD); 
	MPI_Bcast(&TRKILL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&PICKS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
        MPI_Bcast(&DATA_DIR_T0,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);

	DH=fdum[1];
	TIME=fdum[2];
	DT=fdum[3];
	TS=fdum[4];

	TAU=fdum[8];
	TSNAP1=fdum[10];
	TSNAP2=fdum[11];
	TSNAPINC=fdum[12];
	REFREC[1]=fdum[13];
	REFREC[2]=fdum[14];
	PHI=fdum[15];
	XREC1=fdum[16];
	YREC1=fdum[17];

	XREC2=fdum[19];
	YREC2=fdum[20];

	DAMPING=fdum[22];
	REC_ARRAY_DEPTH=fdum[23];
	REC_ARRAY_DIST=fdum[24];
	PLANE_WAVE_DEPTH=fdum[25];

        MUN = fdum[26];
	EPSILON = fdum[27];
	EPSILON_u = fdum[28];
	EPSILON_rho = fdum[29];
	ANGLE = fdum[30];
        FPML = fdum[31];
	
	SRTRADIUS = fdum[32];
	
	VPUPPERLIM = fdum[33];	
	VPLOWERLIM = fdum[34];	
	VSUPPERLIM = fdum[35];	
	VSLOWERLIM = fdum[36];	
	RHOUPPERLIM = fdum[37];	
	RHOLOWERLIM = fdum[38];
	
	npower = fdum[39];
	k_max_PML = fdum[40];
	
	FC_HESSIAN = fdum[41];
	
	FC_START = fdum[42];
	FC_END = fdum[43];
	FC_INCR = fdum[44];
	
	EPS_SCALE = fdum[45];
	SCALEFAC = fdum[46];
	
	PRO = fdum[47];
	WD_DAMP = fdum[48];
	FC_HESS_START = fdum[49];
	FC_HESS_INC = fdum[50];	
	

	NPROCX = idum[1];
	NPROCY = idum[2];
	LOG=idum[3];
	NPROC  = idum[4];
	NX = idum[5];
	NY = idum[6];

	QUELLART = idum[8];
	QUELLTYP = idum[9];
	READMOD = idum[10];
	L = idum[11];
	FREE_SURF = idum[12];
	SNAP = idum[13];
	DRX = idum[14];

	BOUNDARY = idum[16];
	REC_ARRAY = idum[17];
	SRCREC = idum[18];
	IDX = idum[19];
	IDY = idum[20];
	
	
	SNAP_FORMAT = idum[23];
	SEISMO = idum[24];
	READREC = idum[25];
	NGEOPH = idum[26];
	NDT = idum[27];
	SEIS_FORMAT = idum[28];
	CHECKPTREAD = idum[29];
	CHECKPTWRITE = idum[30];

	FDORDER = idum[31];
	MAXRELERROR = idum[32];
        RUN_MULTIPLE_SHOTS = idum[33];
        TAPERLENGTH = idum[34];
	INVTYPE = idum[35]; 
	GRADT1 = idum[36];
	GRADT2 = idum[37];
	GRADT3 = idum[38];
	GRADT4 = idum[39];
	ITERMAX = idum[40];
	INVMAT1 = idum[41];
	     FW = idum[42];
	INVMAT  = idum[43];  
	QUELLTYPB = idum[44];
	
	TESTSHOT_START = idum[45];
	TESTSHOT_END = idum[46];
	TESTSHOT_INCR = idum[47];
	
	SWS_TAPER_GRAD_VERT = idum[48];
	SWS_TAPER_GRAD_HOR = idum[49];
	SWS_TAPER_GRAD_SOURCES = idum[50];
	SWS_TAPER_CIRCULAR_PER_SHOT = idum[51];
	SRTSHAPE = idum[52];
	FILTSIZE = idum[53];
	
	SPATFILTER = idum[54];
	SPAT_FILT_SIZE = idum[55];
	SPAT_FILT_1 = idum[56];
	SPAT_FILT_ITER = idum[57];
	
	INV_RHO_ITER = idum[58];
	nfstart = idum[59];
	nf = idum[60];
	
	nfstart_jac = idum[61];
	nf_jac = idum[62];
	SWS_TAPER_FILE = idum[63];
	HESSIAN = idum[64];
	GRAD_METHOD = idum[65];

	
	MODEL_FILTER = idum[66];
	FILT_SIZE = idum[67];

	
	ORDER_HESSIAN = idum[68];
	
	INV_STF = idum[69];
	N_STF = idum[70];
	N_STF_START = idum[71];

	TIME_FILT = idum[72];
	ORDER = idum[73];
	
	LNORM = idum[74];
	DTINV = idum[75];
	
	STEPMAX = idum[76];
	
	TRKILL = idum[77];

	TIMEWIN = idum[78];
	TWLENGTH_PLUS = fdum[79];
	TWLENGTH_MINUS = fdum[80];
	GAMMA = fdum[81];

	NORMALIZE = idum[82];
	
	INV_VP_ITER = idum[83];
	INV_VS_ITER = idum[84];

	MIN_ITER = idum[85];
	
	GRAD_FILTER = idum[86];
	FILT_SIZE_GRAD = idum[87];
	TIMELAPSE = idum[88];
	NFREQ = idum[89];
	  
	MPI_Bcast(&FL[1],L,MPI_FLOAT,0,MPI_COMM_WORLD);

}
