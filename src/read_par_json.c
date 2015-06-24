/*------------------------------------------------------------------------
 *   program DENISE, reading input-parameters from input-file or stdin
  *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include "fd.h"

char ** varname_list,** value_list;

void read_par_json(FILE *fp, char *fileinp){

	/* declaration of extern variables */
	extern int   NX, NY, FDORDER, MAXRELERROR, QUELLART, QUELLTYP, SNAP, SNAP_FORMAT, L;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, PLANE_WAVE_DEPTH, PHI;
	extern float XREC1, XREC2, YREC1, YREC2, FPML;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST, ANGLE;
	extern int SEISMO, NDT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC, RUN_MULTIPLE_SHOTS;
	extern int BOUNDARY, REC_ARRAY, DRX, LOG, TAPER, TAPERLENGTH, INVTYPE, FW;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], LOG_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE], CHECKPTFILE[STRING_SIZE];
	extern char SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];
	extern char JACOBIAN[STRING_SIZE],DATA_DIR[STRING_SIZE];
	extern int  NPROCX, NPROCY, MYID, IDX, IDY, CHECKPTREAD, CHECKPTWRITE; 
	extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT1, INVMAT, QUELLTYPB;
	extern int HESSIAN, GRAD_METHOD, ORDER_HESSIAN;
	extern float FC_HESSIAN;
	extern int FILT_SIZE, MODEL_FILTER;
	extern int FILT_SIZE_GRAD, GRAD_FILTER;

	extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR; 
	extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
	extern int SWS_TAPER_FILE;
	extern float SRTRADIUS;
	extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
	extern int INV_RHO_ITER, INV_VS_ITER, INV_VP_ITER;
	extern char INV_MODELFILE[STRING_SIZE];
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
	extern int MIN_ITER;

	extern int TRKILL;
	extern char TRKILL_FILE[STRING_SIZE];

	extern int TIMEWIN, NORMALIZE;
	extern float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
	extern char PICKS_FILE[STRING_SIZE];
	
	extern char MISFIT_LOG_FILE[STRING_SIZE];
	
	
	
	/* definition of local variables */

	int number_readobjects=0,fserr=0;
	char errormessage[STRING_SIZE2];

	char ** varname_list, ** value_list;

	if (MYID == 0){

		/* allocate first object in list */
		varname_list = malloc(STRING_SIZE2*sizeof(char*));
		value_list = malloc(STRING_SIZE2*sizeof(char*));

		/* read in objects from file */
		number_readobjects=read_objects_from_intputfile(fp, fileinp, varname_list, value_list);
		fprintf(fp,"\nFrom input file %s, %i objects have been read in. \n",fileinp, number_readobjects);

		/* print objects to screen */
		fprintf(fp, "\n===========================================================");
		fprintf(fp, "\n=   List of Parameters read by the built in Json Parser   =");
		print_objectlist_screen(fp, number_readobjects, varname_list, value_list);

		/* extract variables form object list */

		/*=================================
		section general grid and discretization parameters
		  =================================*/
		if (get_int_from_objectlist("NPROCX",number_readobjects,&NPROCX,varname_list, value_list))
			err("Variable NPROCX could not be retrieved from the json input file!");
		if (get_int_from_objectlist("NPROCY",number_readobjects,&NPROCY,varname_list, value_list))
			err("Variable NPROCY could not be retrieved from the json input file!");
		if (get_int_from_objectlist("FDORDER",number_readobjects,&FDORDER,varname_list, value_list))
			err("Variable FDORDER could not be retrieved from the json input file!");
		if (get_int_from_objectlist("MAXRELERROR",number_readobjects,&MAXRELERROR,varname_list, value_list))
			err("Variable MAXRELERROR could not be retrieved from the json input file!");
		if (get_int_from_objectlist("NX",number_readobjects,&NX,varname_list, value_list))
			err("Variable NX could not be retrieved from the json input file!");
		if (get_int_from_objectlist("NY",number_readobjects,&NY,varname_list, value_list))
			err("Variable NY could not be retrieved from the json input file!");
		if (get_float_from_objectlist("DH",number_readobjects,&DH,varname_list, value_list))
			err("Variable DH could not be retrieved from the json input file!");
		if (get_float_from_objectlist("TIME",number_readobjects,&TIME,varname_list, value_list))
			err("Variable TIME could not be retrieved from the json input file!");
		if (get_float_from_objectlist("DT",number_readobjects,&DT,varname_list, value_list))
			err("Variable DT could not be retrieved from the json input file!");

		/*=================================
		 	 section source parameters
		  =================================*/
		fprintf(fp,"The following default values are set:\n");
		fprintf(fp,"=====================================\n\n");
		
		if (get_int_from_objectlist("QUELLTYP",number_readobjects,&QUELLTYP,varname_list, value_list))
			err("Variable SOURCE_TYPE could not be retrieved from the json input file!");
		else {
			if (QUELLTYP==4) {
				if (get_float_from_objectlist("ANGLE",number_readobjects,&ANGLE,varname_list, value_list))
					fprintf(fp,"\n WARNING: No angle for the rotated force is specified in the input file!\n");
					ANGLE=0.0;
					fprintf(fp,"          ANGLE is set to %f.\n",ANGLE);
			}
		}
		if (get_int_from_objectlist("QUELLART",number_readobjects,&QUELLART,varname_list, value_list))
			err("Variable QUELLART could not be retrieved from the json input file!");
		else {
			if (QUELLART==3) {
				if (get_string_from_objectlist("SIGNAL_FILE",number_readobjects,SIGNAL_FILE,varname_list, value_list))
					err("Variable SIGNAL_FILE could not be retrieved from the json input file!");
			}
		}
		if (get_int_from_objectlist("SRCREC",number_readobjects,&SRCREC,varname_list, value_list)){
			SRCREC=1;
			fprintf(fp,"Variable SRCREC is set to default value %d\n",SRCREC);}
		
		if (SRCREC==1) {
			if (get_string_from_objectlist("SOURCE_FILE",number_readobjects,SOURCE_FILE,varname_list, value_list))
				err("Variable SOURCE_FILE could not be retrieved from the json input file!");
			if (get_int_from_objectlist("RUN_MULTIPLE_SHOTS",number_readobjects,&RUN_MULTIPLE_SHOTS,varname_list, value_list))
				err("Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");
		}
		if (SRCREC==2) {
			if (get_float_from_objectlist("PLANE_WAVE_DEPTH",number_readobjects,&PLANE_WAVE_DEPTH,varname_list, value_list))
				err("Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");
			else {
				if (PLANE_WAVE_DEPTH>0.0) {
					if (get_float_from_objectlist("PHI",number_readobjects,&PHI,varname_list, value_list))
						err("Variable PHI could not be retrieved from the json input file!");
					if (get_float_from_objectlist("TS",number_readobjects,&TS,varname_list, value_list))
						err("Variable TS could not be retrieved from the json input file!");
				}
			}
		}
		


		/*=================================
		 	 section boundary parameters
		  =================================*/

		if (get_int_from_objectlist("FREE_SURF",number_readobjects,&FREE_SURF,varname_list, value_list))
			err("Variable FREE_SURF could not be retrieved from the json input file!");
		if (get_int_from_objectlist("BOUNDARY",number_readobjects,&BOUNDARY,varname_list, value_list))
			err("Variable BOUNDARY could not be retrieved from the json input file!");
		if (get_int_from_objectlist("FW",number_readobjects,&FW,varname_list, value_list))
			err("Variable FW could not be retrieved from the json input file!");
		if (get_float_from_objectlist("DAMPING",number_readobjects,&DAMPING,varname_list, value_list))
			err("Variable DAMPING could not be retrieved from the json input file!");
		if (get_float_from_objectlist("FPML",number_readobjects,&FPML,varname_list, value_list))
			err("Variable FPML could not be retrieved from the json input file!");
		if (get_float_from_objectlist("npower",number_readobjects,&npower,varname_list, value_list))
			err("Variable npower could not be retrieved from the json input file!");
		if (get_float_from_objectlist("k_max_PML",number_readobjects,&k_max_PML,varname_list, value_list))
			err("Variable k_max_PML could not be retrieved from the json input file!");
		

		/*=================================
			 section snapshot parameters
		=================================*/
		if (get_int_from_objectlist("SNAP",number_readobjects,&SNAP,varname_list, value_list)){
			SNAP=0;
			fprintf(fp,"Variable SNAP is set to default value %d.\n",SNAP);}
		else {
			if (SNAP>0) {
				if (get_int_from_objectlist("SNAP_FORMAT",number_readobjects,&SNAP_FORMAT,varname_list, value_list))
					err("Variable SNAP_FORMAT could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TSNAP1",number_readobjects,&TSNAP1,varname_list, value_list))
					err("Variable TSNAP1 could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TSNAP2",number_readobjects,&TSNAP2,varname_list, value_list))
					err("Variable TSNAP2 could not be retrieved from the json input file!");
				if (get_float_from_objectlist("TSNAPINC",number_readobjects,&TSNAPINC,varname_list, value_list))
					err("Variable TSNAPINC could not be retrieved from the json input file!");
				if (get_string_from_objectlist("SNAP_FILE",number_readobjects,SNAP_FILE,varname_list, value_list))
					err("Variable SNAP_FILE could not be retrieved from the json input file!");
			}
		}
		/* increments are read in any case, because they will be also used as increment for model output */
		if (get_int_from_objectlist("IDX",number_readobjects,&IDX,varname_list, value_list)){
			IDX=1;
			fprintf(fp,"Variable IDX is set to default value %d.\n",IDX);}
		if (get_int_from_objectlist("IDY",number_readobjects,&IDY,varname_list, value_list)){
			IDY=1;
			fprintf(fp,"Variable IDY is set to default value %d.\n",IDY);}

		/*=================================
			section seismogramm parameters
		=================================*/
		if (get_int_from_objectlist("SEISMO",number_readobjects,&SEISMO,varname_list, value_list))
			err("Variable SEISMO could not be retrieved from the json input file!");
		else {
			if (SEISMO>0){
				if ((SEISMO==1) || (SEISMO==4)) {
					if (get_string_from_objectlist("SEIS_FILE_VX",number_readobjects,SEIS_FILE_VX,varname_list, value_list))
						err("Variable SEIS_FILE_VX could not be retrieved from the json input file!");
					if (get_string_from_objectlist("SEIS_FILE_VY",number_readobjects,SEIS_FILE_VY,varname_list, value_list))
						err("Variable SEIS_FILE_VY could not be retrieved from the json input file!");
				}
				if ((SEISMO==3) || (SEISMO==4)) {
					if (get_string_from_objectlist("SEIS_FILE_DIV",number_readobjects,SEIS_FILE_DIV,varname_list, value_list))
						err("Variable SEIS_FILE_DIV could not be retrieved from the json input file!");
					if (get_string_from_objectlist("SEIS_FILE_CURL",number_readobjects,SEIS_FILE_CURL,varname_list, value_list))
						err("Variable SEIS_FILE_CURL could not be retrieved from the json input file!");
				}
				if ((SEISMO==2) || (SEISMO==4)) {
					if (get_string_from_objectlist("SEIS_FILE_P",number_readobjects,SEIS_FILE_P,varname_list, value_list))
					err("Variable SEIS_FILE_P could not be retrieved from the json input file!");
				}
				
				if (get_int_from_objectlist("READREC",number_readobjects,&READREC,varname_list, value_list))
					err("Variable READREC could not be retrieved from the json input file!");
				else {
					if (READREC==0) {
						if (get_float_from_objectlist("XREC1",number_readobjects,&XREC1,varname_list, value_list))
							err("Variable XREC1 could not be retrieved from the json input file!");
						if (get_float_from_objectlist("XREC2",number_readobjects,&XREC2,varname_list, value_list))
							err("Variable XREC2T could not be retrieved from the json input file!");
						if (get_float_from_objectlist("YREC1",number_readobjects,&YREC1,varname_list, value_list))
							err("Variable YREC1 could not be retrieved from the json input file!");
						if (get_float_from_objectlist("YREC2",number_readobjects,&YREC2,varname_list, value_list))
							err("Variable YREC2 could not be retrieved from the json input file!");
						if (get_int_from_objectlist("NGEOPH",number_readobjects,&NGEOPH,varname_list, value_list))
							err("Variable NGEOPH could not be retrieved from the json input file!");
					}
					else {
						if (get_string_from_objectlist("REC_FILE",number_readobjects,REC_FILE,varname_list, value_list))
							err("Variable REC_FILE could not be retrieved from the json input file!");
					}
				}
				if (get_int_from_objectlist("NDT",number_readobjects,&NDT,varname_list, value_list)){
					NDT=1;
					fprintf(fp,"Variable NDT is set to default value %d.\n",NDT);}
				if (get_int_from_objectlist("SEIS_FORMAT",number_readobjects,&SEIS_FORMAT,varname_list, value_list))
					err("Variable SEIS_FORMAT could not be retrieved from the json input file!");

				if (get_int_from_objectlist("REC_ARRAY",number_readobjects,&REC_ARRAY,varname_list, value_list)){
					REC_ARRAY=0;
					fprintf(fp,"Variable REC_ARRAY is set to default value %d.\n",REC_ARRAY);}
				else {
					if (REC_ARRAY>0) {
						if (get_int_from_objectlist("DRX",number_readobjects,&DRX,varname_list, value_list))
							err("Variable DRX could not be retrieved from the json input file!");
						if (get_float_from_objectlist("REC_ARRAY_DEPTH",number_readobjects,&REC_ARRAY_DEPTH,varname_list, value_list))
							err("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
						if (get_float_from_objectlist("REC_ARRAY_DIST",number_readobjects,&REC_ARRAY_DIST,varname_list, value_list))
							err("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
					}
				}
				if (get_float_from_objectlist("REFRECX",number_readobjects,&REFREC[1],varname_list, value_list))
					err("Variable REFRECX could not be retrieved from the json input file!");
				if (get_float_from_objectlist("REFRECY",number_readobjects,&REFREC[2],varname_list, value_list))
					err("Variable REFRECY could not be retrieved from the json input file!");
			}
		}

		/*=================================
			section general model and log parameters
		  =================================*/
		if (get_string_from_objectlist("MFILE",number_readobjects,MFILE,varname_list, value_list))
			err("Variable MFILE could not be retrieved from the json input file!");
		if (get_int_from_objectlist("LOG",number_readobjects,&LOG,varname_list, value_list)){
			LOG=1;
			fprintf(fp,"Variable LOG is set to default value %d.\n",LOG);}
		if (get_int_from_objectlist("CHECKPTREAD",number_readobjects,&CHECKPTREAD,varname_list, value_list)){
			CHECKPTREAD=0;
			fprintf(fp,"Variable CHECKPTREAD is set to default value %d.\n",CHECKPTREAD);}
		if (get_int_from_objectlist("CHECKPTWRITE",number_readobjects,&CHECKPTWRITE,varname_list, value_list)){
			CHECKPTWRITE=0;
			fprintf(fp,"Variable CHECKPTWRITE is set to default value %d.\n",CHECKPTWRITE);}
		if (get_string_from_objectlist("LOG_FILE",number_readobjects,LOG_FILE,varname_list, value_list)){
			sprintf(LOG_FILE,"log/LOG_FILE");}
		if (get_string_from_objectlist("CHECKPTFILE",number_readobjects,CHECKPTFILE,varname_list, value_list)){
			sprintf(CHECKPTFILE,"tmp/checkpoint_fdveps");}
		if (get_int_from_objectlist("READMOD",number_readobjects,&READMOD,varname_list, value_list))
			err("Variable READMOD could not be retrieved from the json input file!");
		

		if (get_int_from_objectlist("L",number_readobjects,&L,varname_list, value_list)){
			L=0;
			fprintf(fp,"Variable L is set to default value %d.\n",L);}
		else {
			FL=vector(1,L);
			switch(L) {
			case 0:
				break;
			case 5:
				if (get_float_from_objectlist("FL5",number_readobjects,&FL[5],varname_list, value_list))
					err("Variable FL5 could not be retrieved from the json input file!");
			case 4:
				if (get_float_from_objectlist("FL4",number_readobjects,&FL[4],varname_list, value_list))
					err("Variable FL4 could not be retrieved from the json input file!");
			case 3:
				if (get_float_from_objectlist("FL3",number_readobjects,&FL[3],varname_list, value_list))
					err("Variable FL3 could not be retrieved from the json input file!");
			case 2:
				if (get_float_from_objectlist("FL2",number_readobjects,&FL[2],varname_list, value_list))
					err("Variable FL2 could not be retrieved from the json input file!");
			case 1:
				if (get_float_from_objectlist("FL1",number_readobjects,&FL[1],varname_list, value_list))
					err("Variable FL1 could not be retrieved from the json input file!");
				break;
			default:
				err("More than four relaxation Parameter (L>5) are not implemented yet!");
				break;
			}
			if (get_float_from_objectlist("TAU",number_readobjects,&TAU,varname_list, value_list))
				err("Variable TAU could not be retrieved from the json input file!");
		}
		
		
		/*=================================
			section inversion parameters
		  =================================*/
		if (get_int_from_objectlist("INVMAT1",number_readobjects,&INVMAT1,varname_list, value_list))
					err("Variable INVMAT1 could not be retrieved from the json input file!");
		if (get_int_from_objectlist("INVMAT",number_readobjects,&INVMAT,varname_list, value_list))
			err("Variable INVMAT could not be retrieved from the json input file!");
		else {
			if (INVMAT==0) {	/* FWI is calculated */
			
				/* General inversion parameters */
				if (get_int_from_objectlist("ITERMAX",number_readobjects,&ITERMAX,varname_list, value_list))
					err("Variable ITERMAX could not be retrieved from the json input file!");
				if (get_string_from_objectlist("DATA_DIR",number_readobjects,DATA_DIR,varname_list, value_list))
					err("Variable DATA_DIR could not be retrieved from the json input file!");
				if (get_int_from_objectlist("INVTYPE",number_readobjects,&INVTYPE,varname_list, value_list)){
					INVTYPE=2;
					fprintf(fp,"\nVariable INVTYPE is set to default value %d.\n",INVTYPE);}
				if (get_int_from_objectlist("QUELLTYPB",number_readobjects,&QUELLTYPB,varname_list, value_list))
					err("Variable QUELLTYPB could not be retrieved from the json input file!");
				
				if (get_string_from_objectlist("MISFIT_LOG_FILE",number_readobjects,MISFIT_LOG_FILE,varname_list, value_list)){  	 	 
					strcpy(MISFIT_LOG_FILE,"L2_LOG.dat"); 
					fprintf(fp,"\nVariable MISFIT_LOG_FILE is set to default value %s.\n",MISFIT_LOG_FILE);}

				if (get_int_from_objectlist("TESTSHOT_START",number_readobjects,&TESTSHOT_START,varname_list, value_list))
					err("Variable TESTSHOT_START could not be retrieved from the json input file!");
				if (get_int_from_objectlist("TESTSHOT_END",number_readobjects,&TESTSHOT_END,varname_list, value_list))
					err("Variable TESTSHOT_START could not be retrieved from the json input file!");
				if (get_int_from_objectlist("TESTSHOT_INCR",number_readobjects,&TESTSHOT_INCR,varname_list, value_list))
					err("Variable TESTSHOT_INCR could not be retrieved from the json input file!");
				
				
				/* Cosine taper */
				if (get_int_from_objectlist("TAPER",number_readobjects,&TAPER,varname_list, value_list)){
					TAPER=0;
					fprintf(fp,"Variable TAPER is set to default value %d.\n",TAPER);}
				if (get_int_from_objectlist("TAPERLENGTH",number_readobjects,&TAPERLENGTH,varname_list, value_list)){
					TAPERLENGTH=10;
					fprintf(fp,"Variable TAPERLENGTH is set to default value %d.\n",TAPERLENGTH);}
				
				
				/* Definition of gradient taper geometry */
				if (get_int_from_objectlist("SWS_TAPER_GRAD_VERT",number_readobjects,&SWS_TAPER_GRAD_VERT,varname_list, value_list)){
					SWS_TAPER_GRAD_VERT=0;
					fprintf(fp,"Variable SWS_TAPER_GRAD_VERT is set to default value %d.\n",SWS_TAPER_GRAD_VERT);}
				if (get_int_from_objectlist("SWS_TAPER_GRAD_HOR",number_readobjects,&SWS_TAPER_GRAD_HOR,varname_list, value_list)){
					SWS_TAPER_GRAD_HOR=0;
					fprintf(fp,"Variable SWS_TAPER_GRAD_HOR is set to default value %d.\n",SWS_TAPER_GRAD_HOR);}
				if ((SWS_TAPER_GRAD_VERT==1) || (SWS_TAPER_GRAD_HOR==1)){
					if (get_int_from_objectlist("GRADT1",number_readobjects,&GRADT1,varname_list, value_list))
						err("Variable GRADT1 could not be retrieved from the json input file!");
					if (get_int_from_objectlist("GRADT2",number_readobjects,&GRADT2,varname_list, value_list))
						err("Variable GRADT2 could not be retrieved from the json input file!");
					if (get_int_from_objectlist("GRADT3",number_readobjects,&GRADT3,varname_list, value_list))
						err("Variable GRADT3 could not be retrieved from the json input file!");
					if (get_int_from_objectlist("GRADT4",number_readobjects,&GRADT4,varname_list, value_list))
						err("Variable GRADT4 could not be retrieved from the json input file!");
				}
				if (get_int_from_objectlist("SWS_TAPER_GRAD_SOURCES",number_readobjects,&SWS_TAPER_GRAD_SOURCES,varname_list, value_list)){
					SWS_TAPER_GRAD_SOURCES=0;
					fprintf(fp,"Variable SWS_TAPER_GRAD_SOURCES is set to default value %d.\n",SWS_TAPER_GRAD_SOURCES);}
				if (get_int_from_objectlist("SWS_TAPER_CIRCULAR_PER_SHOT",number_readobjects,&SWS_TAPER_CIRCULAR_PER_SHOT,varname_list, value_list)){
					SWS_TAPER_CIRCULAR_PER_SHOT=0;
					fprintf(fp,"Variable SWS_TAPER_CIRCULAR_PER_SHOT is set to default value %d.\n",SWS_TAPER_CIRCULAR_PER_SHOT);}
				if ((SWS_TAPER_GRAD_SOURCES==1) || (SWS_TAPER_CIRCULAR_PER_SHOT==1)){	
					if (get_int_from_objectlist("SRTSHAPE",number_readobjects,&SRTSHAPE,varname_list, value_list))
						err("Variable SRTSHAPE could not be retrieved from the json input file!");
					if (get_float_from_objectlist("SRTRADIUS",number_readobjects,&SRTRADIUS,varname_list, value_list))
						err("Variable SRTRADIUS could not be retrieved from the json input file!");
					if (get_int_from_objectlist("FILTSIZE",number_readobjects,&FILTSIZE,varname_list, value_list))
						err("Variable FILTSIZE could not be retrieved from the json input file!");
				}
				if (get_int_from_objectlist("SWS_TAPER_FILE",number_readobjects,&SWS_TAPER_FILE,varname_list, value_list)){
					SWS_TAPER_FILE=0;
					fprintf(fp,"Variable SWS_TAPER_FILE is set to default value %d.\n",SWS_TAPER_FILE);}
					
				
				/* Definition of smoothing (spatial filtering) of the gradients */
				if (get_int_from_objectlist("SPATFILTER",number_readobjects,&SPATFILTER,varname_list, value_list)){
					SPATFILTER=0;
					fprintf(fp,"Variable SPATFILTER is set to default value %d.\n",SPATFILTER);}
					
				else {
					if (SPATFILTER==1) {
						if (get_int_from_objectlist("SPAT_FILT_SIZE",number_readobjects,&SPAT_FILT_SIZE,varname_list, value_list))
							err("Variable SPAT_FILT_SIZE could not be retrieved from the json input file!");
						if (get_int_from_objectlist("SPAT_FILT_1",number_readobjects,&SPAT_FILT_1,varname_list, value_list))
							err("Variable SPAT_FILT_1 could not be retrieved from the json input file!");
						if (get_int_from_objectlist("SPAT_FILT_ITER",number_readobjects,&SPAT_FILT_ITER,varname_list, value_list))
							err("Variable SPAT_FILT_ITER could not be retrieved from the json input file!");
					}
				}
				
				/* Definition of 2D-Gaussian filter of the gradients */
				if (get_int_from_objectlist("GRAD_FILTER",number_readobjects,&GRAD_FILTER,varname_list, value_list)){
					GRAD_FILTER=0;
					fprintf(fp,"Variable GRAD_FILTER is set to default value %d.\n",GRAD_FILTER);}
					
				if (get_int_from_objectlist("FILT_SIZE_GRAD",number_readobjects,&FILT_SIZE_GRAD,varname_list, value_list)){
					FILT_SIZE_GRAD=0;
					fprintf(fp,"Variable FILT_SIZE_GRAD is set to default value %d.\n",FILT_SIZE_GRAD);}
				
				
				/* Output of inverted models */
				if (get_string_from_objectlist("INV_MODELFILE",number_readobjects,INV_MODELFILE,varname_list, value_list))
					err("Variable INV_MODELFILE could not be retrieved from the json input file!");
				if (get_int_from_objectlist("nfstart",number_readobjects,&nfstart,varname_list, value_list))
					err("Variable nfstart could not be retrieved from the json input file!");
				if (get_int_from_objectlist("nf",number_readobjects,&nf,varname_list, value_list))
					err("Variable nf could not be retrieved from the json input file!");
				
				
				/* Output of gradients */
				if (get_string_from_objectlist("JACOBIAN",number_readobjects,JACOBIAN,varname_list, value_list))
					err("Variable JACOBIAN could not be retrieved from the json input file!");
				if (get_int_from_objectlist("nfstart_jac",number_readobjects,&nfstart_jac,varname_list, value_list))
					err("Variable nfstart_jac could not be retrieved from the json input file!");
				if (get_int_from_objectlist("nf_jac",number_readobjects,&nf_jac,varname_list, value_list))
					err("Variable nf_jac could not be retrieved from the json input file!");
				
				
				/* Inversion for density */
				if (get_int_from_objectlist("INV_RHO_ITER",number_readobjects,&INV_RHO_ITER,varname_list, value_list)){
					INV_RHO_ITER=0;	
					fprintf(fp,"Variable INV_RHO_ITER is set to default value %d.\n",INV_RHO_ITER);}
					
				/* Inversion for Vp */
				if (get_int_from_objectlist("INV_VP_ITER",number_readobjects,&INV_VP_ITER,varname_list, value_list)){
					INV_VP_ITER=0;	
					fprintf(fp,"Variable INV_VP_ITER is set to default value %d.\n",INV_VP_ITER);}
					
				/* Inversion for Vs */
				if (get_int_from_objectlist("INV_VS_ITER",number_readobjects,&INV_VS_ITER,varname_list, value_list)){
					INV_VS_ITER=0;	
					fprintf(fp,"Variable INV_VS_ITER is set to default value %d.\n",INV_VS_ITER);}
				
				
				/* Upper and lower limits for model parameters */
				if (get_float_from_objectlist("VPUPPERLIM",number_readobjects,&VPUPPERLIM,varname_list, value_list)){
					VPUPPERLIM=25000.0;
					fprintf(fp,"Variable VPUPPERLIM is set to default value %f.\n",VPUPPERLIM);}
				if (get_float_from_objectlist("VPLOWERLIM",number_readobjects,&VPLOWERLIM,varname_list, value_list)){
					VPLOWERLIM=0.0;
					fprintf(fp,"Variable VPLOWERLIM is set to default value %f.\n",VPLOWERLIM);}
				if (get_float_from_objectlist("VSUPPERLIM",number_readobjects,&VSUPPERLIM,varname_list, value_list)){
					VSUPPERLIM=25000.0;
					fprintf(fp,"Variable VSUPPERLIM is set to default value %f.\n",VSUPPERLIM);}
				if (get_float_from_objectlist("VSLOWERLIM",number_readobjects,&VSLOWERLIM,varname_list, value_list)){
					VSLOWERLIM=0.0;
					fprintf(fp,"Variable VSLOWERLIM is set to default value %f.\n",VSLOWERLIM);}
				if (get_float_from_objectlist("RHOUPPERLIM",number_readobjects,&RHOUPPERLIM,varname_list, value_list)){
					RHOUPPERLIM=25000.0;
					fprintf(fp,"Variable RHOUPPERLIM is set to default value %f.\n",RHOUPPERLIM);}
				if (get_float_from_objectlist("RHOLOWERLIM",number_readobjects,&RHOLOWERLIM,varname_list, value_list)){
					RHOLOWERLIM=0.0;
					fprintf(fp,"Variable RHOLOWERLIM is set to default value %f.\n",RHOLOWERLIM);}
					
				
				/* Hessian and Gradient-Method */
				if (get_int_from_objectlist("HESSIAN",number_readobjects,&HESSIAN,varname_list, value_list)){
					HESSIAN=0;
					fprintf(fp,"Variable HESSIAN is set to default value %d.\n",HESSIAN);}
				if (get_int_from_objectlist("GRAD_METHOD",number_readobjects,&GRAD_METHOD,varname_list, value_list))
					err("Variable GRAD_METHOD could not be retrieved from the json input file!");
				if (HESSIAN){
					if (get_float_from_objectlist("FC_HESSIAN",number_readobjects,&FC_HESSIAN,varname_list, value_list))
						err("Variable FC_HESSIAN could not be retrieved from the json input file!");
					if (get_int_from_objectlist("ORDER_HESSIAN",number_readobjects,&ORDER_HESSIAN,varname_list, value_list))
						err("Variable ORDER_HESSIAN could not be retrieved from the json input file!");
				}
				
				
				/* Definition of smoothing the models vp and vs */
				if (get_int_from_objectlist("MODEL_FILTER",number_readobjects,&MODEL_FILTER,varname_list, value_list)){
					MODEL_FILTER=0;
					fprintf(fp,"Variable MODEL_FILTER is set to default value %d.\n",MODEL_FILTER);}
				else {
					if (MODEL_FILTER==1) {
						if (get_int_from_objectlist("FILT_SIZE",number_readobjects,&FILT_SIZE,varname_list, value_list))
							err("Variable FILT_SIZE could not be retrieved from the json input file!");
					}
				}
				
				
				/* Definition of inversion for source time function */
				if (get_int_from_objectlist("INV_STF",number_readobjects,&INV_STF,varname_list, value_list)){
					INV_STF=0;
					fprintf(fp,"Variable INV_STF is set to default value %d.\n",INV_STF);}
				else {
					if (INV_STF==1) {
						if (get_string_from_objectlist("PARA",number_readobjects,PARA,varname_list, value_list))
							err("Variable PARA could not be retrieved from the json input file!");	
						if (get_int_from_objectlist("N_STF",number_readobjects,&N_STF,varname_list, value_list))
							err("Variable N_STF could not be retrieved from the json input file!");
						if (get_int_from_objectlist("N_STF_START",number_readobjects,&N_STF_START,varname_list, value_list))
							err("Variable N_STF_START could not be retrieved from the json input file!");	
					}
				}
				
				
				/* Frequency filtering during inversion */
				if (get_int_from_objectlist("TIME_FILT",number_readobjects,&TIME_FILT,varname_list, value_list)){
					TIME_FILT=0;
					fprintf(fp,"Variable TIME_FILT is set to default value %d.\n",TIME_FILT);}
				else {
					if (TIME_FILT==1) {
						if (get_float_from_objectlist("FC_START",number_readobjects,&FC_START,varname_list, value_list))
							err("Variable FC_START could not be retrieved from the json input file!");
						if (get_float_from_objectlist("FC_END",number_readobjects,&FC_END,varname_list, value_list))
							err("Variable FC_END could not be retrieved from the json input file!");
						if (get_float_from_objectlist("FC_INCR",number_readobjects,&FC_INCR,varname_list, value_list))
							err("Variable FC_INCR could not be retrieved from the json input file!");
						if (get_int_from_objectlist("ORDER",number_readobjects,&ORDER,varname_list, value_list))
							err("Variable ORDER could not be retrieved from the json input file!");
					}
				}
				
				
				/* Gradient calculation */
				if (get_int_from_objectlist("LNORM",number_readobjects,&LNORM,varname_list, value_list))
					err("Variable LNORM could not be retrieved from the json input file!");
				if (get_int_from_objectlist("NORMALIZE",number_readobjects,&NORMALIZE,varname_list, value_list))
					err("Variable NORMALIZE could not be retrieved from the json input file!");
				if (get_int_from_objectlist("DTINV",number_readobjects,&DTINV,varname_list, value_list))
					err("Variable DTINV could not be retrieved from the json input file!");
				
				
				/* Step length estimation */
				if (get_float_from_objectlist("EPS_SCALE",number_readobjects,&EPS_SCALE,varname_list, value_list))
					err("Variable EPS_SCALE could not be retrieved from the json input file!");
				if (get_int_from_objectlist("STEPMAX",number_readobjects,&STEPMAX,varname_list, value_list))
					err("Variable STEPMAX could not be retrieved from the json input file!");
				if (get_float_from_objectlist("SCALEFAC",number_readobjects,&SCALEFAC,varname_list, value_list))
						err("Variable SCALEFAC could not be retrieved from the json input file!");	
				
				
				/* Termination of the program */
				if (get_float_from_objectlist("PRO",number_readobjects,&PRO,varname_list, value_list))
					err("Variable PRO could not be retrieved from the json input file!");
				if (get_int_from_objectlist("MIN_ITER",number_readobjects,&MIN_ITER,varname_list, value_list)){
					MIN_ITER=0;
					fprintf(fp,"Variable MIN_ITER is set to default value %d.\n",MIN_ITER);}	
							
				
				/* Trace killing */
				if (get_int_from_objectlist("TRKILL",number_readobjects,&TRKILL,varname_list, value_list)){
					TRKILL=0;
					fprintf(fp,"Variable TRKILL is set to default value %d.\n",TRKILL);}
				else {
					if (TRKILL==1) {
						if (get_string_from_objectlist("TRKILL_FILE",number_readobjects,TRKILL_FILE,varname_list, value_list))
							err("Variable TRKILL_FILE could not be retrieved from the json input file!");
					}
				}
				
				
				/* Time windowing and damping */
				if (get_int_from_objectlist("TIMEWIN",number_readobjects,&TIMEWIN,varname_list, value_list)){
					TIMEWIN=0;
					fprintf(fp,"Variable TIMEWIN is set to default value %d.\n",TIMEWIN);}
				else {
					if (TIMEWIN==1){
						if (get_string_from_objectlist("PICKS_FILE",number_readobjects,PICKS_FILE,varname_list, value_list))
							err("Variable PICKS_FILE could not be retrieved from the json input file!");
						if (get_float_from_objectlist("TWLENGTH_PLUS",number_readobjects,&TWLENGTH_PLUS,varname_list, value_list))
							err("Variable TWLENGTH_PLUS could not be retrieved from the json input file!");
						if (get_float_from_objectlist("TWLENGTH_MINUS",number_readobjects,&TWLENGTH_MINUS,varname_list, value_list))
							err("Variable TWLENGTH_MINUS could not be retrieved from the json input file!");
						if (get_float_from_objectlist("GAMMA",number_readobjects,&GAMMA,varname_list, value_list))
							err("Variable GAMMA could not be retrieved from the json input file!");
					}
				}
			} /* end if (INVMAT==0) */
			else {
				if (INVMAT==10){	/* only forward modeling is applied */
					ITERMAX=1;
					strcpy(INV_MODELFILE,MFILE);
					DTINV=1;
					INVTYPE=2;
					fprintf(fp,"Setting some default values needed for modeling:\n");
					fprintf(fp,"Variable ITERMAX is set to default value %d.\n",ITERMAX);
					fprintf(fp,"Variable INV_MODELFILE is set to default value %s.\n",MFILE);
					fprintf(fp,"Variable DTINV is set to default value %d.\n",DTINV);
					fprintf(fp,"Variable INVTYPE is set to default value %d.\n",INVTYPE);
					}
			}
					
		}
		
		fprintf(fp,"\nEnd of setting default values\n");
		fprintf(fp,"=====================================\n\n");
				
				
		/********************************************/
		/* Check files and directories if necessary */
		/********************************************/

		/* signal file */
		if (QUELLART == 3)
		{
			if (access(SIGNAL_FILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The source file does not exist!\n");
				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(SIGNAL_FILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The signal file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}

		/* source file */
		if (SRCREC==1)
		{
			if (access(SOURCE_FILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The source file does not exist!\n");
				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(SOURCE_FILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The source file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}

		/* model file -> this is checked in more detail in readmod.c!
		if (READMOD)
		{
			if (access(MFILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The model file does not exist!\n");
				fprintf(fp, "        File name: <%s>", MFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(MFILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The model file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", MFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}*/

		/* receiver file */
		if (READREC)
		{
			if (access(REC_FILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The receiver file does not exist!\n");
				fprintf(fp, "        File name: <%s>", REC_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(REC_FILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The receiver file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", REC_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}

		/* checkpoint file */
		if (CHECKPTREAD || CHECKPTWRITE)
		{
			if (access(CHECKPTFILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The checkpoint file does not exist!\n");
				fprintf(fp, "        File name: <%s>", CHECKPTFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(CHECKPTFILE,6) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The checkpoint file does not have read and/or write access!\n");
				fprintf(fp, "        File name: <%s>", CHECKPTFILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}
		
		/* trace kill file */
		if (TRKILL == 1)
		{
			if (access(TRKILL_FILE,0) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The trace kill file does not exist!\n");
				fprintf(fp, "        File name: <%s>", TRKILL_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
			else if (access(TRKILL_FILE,4) != 0)
			{
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The trace kill file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", TRKILL_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}
		
		/********************************************/
		/* ERROR                                    */
		/********************************************/
		if (fserr)
		{
			fprintf(fp, "\n");
			sprintf(errormessage, "\n  in: <read_par_json.c> \n");
			err(errormessage);
		}


	} /* End of if(MYID==0) */
}
