#include <R.h>
#include <Rdefines.h>

#include <stdlib.h>
#include "ddesolve95.h"
#include "ddeq.h"

#define CH_BUF_SIZE 128

int the_test_phase=0;

/*===========================================================================*/
void PBSerror(char *str)
{
	error(str);
	exit(-1);
}

/*===========================================================================*/
void info(char *str)
{
	warning(str);
}

/*===========================================================================*/
void output(double *s,double t)
{
	/*data.vals[0] is for t,
	  and [1..(no_var+1)] are reserved for s[0..no_var] vars
	  */
	int i;
	data.vals[0][data.vals_ind] = t;
	for(i=0;i<data.no_var;i++)
		data.vals[i+1][data.vals_ind] = s[i];
	for(i=0;i<data.no_otherVars;i++)
		data.vals[1+data.no_var+i][data.vals_ind] = data.tmp_other_vals[i];
	
	data.vals_ind++;
	
	if (data.vals_ind >= data.vals_size) {
		for(i=0;i<(1+data.no_var+data.no_otherVars);i++) {
			data.vals[i] = (double*)realloc(data.vals[i], sizeof(double)*2*data.vals_size);
			if (data.vals[i]==NULL)
				error("memory (re)allocation failed");
		}
		data.vals_size *= 2;
	}
}

/*===========================================================================*/
/* cont is zero for fresh run, and 1 for continuation */
void numerics(double *c,int cont)
{ 
  static double *s;
  double t0,t1,dt,dout;
  int ns,nsw,nhv,nlag,reset=1,fixstep=0;
  static int first=1;
  long hbsize;
  ns=data.no_var;nsw=data.nsw;nhv=data.nhv;nlag=data.nlag;
  t0=data.t0;t1=data.t1;dt=data.dt;dout=data.dout;hbsize=data.hbsize;
  if (cont) reset=0;
  if (!cont)
  { if (!first) { free(s);first=0;}
    s=(double *)calloc(data.no_var,sizeof(double));
    ddeinitstate(s,c,t0);
  }
  dde(s,c,t0,t1,&dt,data.tol,dout,ns,nsw,nhv,hbsize,nlag,reset,fixstep);
  data.dt=dt;
}

/*===========================================================================*/
void setupglobaldata(int no_vars, int no_otherVars, double *settings)
{ 
  int i;

  data.tol=settings[0];
  data.t0=settings[1];        /* start time */
  data.t1=settings[2];        /* stop time */
  
  data.dt=settings[3];        /* initial timestep */
  data.dout = settings[4];    /* approximate average output timestep */

  data.hbsize=settings[5];    /* how many past values to store for each history variable */
  data.no_var=no_vars;
  
  data.no_otherVars=no_otherVars;

  data.nsw=0;          /* number of switch varaibles */  
  data.nhv=no_vars;         /* Number of history (lagged) variables */
  data.nlag=1;        /* Number of lag markers per history variable (set to 1 if unsure)*/
  

	data.vals_size=1000; /* size will grow, this is just initial min size */
	data.vals_ind=0;
	data.vals = (double**)malloc(sizeof(double*)*(data.no_var+no_otherVars+1));
	if (data.vals==NULL)
		error("memory allocation failed");
	for(i=0;i<(data.no_var+no_otherVars+1);i++) {
		data.vals[i]=(double*)malloc(sizeof(double)*data.vals_size);
		if (data.vals[i]==NULL)
			error("memory allocation failed");
	}
	if (data.no_otherVars>0) {
		data.tmp_other_vals = (double*)malloc(sizeof(double)*data.no_otherVars);
		if (data.tmp_other_vals==NULL) {
			error("memory allocation failed");
		}
	} else {
		data.tmp_other_vals=NULL;
	}
}

/*===========================================================================*/
void freeglobaldata()
{
	int i;
	if (data.vals) {
		for(i=0;i<(data.no_var+data.no_otherVars+1);i++) {
			free(data.vals[i]);
		}
		free(data.vals);
		data.vals=NULL;
	}
	if (data.tmp_other_vals) {
		free(data.tmp_other_vals);
		data.tmp_other_vals=NULL;
	}
}

int testFunc(int no_var, double *test_vars, double t, SEXP *names, PROTECT_INDEX *names_ipx)
{
	SEXP fcall, p1, p2, result;
	int len;
	
	/* argument 1 `t' */
	PROTECT(p1=NEW_NUMERIC(1));
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));
    
	/* argument 2 `s' */
	PROTECT(p2=NEW_NUMERIC(no_var));
	memcpy(NUMERIC_POINTER(p2), test_vars, no_var*sizeof(double));
	
	/* call R user function */
	if (r_stuff.useParms)
		PROTECT(fcall = lang4(r_stuff.gradFunc, p1, p2, r_stuff.parms));
	else
    	PROTECT(fcall = lang3(r_stuff.gradFunc, p1, p2));
    PROTECT(result = eval(fcall, r_stuff.env));
    
	if (isReal(result)) {
		r_stuff.gradFuncListReturn=0;
		if (LENGTH(result)!=no_var)
			error("func return error: length of vector (%i) does not match that of initial y values (%i)\n", LENGTH(result), no_var);
		UNPROTECT(4);
    	return 0;
    } else if (isVector(result)) {
		r_stuff.gradFuncListReturn=1;
	} else {
		error("func return error: should return numeric vector or list. (got type \"%i\")\n", TYPEOF(result));
	}
	
	p1 = VECTOR_ELT(result, 0);
	p2 = VECTOR_ELT(result, 1);
	
	if (LENGTH(result)!=2)
		error("func return error: returned list should have length two\n", TYPEOF(p1));
	if (!isReal(p1))
		error("func return error: first element of list should be numeric. (got type \"%i\")\n", TYPEOF(p1));
	if (!isReal(p2) && !isNull(p2))
		error("func return error: second element of list should be numeric or NULL. (got type \"%i\")\n", TYPEOF(p2));
	if (LENGTH(p1)!=no_var)
		error("func return error: length of first element vector (%i) does not match that of initial y values (%i)\n", LENGTH(p1), no_var);

	if (isNull(p2))
		len = 0;
	else {
		len = LENGTH(p2);
		/* return names - which stays protected even after return */
		REPROTECT((*names) = GET_NAMES(p2), *names_ipx);
	}
	UNPROTECT(4);
	return(len);
}

/*===========================================================================*/
SEXP startDDE(SEXP gradFunc, SEXP env, SEXP yinit, SEXP parms, SEXP settings)
{
	SEXP list, vect, extra_names, yinit_names, names;
	PROTECT_INDEX extra_names_ipx;
	double *p;
	int i,j, no_var, no_otherVar;
	char ch_buf[CH_BUF_SIZE];
	
	/* save R global data for later */
	if(!isFunction(gradFunc)) error("‘gradFunc’ must be a function");
	if(!isEnvironment(env)) error("‘env’ should be an environment");
	if(!isNumeric(yinit)) error("‘yinit’ should be a numeric vector");
	if(!isNumeric(settings)) error("‘settings’ should be a numeric vector");
	
	r_stuff.env = env;
	r_stuff.gradFunc = gradFunc;
	r_stuff.parms = parms;
	r_stuff.yinit = yinit;
	
	/* check if supplied function is func(y,t) or func(y,t,parms) */
	list = FORMALS(gradFunc);
	i=0;
	while (list != R_NilValue) {
		i++;
		list = CDR(list);
	}
	if (i!=2 && i!=3) error("‘gradFunc’ must be in the form func(y,t) or func(y,t,parms)");
	r_stuff.useParms = (i==3); /* only use parms if 3 arguments */
	
	no_var = LENGTH(yinit);
	
	/* test the supplied `func` and get number of other vars 
	   ***must set the_test_phase to avoid errors in pastvalue and pastgradient */
	the_test_phase=1;
	PROTECT_WITH_INDEX(extra_names=R_NilValue, &extra_names_ipx);
	no_otherVar=testFunc(LENGTH(yinit), NUMERIC_POINTER(yinit), NUMERIC_POINTER(settings)[1], 
	                     &extra_names, &extra_names_ipx);
	the_test_phase=0;

	
	/* print names to see if it works */
	PROTECT(names = allocVector(STRSXP, no_otherVar + no_var + 1));
	yinit_names = GET_NAMES(yinit);
	
	/* fill up names as <- c("t", names(yinit), names(extra_returned_vector)) */
	SET_STRING_ELT(names, 0, mkChar("t"));
	for(i=0;i<no_var;i++) {
		if (isNull(yinit_names)) {
			/* yuck - possible buffer overflow */
			sprintf(ch_buf, "y%i", i+1);
			SET_STRING_ELT(names, i+1, mkChar(ch_buf));
		} else {
			SET_STRING_ELT(names, i+1, STRING_ELT(yinit_names, i));
		}
	}
	for(i=0;i<no_otherVar;i++) {
		if (isNull(extra_names)) {
			/* yuck - possible buffer overflow */
			sprintf(ch_buf, "extra%i", i+1);
			SET_STRING_ELT(names, i+no_var+1, mkChar(ch_buf));
		} else {
			SET_STRING_ELT(names, i+no_var+1, STRING_ELT(extra_names, i));
		}
	}

	
	setupglobaldata(LENGTH(yinit), no_otherVar, NUMERIC_POINTER(settings));
	
	/* preform dde calculations */
	numerics(NUMERIC_POINTER(yinit), 0);
	
	/* create list which will be base of polyset data.frame */
	PROTECT(list=allocVector(VECSXP, data.no_var+data.no_otherVars+1));
	
	/* room for all data (Y) AND time (T) */
	for(j=0;j<(data.no_var+data.no_otherVars+1);j++) {
		/* create numeric vector */
		PROTECT(vect=NEW_NUMERIC(data.vals_ind));

		/* and fill it up */
		p=NUMERIC_POINTER(vect);
		for(i=0;i<data.vals_ind;i++)
			p[i]=data.vals[j][i];

		SET_VECTOR_ELT(list, j, vect);
		UNPROTECT(1);
	}
	
	setAttrib(list, R_NamesSymbol, names);
	
	freeglobaldata();
	UNPROTECT(3);
	return list;
}

/*===========================================================================*/
SEXP getPastValue(SEXP t, SEXP markno)
{
	SEXP value;
	int i;

	/* return some value for `func`'s test phase
	   --this won't be used during integration */
	if (the_test_phase)
		return r_stuff.yinit;
	
	if (data.vals==NULL) error("pastvalue can only be called from `func` when triggered by dde solver.");
	if (!isNumeric(t)) error("‘t’ should be numeric");
	if (!isInteger(markno)) error("‘markno’ must be an integer");
	if (data.hbsize<=0) error("no history buffer was created. dde(...) "
	                          "should be called with hbsize>0");
	if (INTEGER_POINTER(markno)[0] >= data.nlag || INTEGER_POINTER(markno)[0] < 0) 
		error("markno is out of bounds and should be in 0..data.nlag");
	
	if (NUMERIC_POINTER(t)[0] < data.t0 || NUMERIC_POINTER(t)[0] >= data.current_t) {
		Rprintf("getvalue error: t=%.5f current integration at t=%.5f\n", NUMERIC_POINTER(t)[0], data.current_t);
		error("t is out of bounds and should be >= t0 and < t.\nCalling pastvalue(t) is not allowed.");
	}

	PROTECT(value=NEW_NUMERIC(data.no_var));
	for(i=0;i<data.no_var;i++) {
		NUMERIC_POINTER(value)[i] = pastvalue(i,
		                                      *NUMERIC_POINTER(t),
		                                      *INTEGER_POINTER(markno));
	}
	UNPROTECT(1);
	return(value);
}

/*===========================================================================*/
SEXP getPastGradient(SEXP t, SEXP markno)
{
	SEXP value;
	int i;

	/* return some value for `func`'s test phase
	   --this won't be used during integration */
	if (the_test_phase)
		return r_stuff.yinit;

	if (data.vals==NULL) error("pastgradient can only be called from `func` when triggered by dde solver.");
	if (!isNumeric(t)) error("‘t’ should be numeric");
	if (!isInteger(markno)) error("‘markno’ must be an integer");
	if (data.hbsize<=0) error("no history buffer was created. dde(...) "
	                          "should be called with hbsize>0");
	if (INTEGER_POINTER(markno)[0] >= data.nlag || INTEGER_POINTER(markno)[0] < 0) 
		error("markno is out of bounds and should be in 0..data.nlag");
	
	if (NUMERIC_POINTER(t)[0] < data.t0 || NUMERIC_POINTER(t)[0] >= data.current_t) 
		error("t is out of bounds and should be >= t0 and < t.\nCalling pastvalue(t) is not allowed.");

	PROTECT(value=NEW_NUMERIC(data.no_var));
	for(i=0;i<data.no_var;i++) {
		NUMERIC_POINTER(value)[i] = pastgradient(i,
		                                         *NUMERIC_POINTER(t),
		                                         *INTEGER_POINTER(markno));
	}
	UNPROTECT(1);
	return(value);
}

