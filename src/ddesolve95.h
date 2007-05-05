#include <stdio.h>
#include <R.h>
#include <Rdefines.h>


typedef struct { 
	int no_var, no_otherVars;
	int nhv,nlag,nsw;
	double dt,dout,t0,t1,tol;
	long hbsize;
	char **cname,*initialtext,*initialtitle,**cinfo;
	FILE *file;
	int quit,newrun,cont,*findex,fileno;
	double **vals, *tmp_other_vals;
	int vals_size, vals_ind;
	double current_t;
} globaldatatype;

typedef struct {
	SEXP env;
	SEXP gradFunc;
	SEXP yinit;
	SEXP parms;
	int useParms;
	int gradFuncListReturn;
} globalRdatatype;


globaldatatype data;
globalRdatatype r_stuff;
