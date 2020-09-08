#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include "which.h"
#include "cubic_spline.h"
#include "linear_interpolation.h"

double	wage_f(double agg_cap, double agg_lab, double tfp, double alpha, double delta);
double	rent_f(double agg_cap, double agg_lab, double tfp, double alpha, double delta);
double	marg_util(double cons, double risk_aversion);
void	markov(double *tran, int ng, int ns, int *r);

SEXP	getvar(const unsigned char *name, SEXP env);
SEXP	getListElement(const unsigned char *name, SEXP list);
int	toInt(SEXP a);

#define epsilon 1e-8
#define idx2(e,z) (e + ne*(z))
#define idx3(e,k,z) (e + ne*(k) + ne*nk*(z))
#define idx4(a,e,k,z) (a + na*(e) + na*ne*(k) + na*ne*nk*(z))
#define idxTM(e0,z0,e1,z1) ((e0) + ne*(z0) + ne*nz*(e1) + ne*nz*ne*(z1))

#define DEBUG 0

#if DEBUG==1
FILE	*fileopen(const char *fname, int debug_c){
	FILE *stream;
	char	buf[1024];

	if(NULL==(stream=fopen(fname, 0==debug_c ? "w": "a"))){
		sprintf(buf, "Coudln't open the file: %s\n", fname);
		error(buf);
	}
	return stream;
}
#endif

SEXP	RHS_Euler(SEXP asset, SEXP e_state, SEXP K_agg, SEXP L_agg, SEXP z_state, SEXP consf, SEXP xgrid, SEXP env)
{
	int	alen, glen, i, kloc, na, ne, nk, nz, z, e;
	double	weight, coh, cons, cons0, cons1, alpha, beta, gamma, delta, rent, wage;
	CubicSpline	**cs;
	SEXP	rv1;
#if DEBUG==1
FILE	*stream, *stream2, *stream3, *stream4;
static	int debug_c = 0, k;
#endif
SEXP	gparam, dparam, kgrid, trans_prob, trans_ez, endow, tfp;

/* Rのグローバル領域の変数を引っ張る */
	dparam = getvar("dparam", env);
	gparam = getvar("gparam", env);
	kgrid = getvar("kgrid", env);

	na = toInt(getListElement("na", gparam));
	ne = toInt(getListElement("ne", dparam));
	nk = toInt(getListElement("nk", gparam));
	nz = toInt(getListElement("nz", dparam));

	alpha = REAL(getListElement("alpha", dparam))[0];
	beta = REAL(getListElement("beta", dparam))[0];
	gamma = REAL(getListElement("gamma", dparam))[0];
	delta = REAL(getListElement("delta", dparam))[0];

	trans_prob = getvar("trans_prob", env);
	trans_ez = getListElement("ez", trans_prob);

	endow = getListElement("endow", dparam);
	tfp = getListElement("tfp", dparam);

	glen = length(kgrid);
	alen = length(asset);

	rv1 = PROTECT(allocVector(REALSXP, alen)); /* 計算結果を保持するベクターを作成 */
	cs = Calloc(ne*nk*nz, CubicSpline*); /* スプライン補間用のデータのキャッシング領域 */
	for(i=0; i<ne*nk*nz; i++)
		cs[i] = NULL;

#if DEBUG==1
	stream = fileopen("xgrid-consf_R.txt", debug_c);
	stream2 = fileopen("cons01_R.txt", debug_c);
	stream3 = fileopen("kloc_R.txt", debug_c);
	stream4 = fileopen("rent-wage.txt", debug_c);

	if(1*4*504 <= debug_c && debug_c < 2*4*504){
		for(z=0; z<nz; z++){
			for(k=0; k<nk; k++){
				for(e=0; e<ne; e++){
					for(i=0; i<na; i++){
						fprintf(stream, "%d\t%d\t%d\t%d\t%f\t%f\n", z, k, e, i, REAL(xgrid)[idx4(i, e, k, z)], REAL(consf)[idx4(i, e, k, z)]);
					}
				}
			}
		}
	}
#endif

	for(i=0; i<alen; i++){
		REAL(rv1)[i] = 0;

		kloc = which(REAL(kgrid), glen, REAL(K_agg)[i]);
		if(glen <= kloc + 1)
			kloc = glen - 2;

#if DEBUG==1
		if(1*4*504 <= debug_c && debug_c < 2*4*504){
			fprintf(stream3, "%d\n", kloc + 1);
		}
#endif

		/* use linear interpolation over aggregate capital grid */
		weight = (REAL(K_agg)[i] - REAL(kgrid)[kloc])/(REAL(kgrid)[kloc + 1] - REAL(kgrid)[kloc]);
		if(weight < 0)
			weight = 0;
		else if(weight > 1)
			weight = 1;

		for(z=0; z<nz; z++){ /* next period's TFP level */

			rent = rent_f(REAL(K_agg)[i], REAL(L_agg)[z], REAL(tfp)[z], alpha, delta);
			wage = wage_f(REAL(K_agg)[i], REAL(L_agg)[z], REAL(tfp)[z], alpha, delta);

#if DEBUG==1
			if(0*4*504 <= debug_c && debug_c < 1*4*504){
				fprintf(stream4, "%d\t%f\t%f\n", z, rent, wage);
			}
#endif

			for(e=0; e<ne; e++){ /* next period's employment status */
				coh = wage*REAL(endow)[e] + (1.0 + rent)*REAL(asset)[i];

				/* lower capital grid */
				if (coh < REAL(xgrid)[idx4(0, e, kloc, z)]){
					cons0 = coh;
				} else if (coh > REAL(xgrid)[idx4(na - 1, e, kloc, z)]) { /* スプライン関数の振る舞いの違いから、おそらく不要 */
					LinearInterpolation(&REAL(xgrid)[idx4(0, e, kloc, z)], &REAL(consf)[idx4(0, e, kloc, z)], na, &coh, &cons0, 1);
				} else {
					if(NULL == cs[idx3(e, kloc, z)])
						cs[idx3(e, kloc, z)] = CubicSplineSetup(&REAL(xgrid)[idx4(0, e, kloc, z)], &REAL(consf)[idx4(0, e, kloc, z)], na, BoundaryConditionNatural, 0, BoundaryConditionNatural, 0);
					CubicSplineInterpolation(cs[idx3(e, kloc, z)], &coh, &cons0, 1);
				}

				/* upper capital grid */
				if (coh < REAL(xgrid)[idx4(0, e, kloc + 1, z)]){
					cons1 = coh;
				} else if (coh > REAL(xgrid)[idx4(na - 1, e, kloc + 1, z)]) { /* スプライン関数の振る舞いの違いから、おそらく不要 */
					LinearInterpolation(&REAL(xgrid)[idx4(0, e, kloc + 1, z)], &REAL(consf)[idx4(0, e, kloc + 1, z)], na, &coh, &cons1, 1);
				} else {
					if(NULL == cs[idx3(e, kloc + 1, z)])
						cs[idx3(e, kloc + 1, z)] = CubicSplineSetup(&REAL(xgrid)[idx4(0, e, kloc + 1, z)], &REAL(consf)[idx4(0, e, kloc + 1, z)], na, BoundaryConditionNatural, 0, BoundaryConditionNatural, 0);
					CubicSplineInterpolation(cs[idx3(e, kloc + 1, z)], &coh, &cons1, 1);
				}

				cons = (1.0 - weight)*cons0 + weight*cons1;

#if DEBUG==1
				if(1*4*504 <= debug_c && debug_c < 2*4*504){
					fprintf(stream2, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n", z, kloc, e, REAL(xgrid)[idx4(na - 1, e, kloc + 1, z)], coh, cons0, cons1, cons);
				}
				debug_c++;
#endif


				/* marginal utility x gross interest rate */
				REAL(rv1)[i] += beta*REAL(trans_ez)[
						idxTM(INTEGER(e_state)[i] - 1, INTEGER(z_state)[i] - 1, e, z)
					] * marg_util(cons, gamma)*(1.0 + rent);
			}
		}

	}

	for(i=0; i<ne*nk*nz; i++)
		if(NULL!=cs[i])
			CubicSplineDestroy(cs[i]);
	Free(cs);

#if DEBUG==1
	fclose(stream4);
	fclose(stream3);
	fclose(stream2);
	fclose(stream);
#endif

	UNPROTECT(1);
	return(rv1);
}

double marg_util(double cons, double risk_aversion)
{
	double value, first, second;

	if (cons > epsilon){
		if (risk_aversion == 1.0)
			return 1.0/cons;
		return pow(cons, -risk_aversion);
	}

	if (risk_aversion == 1.0) {
		value  = 1.0/epsilon;
		first  = -1.0 * (pow(epsilon,-2));
		second =  2.0 * (pow(epsilon,-3));
		return value + first*(cons-epsilon) + (second*pow(cons-epsilon,2)) / 2.0;
	}

	value  = pow(epsilon,-risk_aversion);
	first  = -risk_aversion * pow(epsilon,-risk_aversion-1.0);
	second = (-risk_aversion*(-risk_aversion-1.0)) * pow(epsilon,-risk_aversion-2.0);
	return value + first*(cons-epsilon) + pow(second*(cons-epsilon), 2) / 2.0;
}

double wage_f(double agg_cap, double agg_lab, double tfp, double alpha, double delta)
{
	return tfp*(1.0-alpha)*pow(agg_cap, alpha)*pow(agg_lab, -alpha);
}

double rent_f(double agg_cap, double agg_lab, double tfp, double alpha, double delta)
{
	return tfp*alpha*pow(agg_cap,alpha-1.0)*pow(agg_lab,1.0-alpha) - delta;
}

void	markov(double *tran, int ng, int ns, int *r)
{
	int	i, j, current;
	double	next_pr, rnd;

	GetRNGstate(); /* R requires this statement to start generating random numberes. */

	/* set initial state */
	current = r[0] = 1;

	/* main simulation part */
	for(j=1; j<ns; j++){
		next_pr = tran[current];
		for(i=0; i<ng; i++){
			if(runif(0, 1) <= next_pr){
				r[j] = i + 1;
				break;
			} else {
				next_pr += tran[current + ng*(i + 1)]; /* ng×ng行列 */ 
			}
		}
		current = r[j] - 1;
	}

	PutRNGstate(); /* R requires this statement to end generating random numberes. */
}

SEXP	Cmarkov(SEXP tran, SEXP ng, SEXP ns)
{
	SEXP	rv1, dim;
	int	i;

	if(!isMatrix(tran))
		error("A matrix is required for the 1st argument: 'tran'.");

	rv1 = PROTECT(allocVector(INTSXP, INTEGER(ns)[0]));

	markov(REAL(tran), INTEGER(ng)[0], INTEGER(ns)[0], INTEGER(rv1));

	UNPROTECT(1);
	return(rv1);
}


#define pidx(e, k, z) (e-1 + np[1]*(k) + np[1]*np[2]*(z-1))
#define tez_idx(e1, z1, e2, z2) ((e1) - 1 + nez[0]*((z1) - 1) + nez[0]*nez[1]*((e2) - 1) + nez[0]*nez[1]*nez[2]*((z2) - 1))
#define tzz_idx(z1, z2) ((z1) - 1 + nzz[0]*((z2) - 1))

#define IsSplineInterpolation 0

SEXP	law_of_motion_sim(SEXP policy, SEXP env)
{
	SEXP	k_path, u_path, z_path, dim, lst, lstnames;
	double	*cap0, *cap1, savings0, savings1, wght, rndw, count, sum, pr_next;
	int	*emp0, *emp1, i, j, index, it, numi, nums, kloc, glen, *np, *nez, *nzz, e, ne;
#if IsSplineInterpolation==1
	CubicSpline	**cs;
#endif
	SEXP	dparam, gparam, sparam, trans_prob, trans_ez, trans_zz, kgrid, grid;

/* Rのグローバル領域の変数を引っ張る */
	dparam = getvar("dparam", env);
	gparam = getvar("gparam", env);
	sparam = getvar("sparam", env);
	kgrid = getvar("kgrid", env);
	grid = getvar("grid", env);

	ne = toInt(getListElement("ne", dparam));
	trans_prob = getvar("trans_prob", env);
	trans_ez = getListElement("ez", trans_prob);
	trans_zz = getListElement("zz", trans_prob);

	numi = toInt(getListElement("numi", sparam));
	nums = toInt(getListElement("nums", sparam));

	glen = length(kgrid);

	k_path = PROTECT(allocVector(REALSXP, nums));
	u_path = PROTECT(allocVector(REALSXP, nums));
	z_path = PROTECT(allocVector(INTSXP, nums));

	/* 配列policyのサイズを取得しておく */
	dim = getAttrib(policy, R_DimSymbol);
	np = Calloc(length(dim), int);
	for(i=0; i<length(dim); i++){
		np[i] = INTEGER(dim)[i];
	}

	/* 配列trans_ezのサイズを取得しておく */
	dim = getAttrib(trans_ez, R_DimSymbol);
	nez = Calloc(length(dim), int);
	for(i=0; i<length(dim); i++){
		nez[i] = INTEGER(dim)[i];
	}

	/* 配列trans_zzのサイズを取得しておく */
	dim = getAttrib(trans_zz, R_DimSymbol);
	nzz = Calloc(length(dim), int);
	for(i=0; i<length(dim); i++){
		nzz[i] = INTEGER(dim)[i];
	}

	markov(REAL(trans_zz), nzz[0], nums, INTEGER(z_path));

	cap0 = Calloc(numi, double);
	sum = 0;
	for(i=0; i<numi; i++){
		cap0[i] = 35.0;
		sum += cap0[i];
	}

	emp0 = Calloc(numi, int);
	j = 0.9*numi;
	for(i=0;i<j;i++)
		emp0[i] = 1;
	for(i=j;i<numi;i++)
		emp0[i] = 2;

	REAL(k_path)[0] = sum/numi;
	REAL(u_path)[0] = 0.1;

	/* allocate cache for spline interpolation */
#if IsSplineInterpolation==1
	cs = Calloc(np[1]*np[2]*np[3], CubicSpline*);
	for(i=0; i<np[1]*np[2]*np[3]; i++)
		cs[i] = NULL;
#endif

	for(it=0; it<nums-1; it++){

		cap1 = Calloc(numi, double);
		emp1 = Calloc(numi, int);

		kloc = which(REAL(kgrid), glen, REAL(k_path)[it]);
		if(glen <= kloc + 1)
			kloc = glen - 2; /* kloc + 1 < glen so that [kloc + 1] can be used later. */

		wght = (REAL(k_path)[it] - REAL(kgrid)[kloc]) / (REAL(kgrid)[kloc + 1] - REAL(kgrid)[kloc]);
		if (wght > 1.0)
			wght = 1.0;
		else if(wght < 0.0)
			wght = 0.0;

		for(i=0, count=0, sum=0; i<numi; i++){

			GetRNGstate();
			rndw = runif(0, 1);
			PutRNGstate();

			if(1 <= emp0[i]){

				index = pidx(emp0[i], kloc, INTEGER(z_path)[it]);

#if IsSplineInterpolation==1
				if(NULL == cs[index])
					cs[index] = CubicSplineSetup(
						REAL(grid),
						&REAL(policy)[np[0] * index],
						np[0],
						BoundaryConditionNatural, 0, BoundaryConditionNatural, 0);

				CubicSplineInterpolation(cs[index], &cap0[i], &savings0, 1);
#else
				LinearInterpolation(REAL(grid), &REAL(policy)[np[0] * index], np[0], &cap0[i], &savings0, 1);
#endif

				index = pidx(emp0[i], kloc + 1, INTEGER(z_path)[it]);

#if IsSplineInterpolation==1
				if(NULL == cs[index])
					cs[index] = CubicSplineSetup(
						REAL(grid),
						&REAL(policy)[np[0] * index],
						np[0],
						BoundaryConditionNatural, 0, BoundaryConditionNatural, 0);

				CubicSplineInterpolation(cs[index], &cap0[i], &savings1, 1);

#else
				LinearInterpolation(REAL(grid), &REAL(policy)[np[0] * index], np[0], &cap0[i], &savings1, 1);
#endif

				cap1[i] = (1.0-wght)*savings0 + wght*savings1;

				pr_next = REAL(trans_ez)[ tez_idx(emp0[i], INTEGER(z_path)[it], 1, INTEGER(z_path)[it + 1]) ] / REAL(trans_zz)[ tzz_idx(INTEGER(z_path)[it], INTEGER(z_path)[it + 1]) ];

				for(e=1, emp1[i]=1; e<=ne; e++){
					if (rndw <= pr_next) {
						emp1[i] = e;
						break;
					} else {
						pr_next += REAL(trans_ez)[ tez_idx(emp0[i], INTEGER(z_path)[it], e + 1, INTEGER(z_path)[it + 1]) ] / REAL(trans_zz)[ tzz_idx(INTEGER(z_path)[it], INTEGER(z_path)[it + 1]) ];
					}
				}
			} else {
				fprintf(stderr, "emp0[%d]=%d\n", i, emp0[i]);
			}

			if(2 == emp1[i])
				count++;

			sum += cap1[i];
		}

		REAL(u_path)[it + 1] = count/numi;
		REAL(k_path)[it + 1] = sum/numi;

		Free(cap0);
		cap0 = cap1;
		Free(emp0);
		emp0 = emp1;
	}

	Free(cap0);
	Free(emp0);

/* deallocate cache for spline interpolation */
#if IsSplineInterpolation==1
	for(i=0; i<np[1]*np[2]*np[3]; i++)
		if(NULL!=cs[i])
			CubicSplineDestroy(cs[i]);
	Free(cs);
#endif

	Free(np);
	Free(nez);
	Free(nzz);

	lst = PROTECT(allocVector(VECSXP, 3)); /* 戻り値になるリストを作成 */

	/* リスト内に計算結果ベクトルなどを入れる */
	SET_VECTOR_ELT(lst, 0, k_path); 
	SET_VECTOR_ELT(lst, 1, u_path);
	SET_VECTOR_ELT(lst, 2, z_path);

	/* リストの要素の名前をリストの設定 */ 
	PROTECT(lstnames = allocVector(STRSXP, 3));
	SET_STRING_ELT(lstnames, 0,  mkChar("k"));
	SET_STRING_ELT(lstnames, 1,  mkChar("u"));
	SET_STRING_ELT(lstnames, 2,  mkChar("z"));

	setAttrib(lst, R_NamesSymbol, lstnames); 
	UNPROTECT(5);

	return(lst);
}

SEXP getListElement(const unsigned char *name, SEXP list)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	unsigned char	buf[1024];

	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
			elmt = VECTOR_ELT(list, i);
			return elmt;
		}
	
	sprintf(buf, "Couldn't find the variable: %s", name);
	error(buf);
	
	return R_NilValue; /* This part never runs. */
}

SEXP getvar(const unsigned char *name, SEXP env)
{
	SEXP		r;
	unsigned char	buf[1024];

	if(!isEnvironment(env))
		error("env should be an environment");
	
	if(R_NilValue == (r=findVar(install(name), env))){
		sprintf(buf, "Couldn't find the variable: %s", name);
		error(buf);
	}

	return r;
}

int toInt(SEXP a)
{
	if(isInteger(a)){
		return INTEGER(a)[0];
	}

	if(isReal(a)){
		return (int)(REAL(a)[0]);
	}

	error("Couldn't convert the variable to an integer.");
}

