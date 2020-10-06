#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

double pchipend(double h1, double h2, double del1, double del2);
double	pchip(const double *knot_x, const double *knot_y, const int kn, const double x);

/* Piecwise Cubic Hermitean Interpolation Polynomials: Matlab Compatible one */

#define DEBUG 0

double	pchip(NumericVector knot_x, NumericVector knot_y, double x) {
	return pchip((as< std::vector<double> >(knot_x)).data(), (as< std::vector<double> >(knot_y)).data(), knot_x.length(), x);
}

double	pchip(const double *knot_x, const double *knot_y, const int kn, const double x) {
	int	i, j, l, r, m, mp1;
	double	s;
	double	h[3], delta[3], a, b, d[3];
	int	t;
	double	w1, w2;

	/* binary search */
	l = 0;
	r = kn - 1;

	while(l < r){
		m = (l + r)/2;
		if(knot_x[m] <= x){
			mp1 = m + 1;
			if(x < knot_x[mp1]){
				break;
			}
			l = m + 1;
		} else {
			r = m;
		}
	}


// hとdeltaはm-1〜m+1の範囲を計算; h[1],delta[1]が基点
// dはm〜m+1の範囲計算; d[1]が基点
// aとbはmのところを計算

	/* First derivatives */
	t = m;
	for(i=t<1?1:0; i<=2 && t+i<kn; i++){
		h[i] = knot_x[t + i] - knot_x[t + i - 1];
		delta[i] = (knot_y[t + i] - knot_y[t + i - 1]) / h[i];
#if DEBUG==1
		fprintf(stderr, "h[%d]\t%f\t", t + i, h[i]);
		fprintf(stderr, "delta[%d]\t%f\n", t + i, delta[i]);
#endif
	}

	for(i=0;i<=1;i++){
		/* Slopes at endpoints */
		if(0>=t+i){
			d[i] = pchipend(h[0], h[1], delta[0], delta[1]);
			continue;
		}
		if(kn-2<=t+i){
			d[i] = pchipend(h[2], h[1], delta[2], delta[1]);
			continue;
		}
		/* Slopes at interior points */
		d[i] = 0;
		if(0 < delta[i]*delta[i + 1]){
			w1 = 2*h[i + 1] + h[i];
			w2 = h[i + 1] + 2*h[i];
			d[i] = (w1 + w2) / (w1/delta[i] + w2/delta[i + 1]);
		}
	}

#if DEBUG==1
	for(i=0;i<2;i++){
		fprintf(stderr, "d[%d]\t%f\t\n", t + i, d[i]);
	}
#endif

	/* Piecewise polynomial coefficients */
	a = (3*delta[1] - 2*d[0] - d[1]) / h[1];
	b = (d[0] - 2*delta[1] + d[1]) / (h[1]*h[1]);

#if DEBUG==1
	fprintf(stderr, "m:%d\tt:%d\ta:%f\tb:%f\td:%f\n", m, t, a, b, d[0]);
#endif

	/* Evaluate interpolant */
	s = x - knot_x[m];
	return knot_y[m] + s*(d[0] + s*(a + s*b));
}


double pchipend(double h1, double h2, double del1, double del2){
	double d;

	// Noncentered, shape-preserving, three-point formula.
	d = ((2*h1 + h2)*del1 - h1*del2) / (h1 + h2);
	if(d*del1 < 0){
		d = 0;
	}else if((del1*del2< 0) && (fabs(d) > fabs(3*del1))) {
		d = 3*del1;
	}
	return d;
}


// [[Rcpp::export]]
NumericVector pchip(NumericVector knot_x, NumericVector knot_y, NumericVector x) {
	int	i, n = x.length();
	NumericVector y(n);

	if(knot_x.length() != knot_y.length()){
		stop("the lengths of knots_x and knots_y are different.");
	}

	if(knot_x.length() < 4){
		stop("the length of knots is too short.");
	}

	for(i=0; i<n; i++)
		y[i] = pchip(knot_x, knot_y, x[i]);

	return y;
}


