#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cubic_spline.h"

#define TooSmallNumber 1.0e-8

void CubicSplineInterpolation(CubicSpline *s3, double *x, double *y, unsigned n){
	int	i, l, r, m;
	double	a, b, c, d, dx;

	for(i=0; i<n; i++){

		/* binary search */
		l = 0;
		r = s3->n - 1;

		while(l < r){
			m = (l + r)/2;

			if(s3->x[m] <= x[i]){
				if(x[i] < s3->x[m + 1]){ /* m is always less than kn - 1. */
					break;
				}
				l = m + 1;
			} else {
				r = m;
			}
		}

		a = (2.0*(s3->y[m] - s3->y[m+1])/s3->H[m] + s3->rh[m] + s3->rh[m+1]) / (s3->H[m] * s3->H[m]);
		b = (3.0*(s3->y[m+1] - s3->y[m])/s3->H[m] - 2.0*s3->rh[m] - s3->rh[m+1]) / s3->H[m];
		c = s3->rh[m];
		d = s3->y[m]; 
		dx = x[i] - s3->x[m];
		y[i] = dx*(dx*(dx*a + b) + c) + d;
	}
}

void CubicSplineInterpolationSortedInput(CubicSpline *s3, double *x, double *y, unsigned n){
	unsigned	i, j;
	double		a, b, c, d, dx;

	for(j=0, i=0; i < s3->n-1; i++){
		a = (2.0*(s3->y[i] - s3->y[i+1])/s3->H[i] + s3->rh[i] + s3->rh[i+1]) / (s3->H[i] * s3->H[i]);
		b = (3.0*(s3->y[i+1] - s3->y[i])/s3->H[i] - 2.0*s3->rh[i] - s3->rh[i+1]) / s3->H[i];
		c = s3->rh[i];
		d = s3->y[i]; 
/*
	fprintf(stderr, "1. a=%f, b=%f, c=%f, d=%f\n", a, b, c, d);
*/
		while( x[j] < s3->x[i+1] ){
			dx = x[j] - s3->x[i];
			y[j] = dx*(dx*(dx*a + b) + c) + d;

/*
	fprintf(stderr, "2. dx=%f, x[%d]=%f, y[%d]=%f\n", dx, j, x[j], j, y[j]);
*/
			if(n <= ++j){
				return;
			}
		}
	}

	while(j < n){
		y[j++] = s3->y[s3->n - 1];
	}
}

void	CubicSplineDestroy(CubicSpline *s3)
{
	free(s3->H);
	free(s3->rh);
	free(s3);
}

void	sortCoordinate(double x[], double y[], int l, int r)
{
	double tx,ty;
	int s,m,i,j;

	if(l>=r)
		return;

	m = l + rand() % (r - l) + 1;
	s = l;

	tx = x[s]; ty = y[s];
	x[s] = x[m]; y[s] = y[m];
	x[m] = tx; y[m] = ty;

/*
	fprintf(stderr, "recursive l:%d, r:%d\n", l, r);
	for(i=l;i<=r;i++){
		fprintf(stderr, "x[%d]:%f, y[%d]:%f\n", i, x[i], i, y[i]);
	}	
*/

	i = l + 1;
	j = r;
	for(;;){
/*
		fprintf(stderr, "(1)i:%d, j:%d\n", i, j);
*/
		while(x[s] > x[i]){
			if(i >= j){
				goto ESCAPE;
			}
			i++;
		}
		while(x[s] <= x[j]){
			if(i >= j){
				goto ESCAPE;
			}
			j--;
		}
/*
		fprintf(stderr, "(2)i:%d, j:%d\n", i, j);
		fprintf(stderr, "xchg i:%d, j:%d\n", i, j);
*/
		tx = x[i]; ty = y[i];
		x[i] = x[j]; y[i] = y[j];
		x[j] = tx; y[j] = ty;
	}

	ESCAPE:
	if(x[s] > x[j]){
		tx = x[s]; ty = y[s];
		x[s] = x[j]; y[s] = y[j];
		x[j] = tx; y[j] = ty;
	}
	sortCoordinate(x, y, l, j - 1);
	sortCoordinate(x, y, j, r);
}

CubicSpline *CubicSplineSetup(double x[], double y[], int n, int boundary_condition_0, double y_dash_0, int boundary_condition_n, double y_dash_n)
{
	double *a, *b, *c, *dx, *rh;
	int i, j, scode;
	CubicSpline *s3;

	if(3>n)
		return NULL;

	s3 = malloc(sizeof(CubicSpline));
	s3->rh = rh = calloc(n, sizeof(double));
	s3->H = dx = calloc(n, sizeof(double));  

	a = calloc(n, sizeof(double)); /* left series of the tridiagonal matrix */
	b = calloc(n, sizeof(double)); /* center series of the tridiagonal matrix */
	c = calloc(n, sizeof(double)); /* right series of the tridiagonal matrix */

	/* Following lecture notes in the Internet to make a tridiagonal matrix for cubic spline. */
	for(i=0; i<n-1; i++)
		dx[i] = x[i+1] - x[i];

	a[0] = 0.0;

	if(BoundaryConditionNatural == boundary_condition_0){
		b[0] = 2.0*dx[0]; 
		c[0] = dx[0]; 
		rh[0] = 3.0*(y[1] - y[0]);
	} else if(BoundaryConditionFixed == boundary_condition_0){
		b[0] = 1.0;    
		c[0] = 0.0;  
		rh[0] = y_dash_0;
	}

	for(i=1; i<n-1; i++)
	{
		a[i] = dx[i];   
		b[i] = 2.0*(dx[i-1]+dx[i]);   
		c[i] = dx[i-1];
		rh[i] = 3.0*((y[i]-y[i-1])*dx[i]/dx[i-1]+(y[i+1]-y[i])*dx[i-1]/dx[i]);
	}

	c[n-1] = 0.0;
	if(BoundaryConditionNatural == boundary_condition_n){ 
		a[n-1] = dx[n-2];   
		b[n-1] = 2.0*dx[n-2];
		rh[n-1] = 3.0*(y[n-1]-y[n-2]);
	} else if(BoundaryConditionFixed == boundary_condition_n){
		a[n-1] = 0.0; 
		b[n-1] = 1.0; 
		rh[n-1] = y_dash_n;
	}

/*
for(i=0; i<n; i++){
	fprintf(stderr, "a[%d]=%f, b[%d]=%f, c[%d]=%f, dx[%d]=%f, rh[%d]=%f\n", i, a[i], i, b[i], i, c[i], i, dx[i], i, rh[i]);
}
*/

	/* LU-decomposing the tridiagonal matrix */
	for(i=0; i<n; i++){
		if(TooSmallNumber > fabs(b[i])){

			free(c);
			free(b);
			free(a);

			CubicSplineDestroy(s3);

			return NULL;
		}
	}

	c[0] /= b[0];
	for(i=1; i<n-1; i++){
	    b[i] -= c[i-1]*a[i];      
	    c[i] /= b[i];
	}
	b[n-1] -= c[n-2]*a[n-1];

	/* Forward substitution */
	rh[0] /= b[0];  
	for(i=1; i<n; i++)
		rh[i] = (rh[i]-a[i]*rh[i-1])/b[i]; 
	/* Backward substitution */
	for(i=n-2; i>=0; i--)
		rh[i] -= rh[i+1]*c[i]; 

	s3->x = x;
	s3->y = y;
	s3->n = n;

	free(c);
	free(b);
	free(a);

	return s3;
}

