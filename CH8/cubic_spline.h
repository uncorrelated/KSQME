#ifndef _CubicSpline
#define _CubicSpline 1

typedef struct {
	unsigned	n; /* the number of knots */ 
	double		*x, *y; /* cordinates of knots */
	double		*rh; /* right hand of the tridiagonal matrix ---> the result of LU decomposition */
	double		*H;
} CubicSpline;

#define BoundaryConditionNatural 0
#define BoundaryConditionFixed 1

void sortCoordinate(double x[], double y[], int l, int r);
CubicSpline *CubicSplineSetup(double x[], double y[], int n, int boundary_condition_0, double y_dash_0, int boundary_condition_n, double y_dash_n);
void CubicSplineInterpolation(CubicSpline *s3, double *x, double *y, unsigned n);
void CubicSplineInterpolationSortedInput(CubicSpline *s3, double *x, double *y, unsigned n);
void CubicSplineDestroy(CubicSpline *s3);

#endif
