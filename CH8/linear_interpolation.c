#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linear_interpolation.h"

#define OutRangeProcessingType 1

void LinearInterpolation(double *knot_x, double *knot_y, unsigned kn, double *x, double *y, unsigned n){
	int	i, l, r, m, mp1;
	double	p, dx_l, dx_r, dx, dy;

	for(i=0; i<n; i++){

		if(x[i] <= knot_x[0]){
			#if OutRangeProcessingType==1
				y[i] = knot_y[0] + (knot_y[1] - knot_y[0])/(knot_x[1] - knot_x[0])*(x[i] - knot_x[0]);
			#else
				y[i] = knot_y[0];
			#endif
			continue;
		}



		if(x[i] >= knot_x[kn - 1]){
			#if OutRangeProcessingType==1
				y[i] = knot_y[kn - 1] + (knot_y[kn - 1] - knot_y[kn - 2])/(knot_x[kn - 1] - knot_x[kn - 2])*(x[i] - knot_x[kn - 1]);
			#else
				y[i] = knot_y[kn - 1];
			#endif
			continue;
		}

		/* binary search */
		l = 0;
		r = kn - 1;

		while(l < r){
			m = (l + r)/2;

/*			fprintf(stderr, "x[%d]=%f, l:%d, r:%d\n", i, x[i], l, r); */

			if(knot_x[m] <= x[i]){
				mp1 = m + 1; /* m is always less than kn - 1. */
				if(x[i] < knot_x[mp1]){
					break;
				}
				l = m + 1;
			} else {
				r = m;
			}
		}

		mp1 = m + 1;
		dx_l = fabs(knot_x[m] - x[i]);
		dx_r = fabs(knot_x[mp1] - x[i]);
		p = dx_r/(dx_l + dx_r); /* dx_l + dx_r > 0 because x[i] == knot_x[n - 1] doesn't appear. */

		y[i] = p*knot_y[m] + (1-p)*knot_y[mp1];
/*		fprintf(stderr, "x[%d]=%f, y[%d]=%f, l:%d, m:%d, %d, p:%f, r:%d\n", i, x[i], i, y[i], l, m, mp1, p, r);	*/
	}
}

