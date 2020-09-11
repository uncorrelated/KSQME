#include <stdio.h>
#include <stdlib.h>
#include "which.h"

/*
	a ∈ [array[0], array[size-1]]でarray[m] ≤ a < array[m+1]となるmを返す
	a < array[0]のときは0、a > array[size-1]のときはsize-1
*/
int	which(double *array, int size, double a)
{
	int	l, r, m;
	
	l = 0;
	r = size - 1;
	/*
		if(a >= array[size-1]) return size-1; をループ前につけて、
		while(l + 1 < r){ にして、if(a < array[m + 1]){ ... }を無くして、
		l = m + 1;をl = m;にした方がたぶん速い。
	*/
	while(l < r){
		m = (l + r)/2;

		if(a >= array[m]){
			if(a < array[m + 1]){ /* 端数切り捨てによりmは必ずm<size-1 */
				return m;
			}
			l = m + 1;/* int型変数を/2をすると切り捨てになるので、下端を動かすときは+1しないと収束しない */
		} else {
			r = m;
		}
	}
	return (l + r)/2;
}

