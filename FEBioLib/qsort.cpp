#include "stdafx.h"

#define SWAP(a,b) { itemp = a; a = b; b = itemp; }

void qsort(int n, int* arr, int* indx)
{
	const int M = 7;
	const int NSTACK = 50;
	int i, indxt, ir=n-1, itemp, j, k, l=0;
	int jstack = -1, istack[NSTACK];
	int a;

	for (j=0; j<n; ++j) indx[j] = j;
	for (;;)
	{
		if (ir-l<M)
		{
			for (j=l+1; j<=ir; ++j)
			{
				indxt = indx[j];
				a = arr[indxt];
				for (i=j-1; i>=l; --i)
				{
					if (arr[indx[i]] <= a) break;
					indx[i+1] = indx[i];
				}
				indx[i+1] = indxt;
			}
			if (jstack == -1) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		}
		else
		{
			k = ((l+ir+2) >> 1) - 1;
			SWAP(indx[k], indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) SWAP(indx[l], indx[ir]);
			if (arr[indx[l+1]] > arr[indx[ir]]) SWAP(indx[l+1], indx[ir]);
			if (arr[indx[l]] > arr[indx[l+1]]) SWAP(indx[l], indx[l+1]);
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a = arr[indxt];
			for(;;)
			{
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j<i) break;
				SWAP(indx[i], indx[j]);
			}
			indx[l+1] = indx[j];
			indx[j] = indxt;
			jstack += 2;
//			if (jstack < NSTACK) throw "a fit";
			if (ir-i+1 >= j-l)
			{
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j-1;
			}
			else
			{
				istack[jstack] = j-1;
				istack[jstack-1] = l;
				l=i;
			}
		}
	}
}

