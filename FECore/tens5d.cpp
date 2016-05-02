#include "stdafx.h"
#include "tens5d.h"
#include <assert.h>

//-----------------------------------------------------------------------------
double calc_5ds_comp(double K[3][3], double Ri[3], double Rj[3], int i, int j, int k, int l, int m)
{
	i -= 1; j -= 1;  k -= 1; l -= 1; m -= 1; 
		
	return (1/2)*((1/24)*(   Ri[i]*K[j][k]*Rj[l]*Rj[m] + Ri[i]*K[j][m]*Rj[l]*Rj[k] + Ri[i]*K[j][l]*Rj[k]*Rj[m] + Ri[i]*K[j][k]*Rj[m]*Rj[l] + Ri[i]*K[j][l]*Rj[m]*Rj[k] + Ri[i]*K[j][m]*Rj[k]*Rj[l]
					       + Ri[j]*K[i][k]*Rj[l]*Rj[m] + Ri[j]*K[i][m]*Rj[l]*Rj[k] + Ri[j]*K[i][l]*Rj[k]*Rj[m] + Ri[j]*K[i][k]*Rj[m]*Rj[l] + Ri[j]*K[i][l]*Rj[m]*Rj[k] + Ri[j]*K[i][m]*Rj[k]*Rj[l]
					       + Ri[l]*K[m][i]*Rj[j]*Rj[k] + Ri[l]*K[m][k]*Rj[j]*Rj[i] + Ri[l]*K[m][j]*Rj[i]*Rj[k] + Ri[l]*K[m][i]*Rj[k]*Rj[j] + Ri[l]*K[m][j]*Rj[k]*Rj[i] + Ri[l]*K[m][k]*Rj[i]*Rj[j]
					       + Ri[m]*K[l][i]*Rj[j]*Rj[k] + Ri[m]*K[l][k]*Rj[j]*Rj[i] + Ri[m]*K[l][j]*Rj[i]*Rj[k] + Ri[m]*K[l][i]*Rj[k]*Rj[j] + Ri[m]*K[l][j]*Rj[k]*Rj[i] + Ri[m]*K[l][k]*Rj[i]*Rj[j])
		        + (1/24)*(   Ri[i]*Ri[j]*K[k][l]*Rj[m] + Ri[i]*Ri[j]*K[m][l]*Rj[k] + Ri[i]*Ri[j]*K[l][k]*Rj[m] + Ri[i]*Ri[j]*K[k][m]*Rj[l] + Ri[i]*Ri[j]*K[l][m]*Rj[k] + Ri[i]*Ri[j]*K[m][k]*Rj[l]
					       + Ri[j]*Ri[i]*K[k][l]*Rj[m] + Ri[j]*Ri[i]*K[m][l]*Rj[k] + Ri[j]*Ri[i]*K[l][k]*Rj[m] + Ri[j]*Ri[i]*K[k][m]*Rj[l] + Ri[j]*Ri[i]*K[l][m]*Rj[k] + Ri[j]*Ri[i]*K[m][k]*Rj[l]
					       + Ri[l]*Ri[m]*K[i][j]*Rj[k] + Ri[l]*Ri[m]*K[k][j]*Rj[i] + Ri[l]*Ri[m]*K[j][i]*Rj[k] + Ri[l]*Ri[m]*K[i][k]*Rj[j] + Ri[l]*Ri[m]*K[j][k]*Rj[i] + Ri[l]*Ri[m]*K[k][i]*Rj[j]
					       + Ri[m]*Ri[l]*K[i][j]*Rj[k] + Ri[m]*Ri[l]*K[k][j]*Rj[i] + Ri[m]*Ri[l]*K[j][i]*Rj[k] + Ri[m]*Ri[l]*K[i][k]*Rj[j] + Ri[m]*Ri[l]*K[j][k]*Rj[i] + Ri[m]*Ri[l]*K[k][i]*Rj[j]));
}

//-----------------------------------------------------------------------------
void calculate_d2O(tens5ds& d, double K[3][3], double Ri[3], double Rj[3] )
{
	d.d[ 0] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1);
	d.d[ 1] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 1, 2);
	d.d[ 2] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 1, 3);
	d.d[ 3] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 2, 2);
	d.d[ 4] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 2, 3);
	d.d[ 5] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 3, 3);
	d.d[ 6] += calc_5ds_comp(K, Ri, Rj, 1, 1, 2, 2, 2);
	d.d[ 7] += calc_5ds_comp(K, Ri, Rj, 1, 1, 2, 2, 3);
	d.d[ 8] += calc_5ds_comp(K, Ri, Rj, 1, 1, 2, 3, 3);
	d.d[ 9] += calc_5ds_comp(K, Ri, Rj, 1, 1, 3, 3, 3);
	
	d.d[10] += calc_5ds_comp(K, Ri, Rj, 1, 2, 1, 2, 2);
	d.d[11] += calc_5ds_comp(K, Ri, Rj, 1, 2, 1, 2, 3);
	d.d[12] += calc_5ds_comp(K, Ri, Rj, 1, 2, 1, 3, 3);
	d.d[13] += calc_5ds_comp(K, Ri, Rj, 1, 2, 2, 2, 2);
	d.d[14] += calc_5ds_comp(K, Ri, Rj, 1, 2, 2, 2, 3);
	d.d[15] += calc_5ds_comp(K, Ri, Rj, 1, 2, 2, 3, 3);
	d.d[16] += calc_5ds_comp(K, Ri, Rj, 1, 2, 3, 3, 3);

	d.d[17] += calc_5ds_comp(K, Ri, Rj, 1, 3, 1, 3, 3);
	d.d[18] += calc_5ds_comp(K, Ri, Rj, 1, 3, 2, 2, 2);
	d.d[19] += calc_5ds_comp(K, Ri, Rj, 1, 3, 2, 2, 3);
	d.d[20] += calc_5ds_comp(K, Ri, Rj, 1, 3, 2, 3, 3);
	d.d[21] += calc_5ds_comp(K, Ri, Rj, 1, 3, 3, 3, 3);

	d.d[22] += calc_5ds_comp(K, Ri, Rj, 2, 2, 2, 2, 2);
	d.d[23] += calc_5ds_comp(K, Ri, Rj, 2, 2, 2, 2, 3);
	d.d[24] += calc_5ds_comp(K, Ri, Rj, 2, 2, 2, 3, 3);
	d.d[25] += calc_5ds_comp(K, Ri, Rj, 2, 2, 3, 3, 3);

	d.d[26] += calc_5ds_comp(K, Ri, Rj, 2, 3, 2, 3, 3);

	d.d[27] += calc_5ds_comp(K, Ri, Rj, 2, 3, 3, 3, 3);

	d.d[28] += calc_5ds_comp(K, Ri, Rj, 3, 3, 3, 3, 3);
}
