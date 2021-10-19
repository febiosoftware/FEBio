/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "tens6d.h"

//-----------------------------------------------------------------------------
double calc_6ds_comp(double K[3][3], double Ri[3], double Rj[3], int i, int j, int k, int l, int m, int n)
{
	i -= 1; j -= 1;  k -= 1; l -= 1; m -= 1; n -= 1;
		
	return (1/72)*(   Ri[i]*Ri[j]*K[k][l]*Rj[m]*Rj[n] + Ri[i]*Ri[j]*K[k][n]*Rj[m]*Rj[l] + Ri[i]*Ri[j]*K[k][m]*Rj[l]*Rj[n] + Ri[i]*Ri[j]*K[k][l]*Rj[n]*Rj[m] + Ri[i]*Ri[j]*K[k][m]*Rj[n]*Rj[l] + Ri[i]*Ri[j]*K[k][n]*Rj[l]*Rj[m]     
	                + Ri[k]*Ri[j]*K[i][l]*Rj[m]*Rj[n] + Ri[k]*Ri[j]*K[i][n]*Rj[m]*Rj[l] + Ri[k]*Ri[j]*K[i][m]*Rj[l]*Rj[n] + Ri[k]*Ri[j]*K[i][l]*Rj[n]*Rj[m] + Ri[k]*Ri[j]*K[i][m]*Rj[n]*Rj[l] + Ri[k]*Ri[j]*K[i][n]*Rj[l]*Rj[m]
				    + Ri[j]*Ri[i]*K[k][l]*Rj[m]*Rj[n] + Ri[j]*Ri[i]*K[k][n]*Rj[m]*Rj[l] + Ri[j]*Ri[i]*K[k][m]*Rj[l]*Rj[n] + Ri[j]*Ri[i]*K[k][l]*Rj[n]*Rj[m] + Ri[j]*Ri[i]*K[k][m]*Rj[n]*Rj[l] + Ri[j]*Ri[i]*K[k][n]*Rj[l]*Rj[m]
					+ Ri[i]*Ri[k]*K[j][l]*Rj[m]*Rj[n] + Ri[i]*Ri[k]*K[j][n]*Rj[m]*Rj[l] + Ri[i]*Ri[k]*K[j][m]*Rj[l]*Rj[n] + Ri[i]*Ri[k]*K[j][l]*Rj[n]*Rj[m] + Ri[i]*Ri[k]*K[j][m]*Rj[n]*Rj[l] + Ri[i]*Ri[k]*K[j][n]*Rj[l]*Rj[m]
					+ Ri[j]*Ri[k]*K[i][l]*Rj[m]*Rj[n] + Ri[j]*Ri[k]*K[i][n]*Rj[m]*Rj[l] + Ri[j]*Ri[k]*K[i][m]*Rj[l]*Rj[n] + Ri[j]*Ri[k]*K[i][l]*Rj[n]*Rj[m] + Ri[j]*Ri[k]*K[i][m]*Rj[n]*Rj[l] + Ri[j]*Ri[k]*K[i][n]*Rj[l]*Rj[m]
					+ Ri[k]*Ri[i]*K[j][l]*Rj[m]*Rj[n] + Ri[k]*Ri[i]*K[j][n]*Rj[m]*Rj[l] + Ri[k]*Ri[i]*K[j][m]*Rj[l]*Rj[n] + Ri[k]*Ri[i]*K[j][l]*Rj[n]*Rj[m] + Ri[k]*Ri[i]*K[j][m]*Rj[n]*Rj[l] + Ri[k]*Ri[i]*K[j][n]*Rj[l]*Rj[m]
					+ Ri[l]*Ri[m]*K[n][i]*Rj[j]*Rj[k] + Ri[l]*Ri[m]*K[n][k]*Rj[j]*Rj[i] + Ri[l]*Ri[m]*K[n][j]*Rj[i]*Rj[k] + Ri[l]*Ri[m]*K[n][i]*Rj[k]*Rj[j] + Ri[l]*Ri[m]*K[n][j]*Rj[k]*Rj[i] + Ri[l]*Ri[m]*K[n][k]*Rj[i]*Rj[j]
					+ Ri[n]*Ri[m]*K[l][i]*Rj[j]*Rj[k] + Ri[n]*Ri[m]*K[l][k]*Rj[j]*Rj[i] + Ri[n]*Ri[m]*K[l][j]*Rj[i]*Rj[k] + Ri[n]*Ri[m]*K[l][i]*Rj[k]*Rj[j] + Ri[n]*Ri[m]*K[l][j]*Rj[k]*Rj[i] + Ri[n]*Ri[m]*K[l][k]*Rj[i]*Rj[j]
					+ Ri[m]*Ri[l]*K[n][i]*Rj[j]*Rj[k] + Ri[m]*Ri[l]*K[n][k]*Rj[j]*Rj[i] + Ri[m]*Ri[l]*K[n][j]*Rj[i]*Rj[k] + Ri[m]*Ri[l]*K[n][i]*Rj[k]*Rj[j] + Ri[m]*Ri[l]*K[n][j]*Rj[k]*Rj[i] + Ri[m]*Ri[l]*K[n][k]*Rj[i]*Rj[j]
					+ Ri[l]*Ri[n]*K[m][i]*Rj[j]*Rj[k] + Ri[l]*Ri[n]*K[m][k]*Rj[j]*Rj[i] + Ri[l]*Ri[n]*K[m][j]*Rj[i]*Rj[k] + Ri[l]*Ri[n]*K[m][i]*Rj[k]*Rj[j] + Ri[l]*Ri[n]*K[m][j]*Rj[k]*Rj[i] + Ri[l]*Ri[n]*K[m][k]*Rj[i]*Rj[j]
					+ Ri[m]*Ri[n]*K[l][i]*Rj[j]*Rj[k] + Ri[m]*Ri[n]*K[l][k]*Rj[j]*Rj[i] + Ri[m]*Ri[n]*K[l][j]*Rj[i]*Rj[k] + Ri[m]*Ri[n]*K[l][i]*Rj[k]*Rj[j] + Ri[m]*Ri[n]*K[l][j]*Rj[k]*Rj[i] + Ri[m]*Ri[n]*K[l][k]*Rj[i]*Rj[j]
					+ Ri[n]*Ri[l]*K[m][i]*Rj[j]*Rj[k] + Ri[n]*Ri[l]*K[m][k]*Rj[j]*Rj[i] + Ri[n]*Ri[l]*K[m][j]*Rj[i]*Rj[k] + Ri[n]*Ri[l]*K[m][i]*Rj[k]*Rj[j] + Ri[n]*Ri[l]*K[m][j]*Rj[k]*Rj[i] + Ri[n]*Ri[l]*K[m][k]*Rj[i]*Rj[j]);
}

//-----------------------------------------------------------------------------
void calculate_e2O(tens6ds& e, double K[3][3], double Ri[3], double Rj[3] )
{
	e.d[ 0] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1, 1);
	e.d[ 1] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1, 2);
	e.d[ 2] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1, 3);
	e.d[ 3] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 2, 2);
	e.d[ 4] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 2, 3);
	e.d[ 5] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 3, 3);
	e.d[ 6] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 2, 2, 2);
	e.d[ 7] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 2, 2, 3);
	e.d[ 8] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 2, 3, 3);
	e.d [9] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 3, 3, 3);
	
	e.d[10] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 1, 2, 2);
	e.d[11] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 1, 2, 3);
	e.d[12] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 1, 3, 3);
	e.d[13] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 2, 2, 2);
	e.d[14] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 2, 2, 3);
	e.d[15] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 2, 3, 3);
	e.d[16] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 3, 3, 3);

	e.d[17] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 1, 2, 2);
	e.d[18] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 1, 2, 3);
	e.d[19] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 1, 3, 3);
	e.d[20] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 2, 2, 2);
	e.d[21] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 2, 2, 3);
	e.d[22] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 2, 3, 3);
	e.d[23] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 3, 3, 3);

	e.d[24] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 1, 3, 3);
	e.d[25] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 2, 2, 2);
	e.d[26] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 2, 2, 3);
	e.d[27] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 2, 3, 3);
	e.d[28] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 3, 3, 3);

	e.d[29] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 1, 3, 3);
	e.d[30] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 2, 2, 2);
	e.d[31] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 2, 2, 3);
	e.d[32] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 2, 3, 3);
	e.d[33] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 3, 3, 3);

	e.d[34] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 2, 2, 2);
	e.d[35] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 2, 2, 3);
	e.d[36] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 2, 3, 3);
	e.d[37] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 3, 3, 3);

	e.d[38] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 2, 2, 2);
	e.d[39] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 2, 2, 3);
	e.d[40] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 2, 3, 3);
	e.d[41] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 3, 3, 3);

	e.d[42] += calc_6ds_comp(K, Ri, Rj, 2, 2, 3, 2, 3, 3);
	e.d[43] += calc_6ds_comp(K, Ri, Rj, 2, 2, 3, 3, 3, 3);

	e.d[44] += calc_6ds_comp(K, Ri, Rj, 2, 3, 3, 3, 3, 3);

	e.d[45] += calc_6ds_comp(K, Ri, Rj, 3, 3, 3, 3, 3, 3);
}

