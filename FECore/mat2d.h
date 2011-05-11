#pragma once

class mat2d
{
public:
	mat2d(){}
	mat2d(double a00, double a01, double a10, double a11)
	{
		d[0][0] = a00; d[0][1] = a01;
		d[1][0] = a10; d[1][1] = a11;
	}

	double& operator () (int i, int j) { return d[i][j]; }
	double operator () (int i, int j) const { return d[i][j]; }

	double* operator [] (int i) { return d[i]; }

	mat2d inverse()
	{
		double Di = 1/(d[0][0]*d[1][1] - d[0][1]*d[1][0]);
		return mat2d(d[1][1]*Di, -d[0][1]*Di, -d[1][0]*Di, d[0][0]*Di);
	}
	
protected:
	double	d[2][2];
};
