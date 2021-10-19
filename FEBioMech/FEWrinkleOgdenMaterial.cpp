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
#include "FEWrinkleOgdenMaterial.h"
#include "FECore/mat2d.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEWrinkleOgdenMaterial, FEMembraneMaterial)
	ADD_PARAMETER(m_u, FE_RANGE_GREATER(0.0), "mu");
	ADD_PARAMETER(m_a, FE_RANGE_GREATER_OR_EQUAL(2.0), "alpha");
	ADD_PARAMETER(m_bwrinkle, "wrinkle");
	ADD_PARAMETER(m_l0, 2, "prestretch");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
void eig(double *A, double *B, double *C);
void eigen2d(double l[2], double r[4], double m[4]);

//-----------------------------------------------------------------------------
FEWrinkleOgdenMaterial::FEWrinkleOgdenMaterial(FEModel* pfem) : FEMembraneMaterial(pfem)
{
	m_a = 0;
	m_u = 0;
	m_l0[0] = m_l0[1] = 0.0;
	m_bwrinkle = true;
}

//-----------------------------------------------------------------------------
void FEWrinkleOgdenMaterial::principals(FEMaterialPoint& mp, double l[2], double v[4])
{
	FEMembraneMaterialPoint& pt = *mp.ExtractData<FEMembraneMaterialPoint>();

	// get the def gradient
	double* g = pt.g;

	// extract the 2D component
	mat2d Fm(1 + g[0], g[3], g[1], 1 + g[4]);

	// 2D right Cauchy-Green tensor
	mat2d Cm = Fm.transpose()*Fm;

	// get eigenvalues and vectors
	eig(l, v, Cm[0]);

	// get principal stretches 
	l[0] = sqrt(l[0]) + m_l0[0];
	l[1] = sqrt(l[1]) + m_l0[1];

	// lambda0 should be the bigger stretch
	if (l[1] > l[0])
	{
		double t = l[0]; l[0] = l[1]; l[1] = t; 
		t = v[0]; v[0] = v[2]; v[2] = t;
		t = v[1]; v[1] = v[3]; v[3] = t;
	}
}

//-----------------------------------------------------------------------------
void FEWrinkleOgdenMaterial::Stress(FEMaterialPoint& mp, double s[3])
{
	// get the principal strains and vectors
	double l[2], ec[4];
	principals(mp, l, ec);

	// evaluate stress based on principal stretches
	// note that we calculate the 2PK stress in principal directions
	double S[2][2] = {0};
	if (m_bwrinkle && (l[0] < 0.999)) // biaxial compression (wrinkling in both directions)
	{
		// stress is zero
	}
	else if (m_bwrinkle && (l[1] < 1.0/sqrt(l[0]))) // uniaxial tension (wrinkling in one direction)
	{
		S[0][0] = m_u*(pow(l[0], m_a-2.0)-pow(l[0],-2.0-0.5*m_a));
	}
	else // biaxial tension (no wrinkling)
	{
		S[0][0] = m_u*(-pow(l[0], -m_a-2.0)*pow(l[1], -m_a) + pow(l[0], m_a-2.0));
		S[1][1] = m_u*(-pow(l[1], -m_a-2.0)*pow(l[0], -m_a) + pow(l[1], m_a-2.0));
	}

	// setup coordinate transformation matrix to
	// transform from principal coordinates to element coordinates
	double Q[2][2] = {
		{ ec[0],  ec[1]},
		{-ec[1],  ec[0]}};

	// now we convert back to element matrix
	double Sv[2][2];
	for (int i=0; i<2; ++i)
		for (int j=0; j<2; ++j)
		{
			Sv[i][j] = 0;
			for (int k=0; k<2; ++k) 
				for (int l=0; l<2; ++l) Sv[i][j] += Q[k][i]*S[k][l]*Q[l][j];
		}

	// copy stress components
	s[0] = Sv[0][0]; s[1] = Sv[1][1]; s[2] = Sv[0][1];
}

//-----------------------------------------------------------------------------
void FEWrinkleOgdenMaterial::Tangent(FEMaterialPoint &mp, double D[3][3])
{
	// get the principal strains and vectors
	double l[2], ec[4];
	principals(mp, l, ec);

	// evaluate stress based on principal stretches
	// note that we calculate the 2PK stress in principal directions
	double S[2] = {0};
	if (m_bwrinkle && (l[0] < 0.999)) // biaxial compression (wrinkling in both directions)
	{
		// no stress
	}
	else if (m_bwrinkle && (l[1] < 1.0/sqrt(l[0]))) // uniaxial tension (wrinkling in one direction)
	{
		S[0] = m_u*(pow(l[0], m_a-2.0)-pow(l[0],-2.0-0.5*m_a));
	}
	else // biaxial tension (no wrinkling)
	{
		S[0] = m_u*(-pow(l[0], -m_a-2.0)*pow(l[1], -m_a) + pow(l[0], m_a-2.0));
		S[1] = m_u*(-pow(l[1], -m_a-2.0)*pow(l[0], -m_a) + pow(l[1], m_a-2.0));
	}

	// evaluate tangent based on principal stretches
	double Cp[2][2] = {0}, ks = 0;
/*	if (m_bwrinkle && (l[0] < 1)) // biaxial compression (wrinkling in both directions)
	{
		// tangent is zero
	}
	else if (m_bwrinkle && (l[1] < sqrt(l[0])))  // uniaxial tension (wrinkling in uniaxial)
	{
		Cp[0][0] = m_u*((2+m_a*0.5)*pow(l[0], -3-m_a*0.5) + (m_a-2.0)*pow(l[0], m_a-3.0))/l[0];
	}
	else // biaxial tension (no wrinkling)
*/	{
		Cp[0][0] = m_u*((m_a + 2)*pow(l[1], -m_a)*pow(l[0], -m_a-3) + (m_a+2)*pow(l[0], m_a-3))/l[0];
		Cp[1][1] = m_u*((m_a + 2)*pow(l[0], -m_a)*pow(l[1], -m_a-3) + (m_a+2)*pow(l[1], m_a-3))/l[1];
		Cp[0][1] = Cp[1][0] = m_u*m_a*pow(l[0], -m_a-2)*pow(l[1], -m_a-2);

		if (fabs(l[0] - l[1]) < 1e-9)
			ks = (Cp[0][0] - Cp[0][1])/(2*l[0]);
		else
			ks = (S[0] - S[1])/(l[0]*l[0] - l[1]*l[1]);
	}

	 double ct = ec[0], st = ec[1];
	 double T[2][3] = {{ct*ct, st*st, st*ct},{st*st, ct*ct, -st*ct}};

	 double z[3] = {2*ct*st, -2*ct*st, st*st - ct*ct};

	 // transform to element coordinate system
	 for (int i=0; i<3; ++i)
		 for (int j=0; j<3; ++j)
		 {
			 D[i][j] = 0;
			 for (int k=0; k<2; k++)
				 for (int l=0; l<2; ++l)
					 D[i][j] += T[k][i]*Cp[k][l]*T[l][j];
			 D[i][j] += ks*z[i]*z[j];
		 }
}

//-----------------------------------------------------------------------------
void eig(double *A, double *B, double *C)
{
/* calculate eigenvalues A and vectors B of 2x2 matrix C */
double tr,det,t2,pm;
int i,j;
if (((C[1]<1e-6)&&(C[1]>-1e-6))&&((C[2]<1e-6)&&(C[2]>-1e-6))){ /* not equal to zero, plus a tolerance */
A[0]=C[0];A[1]=C[3];
B[0]=1;B[1]=0;B[2]=0;B[3]=1;
}
else{
tr=C[0]+C[3];
det=C[0]*C[3]-C[1]*C[2];
t2=tr/2;
pm=pow((tr*tr/4-det),0.5);
if (pm<1e-4){
A[0]=t2;A[1]=t2;
B[0]=1;B[1]=0;B[2]=0;B[3]=1;
}
else {
A[0]=t2+pm;A[1]=t2-pm;

/* calculate eigenvectors B */
if ((C[1]>1e-9)||(C[1]<-1e-9)){ /* not equal to zero, plus a tolerance for rounding errors */
B[0]=A[0]-C[3];
B[1]=C[1];
B[2]=A[1]-C[3];
B[3]=C[1];
/* normalise eigenvectors to unit length, reusing tr temporarily */
for (i=0;i<2;i++){
   j=i*2;
   tr=pow((B[j]*B[j]+B[j+1]*B[j+1]),0.5);
   B[j]=B[j]/tr;B[j+1]=B[j+1]/tr;
   }
}
else if ((C[2]>1e-9)||(C[2]<-1e-9)){ 
B[0]=C[2];
B[1]=A[0]-C[0];
B[2]=C[2];
B[3]=A[1]-C[0];
/* normalise eigenvectors to unit length, reusing tr temporarily */
for (i=0;i<2;i++){
   j=i*2;
   tr=pow((B[j]*B[j]+B[j+1]*B[j+1]),0.5);
   B[j]=B[j]/tr;B[j+1]=B[j+1]/tr;
   }
}
else {
B[0]=1;B[1]=0;B[2]=0;B[3]=1;
}
}
}
}

#define ROTATE(a, i, j, k, l) g=a[i][j]; h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l] = h + s*(g - h*tau);
void eigen2d(double l[2], double r[4], double m[4])
{
	const int NMAX = 50;
	double sm, tresh, g, h, t, c, tau, s, th;
	int i, j, k;

	// copy the matrix components since we will be overwriting them
	double a[2][2] = {
			{m[0], m[1]},
			{m[2], m[3]}
	};

	// the v matrix contains the eigen vectors
	// intialize to identity
	double v[2][2] = {
		{ 1, 0 },
		{ 0, 1 }
	};

	// initialize b and d to the diagonal of a
	double b[2] = {a[0][0], a[1][1]};
	double d[2] = {a[0][0], a[1][1]};
	double z[2] = {0};

	const double eps = 0;//1.0e-15;

	// loop
	int n, nrot = 0;
	for (n=0; n<NMAX; ++n)
	{
		// sum off-diagonal elements
		sm = fabs(a[0][1]);
		if (sm <= eps) break;

		// set the treshold
		if (n < 2) tresh = 0.2*sm/9.0; else tresh = 0.0;

		// loop over off-diagonal elements
		for (i=0; i<1; ++i)
		{
			for (j=i+1; j<2; ++j)
			{
				g = 100.0*fabs(a[i][j]);

				// after four sweeps, skip the rotation if the off-diagonal element is small
				if ((n > 3) && ((fabs(d[i])+g) == fabs(d[i]))
							&& ((fabs(d[j])+g) == fabs(d[j])))
				{
					a[i][j] = 0.0;
				}
				else if (fabs(a[i][j]) > tresh)
				{
					h = d[j] - d[i];
					if ((fabs(h)+g) == fabs(h))
						t = a[i][j]/h;
					else
					{
						th = 0.5*h/a[i][j];
						t = 1.0/(fabs(th) + sqrt(1+th*th));
						if (th < 0.0) t = -t;
					}

					c = 1.0/sqrt(1.0 + t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[i][j];
					z[i] -= h;
					z[j] += h;
					d[i] -= h;
					d[j] += h;
					a[i][j] = 0;

					for (k=  0; k<=i-1; ++k) { ROTATE(a, k, i, k, j) }
					for (k=i+1; k<=j-1; ++k) { ROTATE(a, i, k, k, j) }
					for (k=j+1; k<   2; ++k) { ROTATE(a, i, k, j, k) }
					for (k=  0; k<   2; ++k) { ROTATE(v, k, i, k, j) }
					++nrot;
				}
			}
		}

		for (i=0; i<2; ++i) 
		{
			b[i] += z[i];
			d[i] = b[i];
			z[i] = 0.0;
		}
	}

	// we sure we converged
	assert(n < NMAX);

	// copy eigenvalues
	l[0] = d[0];
	l[1] = d[1];

	// copy eigenvectors
	if (r)
	{
		r[0] = v[0][0]; r[1] = v[1][0];
		r[2] = v[0][1]; r[3] = v[1][1];
	}
}
