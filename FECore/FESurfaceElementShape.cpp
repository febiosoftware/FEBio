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
#include "FESurfaceElementShape.h"

//=============================================================================
//              Q U A D 4
//=============================================================================

//-----------------------------------------------------------------------------
void FEQuad4::shape_fnc(double* H, double r, double s)
{
	H[0] = 0.25*(1 - r)*(1 - s);
	H[1] = 0.25*(1 + r)*(1 - s);
	H[2] = 0.25*(1 + r)*(1 + s);
	H[3] = 0.25*(1 - r)*(1 + s);
}

//-----------------------------------------------------------------------------
void FEQuad4::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -0.25*(1 - s); Hs[0] = -0.25*(1 - r);
	Hr[1] =  0.25*(1 - s); Hs[1] = -0.25*(1 + r);
	Hr[2] =  0.25*(1 + s); Hs[2] =  0.25*(1 + r);
	Hr[3] = -0.25*(1 + s); Hs[3] =  0.25*(1 - r);
}

//-----------------------------------------------------------------------------
void FEQuad4::shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s)
{
	Hrr[0] = 0; Hrs[0] =  0.25; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = -0.25; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] =  0.25; Hss[2] = 0;
	Hrr[3] = 0; Hrs[3] = -0.25; Hss[3] = 0;
}

//=============================================================================
//              Q U A D 8
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
void FEQuad8::shape_fnc(double* H, double r, double s)
{
	H[4] = 0.5*(1 - r*r)*(1 - s);
	H[5] = 0.5*(1 - s*s)*(1 + r);
	H[6] = 0.5*(1 - r*r)*(1 + s);
	H[7] = 0.5*(1 - s*s)*(1 - r);

	H[0] = 0.25*(1 - r)*(1 - s) - 0.5*(H[4] + H[7]);
	H[1] = 0.25*(1 + r)*(1 - s) - 0.5*(H[4] + H[5]);
	H[2] = 0.25*(1 + r)*(1 + s) - 0.5*(H[5] + H[6]);
	H[3] = 0.25*(1 - r)*(1 + s) - 0.5*(H[6] + H[7]);
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
void FEQuad8::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[4] = -r*(1 - s);
	Hr[5] = 0.5*(1 - s*s);
	Hr[6] = -r*(1 + s);
	Hr[7] = -0.5*(1 - s*s);

	Hr[0] = -0.25*(1 - s) - 0.5*(Hr[4] + Hr[7]);
	Hr[1] = 0.25*(1 - s) - 0.5*(Hr[4] + Hr[5]);
	Hr[2] = 0.25*(1 + s) - 0.5*(Hr[5] + Hr[6]);
	Hr[3] = -0.25*(1 + s) - 0.5*(Hr[6] + Hr[7]);

	Hs[4] = -0.5*(1 - r*r);
	Hs[5] = -s*(1 + r);
	Hs[6] = 0.5*(1 - r*r);
	Hs[7] = -s*(1 - r);

	Hs[0] = -0.25*(1 - r) - 0.5*(Hs[4] + Hs[7]);
	Hs[1] = -0.25*(1 + r) - 0.5*(Hs[4] + Hs[5]);
	Hs[2] = 0.25*(1 + r) - 0.5*(Hs[5] + Hs[6]);
	Hs[3] = 0.25*(1 - r) - 0.5*(Hs[6] + Hs[7]);
}

//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
//! \todo implement this
void FEQuad8::shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s)
{
	Hrr[4] = -(1 - s);
	Hrr[5] = 0.0;
	Hrr[6] = -(1 + s);
	Hrr[7] = 0.0;

	Hrs[4] = r;
	Hrs[5] = -s;
	Hrs[6] = -r;
	Hrs[7] = s;

	Hss[4] = 0.0;
	Hss[5] = -(1 + r);
	Hss[6] = 0.0;
	Hss[7] = -(1 - r);

	Hrr[0] = -0.5*(Hrr[4] + Hrr[7]);
	Hrr[1] = -0.5*(Hrr[4] + Hrr[5]);
	Hrr[2] = -0.5*(Hrr[5] + Hrr[6]);
	Hrr[3] = -0.5*(Hrr[6] + Hrr[7]);

	Hrs[0] = 0.25 - 0.5*(Hrs[4] + Hrs[7]);
	Hrs[1] = -0.25 - 0.5*(Hrs[4] + Hrs[5]);
	Hrs[2] = 0.25 - 0.5*(Hrs[5] + Hrs[6]);
	Hrs[3] = -0.25 - 0.5*(Hrs[6] + Hrs[7]);

	Hss[0] = -0.5*(Hss[4] + Hss[7]);
	Hss[1] = -0.5*(Hss[4] + Hss[5]);
	Hss[2] = -0.5*(Hss[5] + Hss[6]);
	Hss[3] = -0.5*(Hss[6] + Hss[7]);
}

//=============================================================================
//              Q U A D 9
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
void FEQuad9::shape_fnc(double* H, double r, double s)
{
	double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
	double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };

	H[0] = R[0] * S[0];
	H[1] = R[1] * S[0];
	H[2] = R[1] * S[1];
	H[3] = R[0] * S[1];
	H[4] = R[2] * S[0];
	H[5] = R[1] * S[2];
	H[6] = R[2] * S[1];
	H[7] = R[0] * S[2];
	H[8] = R[2] * S[2];
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
void FEQuad9::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
	double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
	double DR[3] = { r - 0.5, r + 0.5, -2.0*r };
	double DS[3] = { s - 0.5, s + 0.5, -2.0*s };

	Hr[0] = DR[0] * S[0];
	Hr[1] = DR[1] * S[0];
	Hr[2] = DR[1] * S[1];
	Hr[3] = DR[0] * S[1];
	Hr[4] = DR[2] * S[0];
	Hr[5] = DR[1] * S[2];
	Hr[6] = DR[2] * S[1];
	Hr[7] = DR[0] * S[2];
	Hr[8] = DR[2] * S[2];

	Hs[0] = R[0] * DS[0];
	Hs[1] = R[1] * DS[0];
	Hs[2] = R[1] * DS[1];
	Hs[3] = R[0] * DS[1];
	Hs[4] = R[2] * DS[0];
	Hs[5] = R[1] * DS[2];
	Hs[6] = R[2] * DS[1];
	Hs[7] = R[0] * DS[2];
	Hs[8] = R[2] * DS[2];
}

//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
void FEQuad9::shape_deriv2(double* Grr, double* Gss, double* Grs, double r, double s)
{
	double R[3] = { 0.5*r*(r - 1.0), 0.5*r*(r + 1.0), 1.0 - r*r };
	double S[3] = { 0.5*s*(s - 1.0), 0.5*s*(s + 1.0), 1.0 - s*s };
	double DR[3] = { r - 0.5, r + 0.5, -2.0*r };
	double DS[3] = { s - 0.5, s + 0.5, -2.0*s };
	double DDR[3] = { 1.0, 1.0, -2.0 };
	double DDS[3] = { 1.0, 1.0, -2.0 };

	Grr[0] = DDR[0] * S[0]; Grs[0] = DR[0] * DS[0]; Gss[0] = R[0] * DDS[0];
	Grr[1] = DDR[1] * S[0]; Grs[1] = DR[1] * DS[0]; Gss[1] = R[1] * DDS[0];
	Grr[2] = DDR[1] * S[1]; Grs[2] = DR[1] * DS[1]; Gss[2] = R[1] * DDS[1];
	Grr[3] = DDR[0] * S[1]; Grs[3] = DR[0] * DS[1]; Gss[3] = R[0] * DDS[1];
	Grr[4] = DDR[2] * S[0]; Grs[4] = DR[2] * DS[0]; Gss[4] = R[2] * DDS[0];
	Grr[5] = DDR[1] * S[2]; Grs[5] = DR[1] * DS[2]; Gss[5] = R[1] * DDS[2];
	Grr[6] = DDR[2] * S[1]; Grs[6] = DR[2] * DS[1]; Gss[6] = R[2] * DDS[1];
	Grr[7] = DDR[0] * S[2]; Grs[7] = DR[0] * DS[2]; Gss[7] = R[0] * DDS[2];
	Grr[8] = DDR[2] * S[2]; Grs[8] = DR[2] * DS[2]; Gss[8] = R[2] * DDS[2];
}

//=============================================================================
//              T R I 3
//=============================================================================

//-----------------------------------------------------------------------------
void FETri3::shape_fnc(double* H, double r, double s)
{
	H[0] = 1.0 - r - s;
	H[1] = r;
	H[2] = s;
}

//-----------------------------------------------------------------------------
void FETri3::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -1; Hs[0] = -1;
	Hr[1] =  1; Hs[1] =  0;
	Hr[2] =  0; Hs[2] =  1;
}

//-----------------------------------------------------------------------------
void FETri3::shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s)
{
	Hrr[0] = 0; Hrs[0] = 0; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = 0; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] = 0; Hss[2] = 0;
}

//=============================================================================
//              T R I 6
//=============================================================================

//-----------------------------------------------------------------------------
void FETri6::shape_fnc(double* H, double r, double s)
{
	double r1 = 1.0 - r - s;
	double r2 = r;
	double r3 = s;

	H[0] = r1*(2.0*r1 - 1.0);
	H[1] = r2*(2.0*r2 - 1.0);
	H[2] = r3*(2.0*r3 - 1.0);
	H[3] = 4.0*r1*r2;
	H[4] = 4.0*r2*r3;
	H[5] = 4.0*r3*r1;
}

//-----------------------------------------------------------------------------
void FETri6::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -3.0 + 4.0*r + 4.0*s;
	Hr[1] = 4.0*r - 1.0;
	Hr[2] = 0.0;
	Hr[3] = 4.0 - 8.0*r - 4.0*s;
	Hr[4] = 4.0*s;
	Hr[5] = -4.0*s;

	Hs[0] = -3.0 + 4.0*s + 4.0*r;
	Hs[1] = 0.0;
	Hs[2] = 4.0*s - 1.0;
	Hs[3] = -4.0*r;
	Hs[4] = 4.0*r;
	Hs[5] = 4.0 - 8.0*s - 4.0*r;
}

//-----------------------------------------------------------------------------
void FETri6::shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s)
{
	Hrr[0] =  4.0; Hrs[0] =  4.0; Hss[0] =  4.0;
	Hrr[1] =  4.0; Hrs[1] =  0.0; Hss[1] =  0.0;
	Hrr[2] =  0.0; Hrs[2] =  0.0; Hss[2] =  4.0;
	Hrr[3] = -8.0; Hrs[3] = -4.0; Hss[3] =  0.0;
	Hrr[4] =  0.0; Hrs[4] =  4.0; Hss[4] =  0.0;
	Hrr[5] =  0.0; Hrs[5] = -4.0; Hss[5] = -8.0;
}

//=============================================================================
//              T R I 7
//=============================================================================

//-----------------------------------------------------------------------------
void FETri7::shape_fnc(double* H, double r, double s)
{
	double r1 = 1.0 - r - s;
	double r2 = r;
	double r3 = s;

	H[6] = 27.0*r1*r2*r3;
	H[0] = r1*(2.0*r1 - 1.0) + H[6] / 9.0;
	H[1] = r2*(2.0*r2 - 1.0) + H[6] / 9.0;
	H[2] = r3*(2.0*r3 - 1.0) + H[6] / 9.0;
	H[3] = 4.0*r1*r2 - 4.0*H[6] / 9.0;
	H[4] = 4.0*r2*r3 - 4.0*H[6] / 9.0;
	H[5] = 4.0*r3*r1 - 4.0*H[6] / 9.0;
}

//-----------------------------------------------------------------------------
void FETri7::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[6] = 27.0*s*(1.0 - 2.0*r - s);
	Hr[0] = -3.0 + 4.0*r + 4.0*s + Hr[6] / 9.0;
	Hr[1] = 4.0*r - 1.0 + Hr[6] / 9.0;
	Hr[2] = 0.0 + Hr[6] / 9.0;
	Hr[3] = 4.0 - 8.0*r - 4.0*s - 4.0*Hr[6] / 9.0;
	Hr[4] = 4.0*s - 4.0*Hr[6] / 9.0;
	Hr[5] = -4.0*s - 4.0*Hr[6] / 9.0;

	Hs[6] = 27.0*r*(1.0 - r - 2.0*s);
	Hs[0] = -3.0 + 4.0*s + 4.0*r + Hs[6] / 9.0;
	Hs[1] = 0.0 + Hs[6] / 9.0;
	Hs[2] = 4.0*s - 1.0 + Hs[6] / 9.0;
	Hs[3] = -4.0*r - 4.0*Hs[6] / 9.0;
	Hs[4] = 4.0*r - 4.0*Hs[6] / 9.0;
	Hs[5] = 4.0 - 8.0*s - 4.0*r - 4.0*Hs[6] / 9.0;
}

//-----------------------------------------------------------------------------
void FETri7::shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s)
{
	Hrr[6] = -54.0*s;
	Hss[6] = -54.0*r;
	Hrs[6] = 27.0*(1.0 - 2.0*r - 2.0*s);

	Hrr[0] = 4.0 + Hrr[6] / 9.0; Hrs[0] = 4.0 + Hrs[6] / 9.0; Hss[0] = 4.0 + Hss[6] / 9.0;
	Hrr[1] = 4.0 + Hrr[6] / 9.0; Hrs[1] = 0.0 + Hrs[6] / 9.0; Hss[1] = 0.0 + Hss[6] / 9.0;
	Hrr[2] = 0.0 + Hrr[6] / 9.0; Hrs[2] = 0.0 + Hrs[6] / 9.0; Hss[2] = 4.0 + Hss[6] / 9.0;
	Hrr[3] = -8.0 - 4.0*Hrr[6] / 9.0; Hrs[3] = -4.0 - 4.0*Hrs[6] / 9.0; Hss[3] = 0.0 - 4.0*Hss[6] / 9.0;
	Hrr[4] = 0.0 - 4.0*Hrr[6] / 9.0; Hrs[4] = 4.0 - 4.0*Hrs[6] / 9.0; Hss[4] = 0.0 - 4.0*Hss[6] / 9.0;
	Hrr[5] = 0.0 - 4.0*Hrr[6] / 9.0; Hrs[5] = -4.0 - 4.0*Hrs[6] / 9.0; Hss[5] = -8.0 - 4.0*Hss[6] / 9.0;
}

//=============================================================================
//              T R I 10
//=============================================================================

//-----------------------------------------------------------------------------
void FETri10::shape_fnc(double* H, double r, double s)
{
	double L1 = 1.0 - r - s;
	double L2 = r;
	double L3 = s;

	H[0] = 0.5*(3 * L1 - 1)*(3 * L1 - 2)*L1;
	H[1] = 0.5*(3 * L2 - 1)*(3 * L2 - 2)*L2;
	H[2] = 0.5*(3 * L3 - 1)*(3 * L3 - 2)*L3;
	H[3] = 4.5*(3 * L1 - 1)*L1*L2;
	H[4] = 4.5*(3 * L2 - 1)*L1*L2;
	H[5] = 4.5*(3 * L2 - 1)*L2*L3;
	H[6] = 4.5*(3 * L3 - 1)*L2*L3;
	H[7] = 4.5*(3 * L1 - 1)*L1*L3;
	H[8] = 4.5*(3 * L3 - 1)*L1*L3;
	H[9] = 27.*L1*L2*L3;
}

//-----------------------------------------------------------------------------
void FETri10::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	double L1 = 1.0 - r - s;
	double L2 = r;
	double L3 = s;

	Hr[0] = -3. / 2.*(3 * L1 - 2)*L1 - 3. / 2.*(3 * L1 - 1)*L1 - 0.5*(3 * L1 - 1)*(3 * L1 - 2);
	Hr[1] = 3. / 2.*(3 * L2 - 2)*L2 + 3. / 2.*(3 * L2 - 1)*L2 + 0.5*(3 * L2 - 1)*(3 * L2 - 2);
	Hr[2] = 0.0;
	Hr[3] = -27. / 2.*L1*L2 - 9. / 2.*(3 * L1 - 1)*L2 + 9. / 2.*(3 * L1 - 1)*L1;
	Hr[4] = 27. / 2.*L1*L2 - 9. / 2.*(3 * L2 - 1)*L2 + 9. / 2.*(3 * L2 - 1)*L1;
	Hr[5] = 27. / 2.*L2*L3 + 9. / 2.*(3 * L2 - 1)*L3;
	Hr[6] = 9. / 2.*(3 * L3 - 1)*L3;
	Hr[7] = -27. / 2.*L1*L3 - 9. / 2.*(3 * L1 - 1)*L3;
	Hr[8] = -9. / 2.*(3 * L3 - 1)*L3;
	Hr[9] = -27.*L2*L3 + 27.*L1*L3;

	Hs[0] = -3. / 2.*(3 * L1 - 2)*L1 - 3. / 2.*(3 * L1 - 1)*L1 - 0.5*(3 * L1 - 1)*(3 * L1 - 2);
	Hs[1] = 0.0;
	Hs[2] = 3. / 2.*(3 * L3 - 2)*L3 + 3. / 2.*(3 * L3 - 1)*L3 + 0.5*(3 * L3 - 1)*(3 * L3 - 2);
	Hs[3] = -27. / 2.*L1*L2 - 9. / 2.*(3 * L1 - 1)*L2;
	Hs[4] = -9. / 2.*(3 * L2 - 1)*L2;
	Hs[5] = 9. / 2.*(3 * L2 - 1)*L2;
	Hs[6] = 27. / 2.*L2*L3 + 9. / 2.*(3 * L3 - 1)*L2;
	Hs[7] = -27. / 2.*L1*L3 - 9. / 2.*(3 * L1 - 1)*L3 + 9. / 2.*(3 * L1 - 1)*L1;
	Hs[8] = 27. / 2.*L1*L3 - 9. / 2.*(3 * L3 - 1)*L3 + 9. / 2.*(3 * L3 - 1)*L1;
	Hs[9] = -27.*L2*L3 + 27.*L1*L2;
}

//-----------------------------------------------------------------------------
void FETri10::shape_deriv2(double* Hrr, double* Hss, double* Hrs, double r, double s)
{
	// TODO: Implement this
}
