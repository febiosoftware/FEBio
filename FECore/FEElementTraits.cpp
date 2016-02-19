// FEElementTraits.cpp: implementation of the FEElementTraits class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElementTraits.h"
#include "FEElement.h"
#include "FEException.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

FEElementTraits::FEElementTraits(int ni, int ne, FE_Element_Class c, FE_Element_Shape s, FE_Element_Type t)
{
	neln = ne;
	nint = ni;
	spec.eclass = c;
	spec.eshape = s;
	spec.etype  = t;
	H.resize(ni, ne);
}

//=============================================================================
FESolidElementTraits::FESolidElementTraits(int ni, int ne, FE_Element_Shape eshape, FE_Element_Type etype) : FEElementTraits(ni, ne, FE_ELEM_SOLID, eshape, etype) 
{
	gr.resize(ni);
	gs.resize(ni);
	gt.resize(ni);
	gw.resize(ni);

	Gr.resize(ni, ne);
	Gs.resize(ni, ne);
	Gt.resize(ni, ne);

	Grr.resize(ni, ne);
	Gsr.resize(ni, ne);
	Gtr.resize(ni, ne);
		
	Grs.resize(ni, ne);
	Gss.resize(ni, ne);
	Gts.resize(ni, ne);
		
	Grt.resize(ni, ne);
	Gst.resize(ni, ne);
	Gtt.resize(ni, ne);
}

//-----------------------------------------------------------------------------
//! initialize element traits data
void FESolidElementTraits::init()
{
	assert(nint > 0);
	assert(neln > 0);
	const int NELN = FEElement::MAX_NODES;

	// calculate shape function values at gauss points
	double N[NELN];
	for (int n=0; n<nint; ++n)
	{
		shape_fnc(N, gr[n], gs[n], gt[n]);
		for (int i=0; i<neln; ++i) H[n][i] = N[i];
	}

	// calculate local derivatives of shape functions at gauss points
	double Hr[NELN], Hs[NELN], Ht[NELN];
	for (int n=0; n<nint; ++n)
	{
		shape_deriv(Hr, Hs, Ht, gr[n], gs[n], gt[n]);
		for (int i=0; i<neln; ++i)
		{
			Gr[n][i] = Hr[i];
			Gs[n][i] = Hs[i];
			Gt[n][i] = Ht[i];
		}
	}
	
	// calculate local second derivatives of shape functions at gauss points
	double Hrr[NELN], Hss[NELN], Htt[NELN], Hrs[NELN], Hst[NELN], Hrt[NELN];
	for (int n=0; n<nint; ++n)
	{
		shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, gr[n], gs[n], gt[n]);
		for (int i=0; i<neln; ++i)
		{
			Grr[n][i] = Hrr[i]; Grs[n][i] = Hrs[i]; Grt[n][i] = Hrt[i]; 
			Gsr[n][i] = Hrs[i]; Gss[n][i] = Hss[i]; Gst[n][i] = Hst[i]; 
			Gtr[n][i] = Hrt[i]; Gts[n][i] = Hst[i]; Gtt[n][i] = Htt[i]; 
		}
	}
}

//=============================================================================
//                                F E H E X 8
//=============================================================================

//-----------------------------------------------------------------------------
void FEHex8_::shape_fnc(double* H, double r, double s, double t)
{
	H[0] = 0.125*(1 - r)*(1 - s)*(1 - t);
	H[1] = 0.125*(1 + r)*(1 - s)*(1 - t);
	H[2] = 0.125*(1 + r)*(1 + s)*(1 - t);
	H[3] = 0.125*(1 - r)*(1 + s)*(1 - t);
	H[4] = 0.125*(1 - r)*(1 - s)*(1 + t);
	H[5] = 0.125*(1 + r)*(1 - s)*(1 + t);
	H[6] = 0.125*(1 + r)*(1 + s)*(1 + t);
	H[7] = 0.125*(1 - r)*(1 + s)*(1 + t);
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FEHex8_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[0] = -0.125*(1 - s)*(1 - t);
	Hr[1] =  0.125*(1 - s)*(1 - t);
	Hr[2] =  0.125*(1 + s)*(1 - t);
	Hr[3] = -0.125*(1 + s)*(1 - t);
	Hr[4] = -0.125*(1 - s)*(1 + t);
	Hr[5] =  0.125*(1 - s)*(1 + t);
	Hr[6] =  0.125*(1 + s)*(1 + t);
	Hr[7] = -0.125*(1 + s)*(1 + t);
		
	Hs[0] = -0.125*(1 - r)*(1 - t);
	Hs[1] = -0.125*(1 + r)*(1 - t);
	Hs[2] =  0.125*(1 + r)*(1 - t);
	Hs[3] =  0.125*(1 - r)*(1 - t);
	Hs[4] = -0.125*(1 - r)*(1 + t);
	Hs[5] = -0.125*(1 + r)*(1 + t);
	Hs[6] =  0.125*(1 + r)*(1 + t);
	Hs[7] =  0.125*(1 - r)*(1 + t);
		
	Ht[0] = -0.125*(1 - r)*(1 - s);
	Ht[1] = -0.125*(1 + r)*(1 - s);
	Ht[2] = -0.125*(1 + r)*(1 + s);
	Ht[3] = -0.125*(1 - r)*(1 + s);
	Ht[4] =  0.125*(1 - r)*(1 - s);
	Ht[5] =  0.125*(1 + r)*(1 - s);
	Ht[6] =  0.125*(1 + r)*(1 + s);
	Ht[7] =  0.125*(1 - r)*(1 + s);
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
void FEHex8_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	Hrr[0] = 0.0; Hss[0] = 0.0; Htt[0] = 0.0;
	Hrr[1] = 0.0; Hss[1] = 0.0; Htt[1] = 0.0;
	Hrr[2] = 0.0; Hss[2] = 0.0; Htt[2] = 0.0;
	Hrr[3] = 0.0; Hss[3] = 0.0; Htt[3] = 0.0;
	Hrr[4] = 0.0; Hss[4] = 0.0; Htt[4] = 0.0;
	Hrr[5] = 0.0; Hss[5] = 0.0; Htt[5] = 0.0;
	Hrr[6] = 0.0; Hss[6] = 0.0; Htt[6] = 0.0;
	Hrr[7] = 0.0; Hss[7] = 0.0; Htt[7] = 0.0;
		
	Hrs[0] =  0.125*(1 - t);
	Hrs[1] = -0.125*(1 - t);
	Hrs[2] =  0.125*(1 - t);
	Hrs[3] = -0.125*(1 - t);
	Hrs[4] =  0.125*(1 + t);
	Hrs[5] = -0.125*(1 + t);
	Hrs[6] =  0.125*(1 + t);
	Hrs[7] = -0.125*(1 + t);
		
	Hrt[0] =  0.125*(1 - s);
	Hrt[1] = -0.125*(1 - s);
	Hrt[2] = -0.125*(1 + s);
	Hrt[3] =  0.125*(1 + s);
	Hrt[4] = -0.125*(1 - s);
	Hrt[5] =  0.125*(1 - s);
	Hrt[6] =  0.125*(1 + s);
	Hrt[7] = -0.125*(1 + s);
		
	Hst[0] =  0.125*(1 - r);
	Hst[1] =  0.125*(1 + r);
	Hst[2] = -0.125*(1 + r);
	Hst[3] = -0.125*(1 - r);
	Hst[4] = -0.125*(1 - r);
	Hst[5] = -0.125*(1 + r);
	Hst[6] =  0.125*(1 + r);
	Hst[7] =  0.125*(1 - r);
}

//*****************************************************************************
//                          H E X 8 G 8 
//*****************************************************************************

FEHex8G8::FEHex8G8() : FEHex8_(NINT, FE_HEX8G8)
{
	// integration point coordinates
	const double a = 1.0 / sqrt(3.0);
	gr[0] = -a; gs[0] = -a; gt[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gt[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gt[2] = -a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gt[3] = -a; gw[3] = 1;
	gr[4] = -a; gs[4] = -a; gt[4] =  a; gw[4] = 1;
	gr[5] =  a; gs[5] = -a; gt[5] =  a; gw[5] = 1;
	gr[6] =  a; gs[6] =  a; gt[6] =  a; gw[6] = 1;
	gr[7] = -a; gs[7] =  a; gt[7] =  a; gw[7] = 1;
	init();
	Hi = H.inverse();
}

//-----------------------------------------------------------------------------
void FEHex8G8::project_to_nodes(double* ai, double* ao)
{
	for (int j=0; j<NELN; ++j)
	{
		ao[j] = 0;
		for (int k=0; k<NINT; ++k) 
		{
			ao[j] += Hi[j][k]*ai[k];
		}
	}
}

//*****************************************************************************
//                          F E H E X R I 
//*****************************************************************************

FEHex8RI::FEHex8RI(): FEHex8_(NINT, FE_HEX8RI)
{
	// This is for a six point integration rule
	// integration point coordinates
	const double a = 8.0 / 6.0;
	gr[0] = -1; gs[0] = 0; gt[0] = 0; gw[0] = a;
	gr[1] =  1; gs[1] = 0; gt[1] = 0; gw[1] = a;
	gr[2] =  0; gs[2] =-1; gt[2] = 0; gw[2] = a;
	gr[3] =  0; gs[3] = 1; gt[3] = 0; gw[3] = a;
	gr[4] =  0; gs[4] = 0; gt[4] =-1; gw[4] = a;
	gr[5] =  0; gs[5] = 0; gt[5] = 1; gw[5] = a;
	
	init();
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FEHex8RI::project_to_nodes(double* ai, double* ao)
{
	
}

//*****************************************************************************
//                          F E H E X G 1
//*****************************************************************************

FEHex8G1::FEHex8G1() : FEHex8_(NINT, FE_HEX8G1)
{
	// single gauss-point integration rule
	gr[0] = 0; gs[0] = 0; gt[0] = 0; gw[0] = 8.0;
	init();
}

//-----------------------------------------------------------------------------
void FEHex8G1::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
	ao[3] = ai[0];
	ao[4] = ai[0];
	ao[5] = ai[0];
	ao[6] = ai[0];
	ao[7] = ai[0];
}

//=============================================================================
//                              F E T E T 4
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FETet4_::shape_fnc(double* H, double r, double s, double t)
{
	H[0] = 1 - r - s - t;
	H[1] = r;
	H[2] = s;
	H[3] = t;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FETet4_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[0] = -1; Hs[0] = -1; Ht[0] = -1;
	Hr[1] =  1;	Hs[1] =  0; Ht[1] =  0;
	Hr[2] =  0;	Hs[2] =  1; Ht[2] =  0;
	Hr[3] =  0;	Hs[3] =  0; Ht[3] =  1;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
void FETet4_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	Hrr[0] =  0.0; Hss[0] =  0.0; Htt[0] =  0.0;
	Hrr[1] =  0.0; Hss[1] =  0.0; Htt[1] =  0.0;
	Hrr[2] =  0.0; Hss[2] =  0.0; Htt[2] =  0.0;
	Hrr[3] =  0.0; Hss[3] =  0.0; Htt[3] =  0.0;

	Hrs[0] =  0.0; Hst[0] =  0.0; Hrt[0] =  0.0;
	Hrs[1] =  0.0; Hst[1] =  0.0; Hrt[1] =  0.0;
	Hrs[2] =  0.0; Hst[2] =  0.0; Hrt[2] =  0.0;
	Hrs[3] =  0.0; Hst[3] =  0.0; Hrt[3] =  0.0;
}

//=============================================================================
//                          T E T 4
//=============================================================================

FETet4G4::FETet4G4() : FETet4_(NINT, FE_TET4G4)
{
	// gaussian integration for tetrahedral elements
	const double a = 0.58541020;
	const double b = 0.13819660;
	const double w = 1.0 / 24.0;
	
	gr[0] = b; gs[0] = b; gt[0] = b; gw[0] = w;
	gr[1] = a; gs[1] = b; gt[1] = b; gw[1] = w;
	gr[2] = b; gs[2] = a; gt[2] = b; gw[2] = w;
	gr[3] = b; gs[3] = b; gt[3] = a; gw[3] = w;

	init();
	Hi = H.inverse();
}

//-----------------------------------------------------------------------------
void FETet4G4::project_to_nodes(double* ai, double* ao)
{
	for (int j=0; j<NELN; ++j)
	{
		ao[j] = 0;
		for (int k=0; k<NINT; ++k) 
		{
			ao[j] += Hi[j][k]*ai[k];
		}
	}
}

//=============================================================================
//                          F E G 1 T E T E L E M E N T
//=============================================================================

FETet4G1::FETet4G1() : FETet4_(NINT, FE_TET4G1)
{
	// gaussian integration for tetrahedral elements
	const double a = 0.25;
	const double w = 1.0 / 6.0;
	
	gr[0] = a; gs[0] = a; gt[0] = a; gw[0] = w;
	init();
}

//-----------------------------------------------------------------------------
void FETet4G1::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
	ao[3] = ai[0];
}

//=============================================================================
//                       P E N T A 6
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FEPenta6_::shape_fnc(double* H, double r, double s, double t)
{
	H[0] = 0.5*(1 - t)*(1 - r - s);
	H[1] = 0.5*(1 - t)*r;
	H[2] = 0.5*(1 - t)*s;
	H[3] = 0.5*(1 + t)*(1 - r - s);
	H[4] = 0.5*(1 + t)*r;
	H[5] = 0.5*(1 + t)*s;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FEPenta6_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[0] = -0.5*(1 - t);
	Hr[1] =  0.5*(1 - t);
	Hr[2] =  0.0;
	Hr[3] = -0.5*(1 + t);
	Hr[4] =  0.5*(1 + t);
	Hr[5] =  0.0;
		
	Hs[0] = -0.5*(1 - t);
	Hs[1] =  0.0;
	Hs[2] =  0.5*(1 - t);
	Hs[3] = -0.5*(1 + t);
	Hs[4] =  0.0;
	Hs[5] =  0.5*(1 + t);
		
	Ht[0] = -0.5*(1 - r - s);
	Ht[1] = -0.5*r;
	Ht[2] = -0.5*s;
	Ht[3] =  0.5*(1 - r - s);
	Ht[4] =  0.5*r;
	Ht[5] =  0.5*s;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
void FEPenta6_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	Hrr[0] =  0.0; Hss[0] =  0.0; Htt[0] =  0.0;
	Hrr[1] =  0.0; Hss[1] =  0.0; Htt[1] =  0.0;
	Hrr[2] =  0.0; Hss[2] =  0.0; Htt[2] =  0.0;
	Hrr[3] =  0.0; Hss[3] =  0.0; Htt[3] =  0.0;
	Hrr[4] =  0.0; Hss[4] =  0.0; Htt[4] =  0.0;
	Hrr[5] =  0.0; Hss[5] =  0.0; Htt[5] =  0.0;
		
	Hrs[0] =  0.0; Hst[0] =  0.5; Hrt[0] =  0.5;
	Hrs[1] =  0.0; Hst[1] =  0.0; Hrt[1] = -0.5;
	Hrs[2] =  0.0; Hst[2] = -0.5; Hrt[2] =  0.0;
	Hrs[3] =  0.0; Hst[3] = -0.5; Hrt[3] = -0.5;
	Hrs[4] =  0.0; Hst[4] =  0.0; Hrt[4] =  0.5;
	Hrs[5] =  0.0; Hst[5] =  0.5; Hrt[5] =  0.0;
}

//=============================================================================
//                         P E N T A 6 G 6
//=============================================================================

FEPenta6G6::FEPenta6G6(): FEPenta6_(NINT, FE_PENTA6G6)
{
	//gauss intergration points
	const double a = 1.0/6.0;
	const double b = 2.0/3.0;
	const double c = 1.0 / sqrt(3.0);
	const double w = 1.0 / 6.0;
	
	gr[0] = a; gs[0] = a; gt[0] = -c; gw[0] = w;
	gr[1] = b; gs[1] = a; gt[1] = -c; gw[1] = w;
	gr[2] = a; gs[2] = b; gt[2] = -c; gw[2] = w;
	gr[3] = a; gs[3] = a; gt[3] =  c; gw[3] = w;
	gr[4] = b; gs[4] = a; gt[4] =  c; gw[4] = w;
	gr[5] = a; gs[5] = b; gt[5] =  c; gw[5] = w;

	init();

	Hi = H.inverse();
}

//-----------------------------------------------------------------------------
void FEPenta6G6::project_to_nodes(double* ai, double* ao)
{
	for (int j=0; j<NELN; ++j)
	{
		ao[j] = 0;
		for (int k=0; k<NINT; ++k) 
		{
			ao[j] += Hi[j][k]*ai[k];
		}
	}
}

//=============================================================================
//                           T E T 1 0
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FETet10_::shape_fnc(double* H, double r, double s, double t)
{
	double r1 = 1.0 - r - s - t;
	double r2 = r;
	double r3 = s;
	double r4 = t;

	H[0] = r1*(2.0*r1 - 1.0);
	H[1] = r2*(2.0*r2 - 1.0);
	H[2] = r3*(2.0*r3 - 1.0);
	H[3] = r4*(2.0*r4 - 1.0);
	H[4] = 4.0*r1*r2;
	H[5] = 4.0*r2*r3;
	H[6] = 4.0*r3*r1;
	H[7] = 4.0*r1*r4;
	H[8] = 4.0*r2*r4;
	H[9] = 4.0*r3*r4;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FETet10_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[0] = -3.0 + 4.0*r + 4.0*(s + t);
	Hr[1] =  4.0*r - 1.0;
	Hr[2] =  0.0;
	Hr[3] =  0.0;
	Hr[4] =  4.0 - 8.0*r - 4.0*(s + t);
	Hr[5] =  4.0*s;
	Hr[6] = -4.0*s;
	Hr[7] = -4.0*t;
	Hr[8] =  4.0*t;
	Hr[9] =  0.0;

	Hs[0] = -3.0 + 4.0*s + 4.0*(r + t);
	Hs[1] =  0.0;
	Hs[2] =  4.0*s - 1.0;
	Hs[3] =  0.0;
	Hs[4] = -4.0*r;
	Hs[5] =  4.0*r;
	Hs[6] =  4.0 - 8.0*s - 4.0*(r + t);
	Hs[7] = -4.0*t;
	Hs[8] =  0.0;
	Hs[9] =  4.0*t;

	Ht[0] = -3.0 + 4.0*t + 4.0*(r + s);
	Ht[1] =  0.0;
	Ht[2] =  0.0;
	Ht[3] =  4.0*t - 1.0;
	Ht[4] = -4.0*r;
	Ht[5] =  0.0;
	Ht[6] = -4.0*s;
	Ht[7] =  4.0 - 8.0*t - 4.0*(r + s);
	Ht[8] =  4.0*r;
	Ht[9] =  4.0*s;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
void FETet10_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	Hrr[0] =  4.0; Hss[0] =  4.0; Htt[0] =  4.0;
	Hrr[1] =  4.0; Hss[1] =  0.0; Htt[1] =  0.0;
	Hrr[2] =  0.0; Hss[2] =  4.0; Htt[2] =  0.0;
	Hrr[3] =  0.0; Hss[3] =  0.0; Htt[3] =  4.0;
	Hrr[4] = -8.0; Hss[4] =  0.0; Htt[4] =  0.0;
	Hrr[5] =  0.0; Hss[5] =  0.0; Htt[5] =  0.0;
	Hrr[6] =  0.0; Hss[6] = -8.0; Htt[6] =  0.0;
	Hrr[7] =  0.0; Hss[7] =  0.0; Htt[7] =  8.0;
	Hrr[8] =  0.0; Hss[8] =  0.0; Htt[8] =  0.0;
	Hrr[9] =  0.0; Hss[9] =  0.0; Htt[9] =  0.0;

	Hrs[0] =  4.0; Hst[0] =  4.0; Hrt[0] =  4.0;
	Hrs[1] =  0.0; Hst[1] =  0.0; Hrt[1] =  0.0;
	Hrs[2] =  0.0; Hst[2] =  0.0; Hrt[2] =  0.0;
	Hrs[3] =  0.0; Hst[3] =  0.0; Hrt[3] =  0.0;
	Hrs[4] = -4.0; Hst[4] =  0.0; Hrt[4] = -4.0;
	Hrs[5] =  4.0; Hst[5] =  0.0; Hrt[5] =  0.0;
	Hrs[6] = -4.0; Hst[6] = -4.0; Hrt[6] =  0.0;
	Hrs[7] =  0.0; Hst[7] = -4.0; Hrt[7] = -4.0;
	Hrs[8] =  0.0; Hst[8] =  0.0; Hrt[8] =  4.0;
	Hrs[9] =  0.0; Hst[9] =  4.0; Hrt[9] =  0.0;
}

//=============================================================================
//                          T E T 1 0 G 1
//=============================================================================

FETet10G1::FETet10G1() : FETet10_(NINT, FE_TET10G1)
{
	// gaussian integration for tetrahedral elements
	const double a = 0.25;
	const double w = 1.0 / 6.0;
	
	gr[0] = a; gs[0] = a; gt[0] = a; gw[0] = w;
	init();
}

//-----------------------------------------------------------------------------
void FETet10G1::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
	ao[3] = ai[0];
	ao[4] = ai[0];
	ao[5] = ai[0];
	ao[6] = ai[0];
	ao[7] = ai[0];
	ao[8] = ai[0];
	ao[9] = ai[0];
}

//*****************************************************************************
//                          F E T E T 1 0 E L E M E N T
//*****************************************************************************
// I think this assumes that the tetrahedron in natural space is essentially
// a constant metric tet with the edge nodes at the center of the edges.
FETet10G4::FETet10G4() : FETet10_(NINT, FE_TET10G4)
{
	// integration point coordinates
	const double a = 0.58541020;
	const double b = 0.13819660;
	const double w = 0.25 / 6.0;
	gr[ 0] = a; gs[ 0] = b; gt[ 0] = b; gw[ 0] = w;
	gr[ 1] = b; gs[ 1] = a; gt[ 1] = b; gw[ 1] = w;
	gr[ 2] = b; gs[ 2] = b; gt[ 2] = a; gw[ 2] = w;
	gr[ 3] = b; gs[ 3] = b; gt[ 3] = b; gw[ 3] = w;

	init();

	// setup the shape function matrix
	matrix A(4,4);
	for (int i=0; i<4; ++i)
	{
		double r = gr[i];
		double s = gs[i];
		double t = gt[i];

		A[i][0] = 1 - r - s - t;
		A[i][1] = r;
		A[i][2] = s;
		A[i][3] = t;
	}

	// calculate inverse matrix
	Ai.resize(4, 4);
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETet10G4::project_to_nodes(double* ai, double* ao)
{
	ao[0] = Ai[0][0]*ai[0] + Ai[0][1]*ai[1] + Ai[0][2]*ai[2] + Ai[0][3]*ai[3];
	ao[1] = Ai[1][0]*ai[0] + Ai[1][1]*ai[1] + Ai[1][2]*ai[2] + Ai[1][3]*ai[3];
	ao[2] = Ai[2][0]*ai[0] + Ai[2][1]*ai[1] + Ai[2][2]*ai[2] + Ai[2][3]*ai[3];
	ao[3] = Ai[3][0]*ai[0] + Ai[3][1]*ai[1] + Ai[3][2]*ai[2] + Ai[3][3]*ai[3];

	ao[4] = 0.5*(ao[0] + ao[1]);
	ao[5] = 0.5*(ao[1] + ao[2]);
	ao[6] = 0.5*(ao[2] + ao[0]);
	ao[7] = 0.5*(ao[0] + ao[3]);
	ao[8] = 0.5*(ao[1] + ao[3]);
	ao[9] = 0.5*(ao[2] + ao[3]);
}

//=============================================================================
//                          T E T 1 0 G 8
//=============================================================================
// I think this assumes that the tetrahedron in natural space is essentially
// a constant metric tet with the edge nodes at the center of the edges.
FETet10G8::FETet10G8() : FETet10_(NINT, FE_TET10G8)
{
	const double w = 1.0/6.0;
    gr[0] = 0.0158359099; gs[0] = 0.3280546970; gt[0] = 0.3280546970; gw[0] = 0.138527967*w;
    gr[1] = 0.3280546970; gs[1] = 0.0158359099; gt[1] = 0.3280546970; gw[1] = 0.138527967*w;
    gr[2] = 0.3280546970; gs[2] = 0.3280546970; gt[2] = 0.0158359099; gw[2] = 0.138527967*w;
    gr[3] = 0.3280546970; gs[3] = 0.3280546970; gt[3] = 0.3280546970; gw[3] = 0.138527967*w;
    gr[4] = 0.6791431780; gs[4] = 0.1069522740; gt[4] = 0.1069522740; gw[4] = 0.111472033*w;
    gr[5] = 0.1069522740; gs[5] = 0.6791431780; gt[5] = 0.1069522740; gw[5] = 0.111472033*w;
    gr[6] = 0.1069522740; gs[6] = 0.1069522740; gt[6] = 0.6791431780; gw[6] = 0.111472033*w;
    gr[7] = 0.1069522740; gs[7] = 0.1069522740; gt[7] = 0.1069522740; gw[7] = 0.111472033*w;

	init();

	// setup the shape function matrix
	N.resize(8, 4);
	for (int i=0; i<8; ++i)
	{
		N[i][0] = 1.0 - gr[i] - gs[i] - gt[i];
		N[i][1] = gr[i];
		N[i][2] = gs[i];
		N[i][3] = gt[i];
	}

	matrix A(4, 4);
	A = N.transpose()*N;

	// calculate inverse matrix
	Ai.resize(4, 4);
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETet10G8::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(4);
	for (int i=0; i<4; ++i)
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += N[j][i]*ai[j];
	}

	for (int i=0; i<4; ++i)
	{
		ao[i] = 0.0;
		for (int j=0; j<4; ++j) ao[i] += Ai[i][j]*b[j];
	}

	ao[4] = 0.5*(ao[0] + ao[1]);
	ao[5] = 0.5*(ao[1] + ao[2]);
	ao[6] = 0.5*(ao[2] + ao[0]);
	ao[7] = 0.5*(ao[0] + ao[3]);
	ao[8] = 0.5*(ao[1] + ao[3]);
	ao[9] = 0.5*(ao[2] + ao[3]);
}

//=============================================================================
//                          T E T 1 0 G 4 R I 1
//=============================================================================

FETet10G4RI1::FETet10G4RI1()
{
	m_pTRI = new FETet10G1;
}

//=============================================================================
//                          T E T 1 0 G 8 R I 4
//=============================================================================

FETet10G8RI4::FETet10G8RI4()
{
	m_pTRI = new FETet10G4;
}

//=============================================================================
//                             T E T 1 0 G L 1 1
//=============================================================================
// I think this assumes that the tetrahedron in natural space is essentially
// a constant metric tet with the edge nodes at the center of the edges.
FETet10GL11::FETet10GL11() : FETet10_(NINT, FE_TET10GL11)
{
	const double w = 1.0/6.0;
	const double a = w*1.0/60.0;
	const double b = w*4.0/60.0;
	gr[ 0] = 0.0; gs[ 0] = 0.0; gt[ 0] = 0.0; gw[ 0] = a;
	gr[ 1] = 1.0; gs[ 1] = 0.0; gt[ 1] = 0.0; gw[ 1] = a;
	gr[ 2] = 0.0; gs[ 2] = 1.0; gt[ 2] = 0.0; gw[ 2] = a;
	gr[ 3] = 0.0; gs[ 3] = 0.0; gt[ 3] = 1.0; gw[ 3] = a;
	gr[ 4] = 0.5; gs[ 4] = 0.0; gt[ 4] = 0.0; gw[ 4] = b;
	gr[ 5] = 0.5; gs[ 5] = 0.5; gt[ 5] = 0.0; gw[ 5] = b;
	gr[ 6] = 0.0; gs[ 6] = 0.5; gt[ 6] = 0.0; gw[ 6] = b;
	gr[ 7] = 0.0; gs[ 7] = 0.0; gt[ 7] = 0.5; gw[ 7] = b;
	gr[ 8] = 0.5; gs[ 8] = 0.0; gt[ 8] = 0.5; gw[ 8] = b;
	gr[ 9] = 0.0; gs[ 9] = 0.5; gt[ 9] = 0.5; gw[ 9] = b;
	gr[10] = 0.25; gs[10] = 0.25; gt[10] = 0.25; gw[10] = 32*a;
	init();
}

//-----------------------------------------------------------------------------
void FETet10GL11::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0]; ao[1] = ai[1]; ao[2] = ai[2]; ao[3] = ai[3];
	ao[4] = ai[4]; ao[5] = ai[5]; ao[6] = ai[6]; ao[7] = ai[7]; ao[8] = ai[8]; ao[9] = ai[9];
}

//=============================================================================
//                           T E T 1 5
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FETet15_::shape_fnc(double* H, double r, double s, double t)
{
	double r1 = 1.0 - r - s - t;
	double r2 = r;
	double r3 = s;
	double r4 = t;

	H[14] = 256*r1*r2*r3*r4;

	H[10] = 27.0*r1*r2*r3;
	H[11] = 27.0*r1*r2*r4;
	H[12] = 27.0*r2*r3*r4;
	H[13] = 27.0*r3*r1*r4;

	H[0] = r1*(2.0*r1 - 1.0) + (H[10] + H[11] + H[13])/9.0 - H[14]/64.0;
	H[1] = r2*(2.0*r2 - 1.0) + (H[10] + H[11] + H[12])/9.0 - H[14]/64.0;
	H[2] = r3*(2.0*r3 - 1.0) + (H[10] + H[12] + H[13])/9.0 - H[14]/64.0;
	H[3] = r4*(2.0*r4 - 1.0) + (H[11] + H[12] + H[13])/9.0 - H[14]/64.0;

	H[4] = 4.0*r1*r2 - 4.0*(H[10] + H[11])/9.0 + H[14]/8.0;
	H[5] = 4.0*r2*r3 - 4.0*(H[10] + H[12])/9.0 + H[14]/8.0;
	H[6] = 4.0*r3*r1 - 4.0*(H[10] + H[13])/9.0 + H[14]/8.0;
	H[7] = 4.0*r1*r4 - 4.0*(H[11] + H[13])/9.0 + H[14]/8.0;
	H[8] = 4.0*r2*r4 - 4.0*(H[11] + H[12])/9.0 + H[14]/8.0;
	H[9] = 4.0*r3*r4 - 4.0*(H[12] + H[13])/9.0 + H[14]/8.0;

	H[10] -= 27.0*H[14]/64.0;
	H[11] -= 27.0*H[14]/64.0;
	H[12] -= 27.0*H[14]/64.0;
	H[13] -= 27.0*H[14]/64.0;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FETet15_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[14] = 256.0*s*t*(1.0 - 2.0*r - s - t);
	Hs[14] = 256.0*r*t*(1.0 - r - 2.0*s - t);
	Ht[14] = 256.0*r*s*(1.0 - r - s - 2.0*t);

	Hr[10] =  27.0*s*(1.0 - 2.0*r - s - t);
	Hr[11] =  27.0*t*(1.0 - 2.0*r - s - t);
	Hr[12] =  27.0*s*t;
	Hr[13] = -27.0*s*t;

	Hs[10] =  27.0*r*(1.0 - r - 2.0*s - t);
	Hs[11] = -27.0*r*t;
	Hs[12] =  27.0*r*t;
	Hs[13] =  27.0*t*(1.0 - r - 2.0*s - t);

	Ht[10] = -27.0*r*s;
	Ht[11] =  27.0*r*(1.0 - r - s - 2.0*t);
	Ht[12] =  27.0*r*s;
	Ht[13] =  27.0*s*(1.0 - r - s - 2.0*t);

	Hr[0] = -3.0 + 4.0*r + 4.0*(s + t) + (Hr[10] + Hr[11] + Hr[13])/9.0 - Hr[14]/64.0;
	Hr[1] =  4.0*r - 1.0			   + (Hr[10] + Hr[11] + Hr[12])/9.0 - Hr[14]/64.0;
	Hr[2] =  0.0					   + (Hr[10] + Hr[12] + Hr[13])/9.0 - Hr[14]/64.0;
	Hr[3] =  0.0					   + (Hr[11] + Hr[12] + Hr[13])/9.0 - Hr[14]/64.0;
	Hr[4] =  4.0 - 8.0*r - 4.0*(s + t) - 4.0*(Hr[10] + Hr[11])/9.0 + Hr[14]/8.0;
	Hr[5] =  4.0*s					   - 4.0*(Hr[10] + Hr[12])/9.0 + Hr[14]/8.0;
	Hr[6] = -4.0*s					   - 4.0*(Hr[10] + Hr[13])/9.0 + Hr[14]/8.0;
	Hr[7] = -4.0*t					   - 4.0*(Hr[11] + Hr[13])/9.0 + Hr[14]/8.0;
	Hr[8] =  4.0*t					   - 4.0*(Hr[11] + Hr[12])/9.0 + Hr[14]/8.0;
	Hr[9] =  0.0					   - 4.0*(Hr[12] + Hr[13])/9.0 + Hr[14]/8.0;

	Hs[0] = -3.0 + 4.0*s + 4.0*(r + t) + (Hs[10] + Hs[11] + Hs[13])/9.0 - Hs[14]/64.0;
	Hs[1] =  0.0					   + (Hs[10] + Hs[11] + Hs[12])/9.0 - Hs[14]/64.0;
	Hs[2] =  4.0*s - 1.0			   + (Hs[10] + Hs[12] + Hs[13])/9.0 - Hs[14]/64.0;
	Hs[3] =  0.0					   + (Hs[11] + Hs[12] + Hs[13])/9.0 - Hs[14]/64.0;
	Hs[4] = -4.0*r					   - 4.0*(Hs[10] + Hs[11])/9.0 + Hs[14]/8.0;
	Hs[5] =  4.0*r					   - 4.0*(Hs[10] + Hs[12])/9.0 + Hs[14]/8.0;
	Hs[6] =  4.0 - 8.0*s - 4.0*(r + t) - 4.0*(Hs[10] + Hs[13])/9.0 + Hs[14]/8.0;
	Hs[7] = -4.0*t					   - 4.0*(Hs[11] + Hs[13])/9.0 + Hs[14]/8.0;
	Hs[8] =  0.0					   - 4.0*(Hs[11] + Hs[12])/9.0 + Hs[14]/8.0;
	Hs[9] =  4.0*t					   - 4.0*(Hs[12] + Hs[13])/9.0 + Hs[14]/8.0;

	Ht[0] = -3.0 + 4.0*t + 4.0*(r + s) + (Ht[10] + Ht[11] + Ht[13])/9.0 - Ht[14]/64.0;
	Ht[1] =  0.0					   + (Ht[10] + Ht[11] + Ht[12])/9.0 - Ht[14]/64.0;
	Ht[2] =  0.0					   + (Ht[10] + Ht[12] + Ht[13])/9.0 - Ht[14]/64.0;
	Ht[3] =  4.0*t - 1.0			   + (Ht[11] + Ht[12] + Ht[13])/9.0 - Ht[14]/64.0;
	Ht[4] = -4.0*r					   - 4.0*(Ht[10] + Ht[11])/9.0 + Ht[14]/8.0;
	Ht[5] =  0.0					   - 4.0*(Ht[10] + Ht[12])/9.0 + Ht[14]/8.0;
	Ht[6] = -4.0*s					   - 4.0*(Ht[10] + Ht[13])/9.0 + Ht[14]/8.0;
	Ht[7] =  4.0 - 8.0*t - 4.0*(r + s) - 4.0*(Ht[11] + Ht[13])/9.0 + Ht[14]/8.0;
	Ht[8] =  4.0*r					   - 4.0*(Ht[11] + Ht[12])/9.0 + Ht[14]/8.0;
	Ht[9] =  4.0*s					   - 4.0*(Ht[12] + Ht[13])/9.0 + Ht[14]/8.0;

	Hr[10] -= 27.0*Hr[14]/64.0;
	Hr[11] -= 27.0*Hr[14]/64.0;
	Hr[12] -= 27.0*Hr[14]/64.0;
	Hr[13] -= 27.0*Hr[14]/64.0;

	Hs[10] -= 27.0*Hs[14]/64.0;
	Hs[11] -= 27.0*Hs[14]/64.0;
	Hs[12] -= 27.0*Hs[14]/64.0;
	Hs[13] -= 27.0*Hs[14]/64.0;

	Ht[10] -= 27.0*Ht[14]/64.0;
	Ht[11] -= 27.0*Ht[14]/64.0;
	Ht[12] -= 27.0*Ht[14]/64.0;
	Ht[13] -= 27.0*Ht[14]/64.0;
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
//! \todo implement this
void FETet15_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{

}

//=============================================================================
//                          T E T 1 5 G 4
//=============================================================================
FETet15G4::FETet15G4() : FETet15_(NINT, FE_TET15G4)
{
	// integration point coordinates
	const double a = 0.58541020;
	const double b = 0.13819660;
	const double w = 0.25 / 6.0;
	gr[ 0] = a; gs[ 0] = b; gt[ 0] = b; gw[ 0] = w;
	gr[ 1] = b; gs[ 1] = a; gt[ 1] = b; gw[ 1] = w;
	gr[ 2] = b; gs[ 2] = b; gt[ 2] = a; gw[ 2] = w;
	gr[ 3] = b; gs[ 3] = b; gt[ 3] = b; gw[ 3] = w;

	init();

	// setup the shape function matrix
	matrix A(4,4);
	for (int i=0; i<4; ++i)
	{
		double r = gr[i];
		double s = gs[i];
		double t = gt[i];

		A[i][0] = 1 - r - s - t;
		A[i][1] = r;
		A[i][2] = s;
		A[i][3] = t;
	}

	// calculate inverse matrix
	Ai.resize(4, 4);
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETet15G4::project_to_nodes(double* ai, double* ao)
{
	ao[0] = Ai[0][0]*ai[0] + Ai[0][1]*ai[1] + Ai[0][2]*ai[2] + Ai[0][3]*ai[3];
	ao[1] = Ai[1][0]*ai[0] + Ai[1][1]*ai[1] + Ai[1][2]*ai[2] + Ai[1][3]*ai[3];
	ao[2] = Ai[2][0]*ai[0] + Ai[2][1]*ai[1] + Ai[2][2]*ai[2] + Ai[2][3]*ai[3];
	ao[3] = Ai[3][0]*ai[0] + Ai[3][1]*ai[1] + Ai[3][2]*ai[2] + Ai[3][3]*ai[3];

	ao[4] = 0.5*(ao[0] + ao[1]);
	ao[5] = 0.5*(ao[1] + ao[2]);
	ao[6] = 0.5*(ao[2] + ao[0]);
	ao[7] = 0.5*(ao[0] + ao[3]);
	ao[8] = 0.5*(ao[1] + ao[3]);
	ao[9] = 0.5*(ao[2] + ao[3]);

	ao[10] = (ao[0] + ao[1] + ao[2])/3.0;
	ao[11] = (ao[0] + ao[1] + ao[3])/3.0;
	ao[12] = (ao[1] + ao[2] + ao[3])/3.0;
	ao[13] = (ao[0] + ao[2] + ao[3])/3.0;

	ao[14] = 0.25*(ao[0] + ao[1] + ao[2] + ao[3]);
}


//=============================================================================
//                          T E T 1 5 G 8
//=============================================================================
// I think this assumes that the tetrahedron in natural space is essentially
// a constant metric tet with the edge nodes at the center of the edges.
FETet15G8::FETet15G8() : FETet15_(NINT, FE_TET15G8)
{
	const double w = 1.0/6.0;
    gr[0] = 0.0158359099; gs[0] = 0.3280546970; gt[0] = 0.3280546970; gw[0] = 0.138527967*w;
    gr[1] = 0.3280546970; gs[1] = 0.0158359099; gt[1] = 0.3280546970; gw[1] = 0.138527967*w;
    gr[2] = 0.3280546970; gs[2] = 0.3280546970; gt[2] = 0.0158359099; gw[2] = 0.138527967*w;
    gr[3] = 0.3280546970; gs[3] = 0.3280546970; gt[3] = 0.3280546970; gw[3] = 0.138527967*w;
    gr[4] = 0.6791431780; gs[4] = 0.1069522740; gt[4] = 0.1069522740; gw[4] = 0.111472033*w;
    gr[5] = 0.1069522740; gs[5] = 0.6791431780; gt[5] = 0.1069522740; gw[5] = 0.111472033*w;
    gr[6] = 0.1069522740; gs[6] = 0.1069522740; gt[6] = 0.6791431780; gw[6] = 0.111472033*w;
    gr[7] = 0.1069522740; gs[7] = 0.1069522740; gt[7] = 0.1069522740; gw[7] = 0.111472033*w;

	init();

	// setup the shape function matrix
	N.resize(8, 4);
	for (int i=0; i<8; ++i)
	{
		N[i][0] = 1.0 - gr[i] - gs[i] - gt[i];
		N[i][1] = gr[i];
		N[i][2] = gs[i];
		N[i][3] = gt[i];
	}

	matrix A(4, 4);
	A = N.transpose()*N;

	// calculate inverse matrix
	Ai.resize(4, 4);
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETet15G8::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(4);
	for (int i=0; i<4; ++i)
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += N[j][i]*ai[j];
	}

	for (int i=0; i<4; ++i)
	{
		ao[i] = 0.0;
		for (int j=0; j<4; ++j) ao[i] += Ai[i][j]*b[j];
	}

	ao[4] = 0.5*(ao[0] + ao[1]);
	ao[5] = 0.5*(ao[1] + ao[2]);
	ao[6] = 0.5*(ao[2] + ao[0]);
	ao[7] = 0.5*(ao[0] + ao[3]);
	ao[8] = 0.5*(ao[1] + ao[3]);
	ao[9] = 0.5*(ao[2] + ao[3]);

	ao[10] = (ao[0] + ao[1] + ao[2])/3.0;
	ao[11] = (ao[0] + ao[1] + ao[3])/3.0;
	ao[12] = (ao[1] + ao[2] + ao[3])/3.0;
	ao[13] = (ao[0] + ao[2] + ao[3])/3.0;

	ao[14] = (ao[0] + ao[1] + ao[2] + ao[3])*0.25;
}

//=============================================================================
//                          T E T 1 5 G 1 1
//=============================================================================
FETet15G11::FETet15G11() : FETet15_(NINT, FE_TET15G11)
{
    gr[0] = 0.25; gs[0] = 0.25; gt[0] = 0.25; gw[0] = -0.01315555556;

    gr[1] = 0.071428571428571; gs[1] = 0.071428571428571; gt[1] = 0.071428571428571; gw[1] = 0.007622222222;
    gr[2] = 0.785714285714286; gs[2] = 0.071428571428571; gt[2] = 0.071428571428571; gw[2] = 0.007622222222;
    gr[3] = 0.071428571428571; gs[3] = 0.785714285714286; gt[3] = 0.071428571428571; gw[3] = 0.007622222222;
    gr[4] = 0.071428571428571; gs[4] = 0.071428571428571; gt[4] = 0.785714285714286; gw[4] = 0.007622222222;

    gr[ 5] = 0.399403576166799; gs[ 5] = 0.100596423833201; gt[ 5] = 0.100596423833201; gw[ 5] = 0.024888888889;
    gr[ 6] = 0.100596423833201; gs[ 6] = 0.399403576166799; gt[ 6] = 0.100596423833201; gw[ 6] = 0.024888888889;
    gr[ 7] = 0.100596423833201; gs[ 7] = 0.100596423833201; gt[ 7] = 0.399403576166799; gw[ 7] = 0.024888888889;
    gr[ 8] = 0.399403576166799; gs[ 8] = 0.399403576166799; gt[ 8] = 0.100596423833201; gw[ 8] = 0.024888888889;
    gr[ 9] = 0.399403576166799; gs[ 9] = 0.100596423833201; gt[ 9] = 0.399403576166799; gw[ 9] = 0.024888888889;
    gr[10] = 0.100596423833201; gs[10] = 0.399403576166799; gt[10] = 0.399403576166799; gw[10] = 0.024888888889;

	init();

	// setup the shape function matrix
	N.resize(11, 4);
	for (int i=0; i<11; ++i)
	{
		N[i][0] = 1.0 - gr[i] - gs[i] - gt[i];
		N[i][1] = gr[i];
		N[i][2] = gs[i];
		N[i][3] = gt[i];
	}

	matrix A(4, 4);
	A = N.transpose()*N;

	// calculate inverse matrix
	Ai.resize(4, 4);
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETet15G11::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(4);
	for (int i=0; i<4; ++i)
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += N[j][i]*ai[j];
	}

	for (int i=0; i<4; ++i)
	{
		ao[i] = 0.0;
		for (int j=0; j<4; ++j) ao[i] += Ai[i][j]*b[j];
	}

	ao[4] = 0.5*(ao[0] + ao[1]);
	ao[5] = 0.5*(ao[1] + ao[2]);
	ao[6] = 0.5*(ao[2] + ao[0]);
	ao[7] = 0.5*(ao[0] + ao[3]);
	ao[8] = 0.5*(ao[1] + ao[3]);
	ao[9] = 0.5*(ao[2] + ao[3]);

	ao[10] = (ao[0] + ao[1] + ao[2])/3.0;
	ao[11] = (ao[0] + ao[1] + ao[3])/3.0;
	ao[12] = (ao[1] + ao[2] + ao[3])/3.0;
	ao[13] = (ao[0] + ao[2] + ao[3])/3.0;

	ao[14] = (ao[0] + ao[1] + ao[2] + ao[3])*0.25;
}

//=============================================================================
//                          T E T 1 5 G 1 5
//=============================================================================
FETet15G15::FETet15G15() : FETet15_(NINT, FE_TET15G15)
{
    gr[0] = 0.25; gs[0] = 0.25; gt[0] = 0.25; gw[0] = 0.030283678097089;

    gr[1] = 0.333333333333333; gs[1] = 0.333333333333333; gt[1] = 0.333333333333333; gw[1] = 0.006026785714286;
    gr[2] = 0.000000000000000; gs[2] = 0.333333333333333; gt[2] = 0.333333333333333; gw[2] = 0.006026785714286;
    gr[3] = 0.333333333333333; gs[3] = 0.000000000000000; gt[3] = 0.333333333333333; gw[3] = 0.006026785714286;
    gr[4] = 0.333333333333333; gs[4] = 0.333333333333333; gt[4] = 0.000000000000000; gw[4] = 0.006026785714286;

    gr[ 5] = 0.090909090909091; gs[ 5] = 0.090909090909091; gt[ 5] = 0.090909090909091; gw[ 5] = 0.011645249086029;
    gr[ 6] = 0.727272727272727; gs[ 6] = 0.090909090909091; gt[ 6] = 0.090909090909091; gw[ 6] = 0.011645249086029;
    gr[ 7] = 0.090909090909091; gs[ 7] = 0.727272727272727; gt[ 7] = 0.090909090909091; gw[ 7] = 0.011645249086029;
    gr[ 8] = 0.090909090909091; gs[ 8] = 0.090909090909091; gt[ 8] = 0.727272727272727; gw[ 8] = 0.011645249086029;

    gr[ 9] = 0.433449846426336; gs[ 9] = 0.066550153573664; gt[ 9] = 0.066550153573664; gw[ 9] = 0.010949141561386;
    gr[10] = 0.066550153573664; gs[10] = 0.433449846426336; gt[10] = 0.066550153573664; gw[10] = 0.010949141561386;
    gr[11] = 0.066550153573664; gs[11] = 0.066550153573664; gt[11] = 0.433449846426336; gw[11] = 0.010949141561386;
    gr[12] = 0.066550153573664; gs[12] = 0.433449846426336; gt[12] = 0.433449846426336; gw[12] = 0.010949141561386;
    gr[13] = 0.433449846426336; gs[13] = 0.066550153573664; gt[13] = 0.433449846426336; gw[13] = 0.010949141561386;
    gr[14] = 0.433449846426336; gs[14] = 0.433449846426336; gt[14] = 0.066550153573664; gw[14] = 0.010949141561386;

	init();

	// setup the shape function matrix
	N.resize(NINT, 4);
	for (int i=0; i<NINT; ++i)
	{
		N[i][0] = 1.0 - gr[i] - gs[i] - gt[i];
		N[i][1] = gr[i];
		N[i][2] = gs[i];
		N[i][3] = gt[i];
	}

	matrix A(4, 4);
	A = N.transpose()*N;

	// calculate inverse matrix
	Ai.resize(4, 4);
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETet15G15::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(4);
	for (int i=0; i<4; ++i)
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += N[j][i]*ai[j];
	}

	for (int i=0; i<4; ++i)
	{
		ao[i] = 0.0;
		for (int j=0; j<4; ++j) ao[i] += Ai[i][j]*b[j];
	}

	ao[4] = 0.5*(ao[0] + ao[1]);
	ao[5] = 0.5*(ao[1] + ao[2]);
	ao[6] = 0.5*(ao[2] + ao[0]);
	ao[7] = 0.5*(ao[0] + ao[3]);
	ao[8] = 0.5*(ao[1] + ao[3]);
	ao[9] = 0.5*(ao[2] + ao[3]);

	ao[10] = (ao[0] + ao[1] + ao[2])/3.0;
	ao[11] = (ao[0] + ao[1] + ao[3])/3.0;
	ao[12] = (ao[1] + ao[2] + ao[3])/3.0;
	ao[13] = (ao[0] + ao[2] + ao[3])/3.0;

	ao[14] = (ao[0] + ao[1] + ao[2] + ao[3])*0.25;
}

//=============================================================================
//                          T E T 1 5 G 1 5 R I 4
//=============================================================================

FETet15G15RI4::FETet15G15RI4()
{
	m_pTRI = new FETet15G4;
}

//=============================================================================
//              H E X 2 0
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FEHex20_::shape_fnc(double* H, double r, double s, double t)
{
	H[ 8] = 0.25*(1 - r*r)*(1 - s)*(1 - t);
	H[ 9] = 0.25*(1 - s*s)*(1 + r)*(1 - t);
	H[10] = 0.25*(1 - r*r)*(1 + s)*(1 - t);
	H[11] = 0.25*(1 - s*s)*(1 - r)*(1 - t);
	H[12] = 0.25*(1 - r*r)*(1 - s)*(1 + t);
	H[13] = 0.25*(1 - s*s)*(1 + r)*(1 + t);
	H[14] = 0.25*(1 - r*r)*(1 + s)*(1 + t);
	H[15] = 0.25*(1 - s*s)*(1 - r)*(1 + t);
	H[16] = 0.25*(1 - t*t)*(1 - r)*(1 - s);
	H[17] = 0.25*(1 - t*t)*(1 + r)*(1 - s);
	H[18] = 0.25*(1 - t*t)*(1 + r)*(1 + s);
	H[19] = 0.25*(1 - t*t)*(1 - r)*(1 + s);

	H[0] = 0.125*(1 - r)*(1 - s)*(1 - t) - 0.5*(H[ 8] + H[11] + H[16]);
	H[1] = 0.125*(1 + r)*(1 - s)*(1 - t) - 0.5*(H[ 8] + H[ 9] + H[17]);
	H[2] = 0.125*(1 + r)*(1 + s)*(1 - t) - 0.5*(H[ 9] + H[10] + H[18]);
	H[3] = 0.125*(1 - r)*(1 + s)*(1 - t) - 0.5*(H[10] + H[11] + H[19]);
	H[4] = 0.125*(1 - r)*(1 - s)*(1 + t) - 0.5*(H[12] + H[15] + H[16]);
	H[5] = 0.125*(1 + r)*(1 - s)*(1 + t) - 0.5*(H[12] + H[13] + H[17]);
	H[6] = 0.125*(1 + r)*(1 + s)*(1 + t) - 0.5*(H[13] + H[14] + H[18]);
	H[7] = 0.125*(1 - r)*(1 + s)*(1 + t) - 0.5*(H[14] + H[15] + H[19]);
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FEHex20_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	Hr[ 8] = -0.5*r*(1 - s)*(1 - t);
	Hr[ 9] =  0.25*(1 - s*s)*(1 - t);
	Hr[10] = -0.5*r*(1 + s)*(1 - t);
	Hr[11] = -0.25*(1 - s*s)*(1 - t);
	Hr[12] = -0.5*r*(1 - s)*(1 + t);
	Hr[13] =  0.25*(1 - s*s)*(1 + t);
	Hr[14] = -0.5*r*(1 + s)*(1 + t);
	Hr[15] = -0.25*(1 - s*s)*(1 + t);
	Hr[16] = -0.25*(1 - t*t)*(1 - s);
	Hr[17] =  0.25*(1 - t*t)*(1 - s);
	Hr[18] =  0.25*(1 - t*t)*(1 + s);
	Hr[19] = -0.25*(1 - t*t)*(1 + s);

	Hr[0] = -0.125*(1 - s)*(1 - t) - 0.5*(Hr[ 8] + Hr[11] + Hr[16]);
	Hr[1] =  0.125*(1 - s)*(1 - t) - 0.5*(Hr[ 8] + Hr[ 9] + Hr[17]);
	Hr[2] =  0.125*(1 + s)*(1 - t) - 0.5*(Hr[ 9] + Hr[10] + Hr[18]);
	Hr[3] = -0.125*(1 + s)*(1 - t) - 0.5*(Hr[10] + Hr[11] + Hr[19]);
	Hr[4] = -0.125*(1 - s)*(1 + t) - 0.5*(Hr[12] + Hr[15] + Hr[16]);
	Hr[5] =  0.125*(1 - s)*(1 + t) - 0.5*(Hr[12] + Hr[13] + Hr[17]);
	Hr[6] =  0.125*(1 + s)*(1 + t) - 0.5*(Hr[13] + Hr[14] + Hr[18]);
	Hr[7] = -0.125*(1 + s)*(1 + t) - 0.5*(Hr[14] + Hr[15] + Hr[19]);
		
	Hs[ 8] = -0.25*(1 - r*r)*(1 - t);
	Hs[ 9] = -0.5*s*(1 + r)*(1 - t);
	Hs[10] = 0.25*(1 - r*r)*(1 - t);
	Hs[11] = -0.5*s*(1 - r)*(1 - t);
	Hs[12] = -0.25*(1 - r*r)*(1 + t);
	Hs[13] = -0.5*s*(1 + r)*(1 + t);
	Hs[14] = 0.25*(1 - r*r)*(1 + t);
	Hs[15] = -0.5*s*(1 - r)*(1 + t);
	Hs[16] = -0.25*(1 - t*t)*(1 - r);
	Hs[17] = -0.25*(1 - t*t)*(1 + r);
	Hs[18] =  0.25*(1 - t*t)*(1 + r);
	Hs[19] =  0.25*(1 - t*t)*(1 - r);

	Hs[0] = -0.125*(1 - r)*(1 - t) - 0.5*(Hs[ 8] + Hs[11] + Hs[16]);
	Hs[1] = -0.125*(1 + r)*(1 - t) - 0.5*(Hs[ 8] + Hs[ 9] + Hs[17]);
	Hs[2] =  0.125*(1 + r)*(1 - t) - 0.5*(Hs[ 9] + Hs[10] + Hs[18]);
	Hs[3] =  0.125*(1 - r)*(1 - t) - 0.5*(Hs[10] + Hs[11] + Hs[19]);
	Hs[4] = -0.125*(1 - r)*(1 + t) - 0.5*(Hs[12] + Hs[15] + Hs[16]);
	Hs[5] = -0.125*(1 + r)*(1 + t) - 0.5*(Hs[12] + Hs[13] + Hs[17]);
	Hs[6] =  0.125*(1 + r)*(1 + t) - 0.5*(Hs[13] + Hs[14] + Hs[18]);
	Hs[7] =  0.125*(1 - r)*(1 + t) - 0.5*(Hs[14] + Hs[15] + Hs[19]);

	Ht[ 8] = -0.25*(1 - r*r)*(1 - s);
	Ht[ 9] = -0.25*(1 - s*s)*(1 + r);
	Ht[10] = -0.25*(1 - r*r)*(1 + s);
	Ht[11] = -0.25*(1 - s*s)*(1 - r);
	Ht[12] =  0.25*(1 - r*r)*(1 - s);
	Ht[13] =  0.25*(1 - s*s)*(1 + r);
	Ht[14] =  0.25*(1 - r*r)*(1 + s);
	Ht[15] =  0.25*(1 - s*s)*(1 - r);
	Ht[16] = -0.5*t*(1 - r)*(1 - s);
	Ht[17] = -0.5*t*(1 + r)*(1 - s);
	Ht[18] = -0.5*t*(1 + r)*(1 + s);
	Ht[19] = -0.5*t*(1 - r)*(1 + s);
		
	Ht[0] = -0.125*(1 - r)*(1 - s) - 0.5*(Ht[ 8] + Ht[11] + Ht[16]);
	Ht[1] = -0.125*(1 + r)*(1 - s) - 0.5*(Ht[ 8] + Ht[ 9] + Ht[17]);
	Ht[2] = -0.125*(1 + r)*(1 + s) - 0.5*(Ht[ 9] + Ht[10] + Ht[18]);
	Ht[3] = -0.125*(1 - r)*(1 + s) - 0.5*(Ht[10] + Ht[11] + Ht[19]);
	Ht[4] =  0.125*(1 - r)*(1 - s) - 0.5*(Ht[12] + Ht[15] + Ht[16]);
	Ht[5] =  0.125*(1 + r)*(1 - s) - 0.5*(Ht[12] + Ht[13] + Ht[17]);
	Ht[6] =  0.125*(1 + r)*(1 + s) - 0.5*(Ht[13] + Ht[14] + Ht[18]);
	Ht[7] =  0.125*(1 - r)*(1 + s) - 0.5*(Ht[14] + Ht[15] + Ht[19]);
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
//! \todo implement this (needed for biphasic problems)
void FEHex20_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	// Hrr
	Hrr[ 8] = -0.5*(1 - s)*(1 - t);
	Hrr[ 9] =  0.;
	Hrr[10] = -0.5*(1 + s)*(1 - t);
	Hrr[11] =  0.;
	Hrr[12] = -0.5*(1 - s)*(1 + t);
	Hrr[13] =  0.;
	Hrr[14] = -0.5*(1 + s)*(1 + t);
	Hrr[15] =  0.;
	Hrr[16] =  0.;
	Hrr[17] =  0.;
	Hrr[18] =  0.;
	Hrr[19] =  0.;

	Hrr[0] = - 0.5*(Hrr[ 8] + Hrr[11] + Hrr[16]);
	Hrr[1] = - 0.5*(Hrr[ 8] + Hrr[ 9] + Hrr[17]);
	Hrr[2] = - 0.5*(Hrr[ 9] + Hrr[10] + Hrr[18]);
	Hrr[3] = - 0.5*(Hrr[10] + Hrr[11] + Hrr[19]);
	Hrr[4] = - 0.5*(Hrr[12] + Hrr[15] + Hrr[16]);
	Hrr[5] = - 0.5*(Hrr[12] + Hrr[13] + Hrr[17]);
	Hrr[6] = - 0.5*(Hrr[13] + Hrr[14] + Hrr[18]);
	Hrr[7] = - 0.5*(Hrr[14] + Hrr[15] + Hrr[19]);

	// Hss
	Hss[ 8] =  0.;
	Hss[ 9] = -0.5*(1 + r)*(1 - t);
	Hss[10] =  0.;
	Hss[11] = -0.5*(1 - r)*(1 - t);
	Hss[12] =  0.;
	Hss[13] = -0.5*(1 + r)*(1 + t);
	Hss[14] =  0.;
	Hss[15] = -0.5*(1 - r)*(1 + t);
	Hss[16] =  0.;
	Hss[17] =  0.;
	Hss[18] =  0.;
	Hss[19] =  0.;

	Hss[0] = - 0.5*(Hss[ 8] + Hss[11] + Hss[16]);
	Hss[1] = - 0.5*(Hss[ 8] + Hss[ 9] + Hss[17]);
	Hss[2] = - 0.5*(Hss[ 9] + Hss[10] + Hss[18]);
	Hss[3] = - 0.5*(Hss[10] + Hss[11] + Hss[19]);
	Hss[4] = - 0.5*(Hss[12] + Hss[15] + Hss[16]);
	Hss[5] = - 0.5*(Hss[12] + Hss[13] + Hss[17]);
	Hss[6] = - 0.5*(Hss[13] + Hss[14] + Hss[18]);
	Hss[7] = - 0.5*(Hss[14] + Hss[15] + Hss[19]);

	// Htt
	Htt[ 8] =  0.;
	Htt[ 9] =  0.;
	Htt[10] =  0.;
	Htt[11] =  0.;
	Htt[12] =  0.;
	Htt[13] =  0.;
	Htt[14] =  0.;
	Htt[15] =  0.;
	Htt[16] = -0.5*(1 - r)*(1 - s);
	Htt[17] = -0.5*(1 + r)*(1 - s);
	Htt[18] = -0.5*(1 + r)*(1 + s);
	Htt[19] = -0.5*(1 - r)*(1 + s);

	Htt[0] = - 0.5*(Htt[ 8] + Htt[11] + Htt[16]);
	Htt[1] = - 0.5*(Htt[ 8] + Htt[ 9] + Htt[17]);
	Htt[2] = - 0.5*(Htt[ 9] + Htt[10] + Htt[18]);
	Htt[3] = - 0.5*(Htt[10] + Htt[11] + Htt[19]);
	Htt[4] = - 0.5*(Htt[12] + Htt[15] + Htt[16]);
	Htt[5] = - 0.5*(Htt[12] + Htt[13] + Htt[17]);
	Htt[6] = - 0.5*(Htt[13] + Htt[14] + Htt[18]);
	Htt[7] = - 0.5*(Htt[14] + Htt[15] + Htt[19]);

	// Hrs
	Hrs[ 8] =  0.5*r*(1 - t);
	Hrs[ 9] = -0.5*s*(1 - t);
	Hrs[10] = -0.5*r*(1 - t);
	Hrs[11] =  0.5*s*(1 - t);
	Hrs[12] =  0.5*r*(1 + t);
	Hrs[13] = -0.5*s*(1 + t);
	Hrs[14] = -0.5*r*(1 + t);
	Hrs[15] =  0.5*s*(1 + t);
	Hrs[16] =  0.25*(1 - t*t);
	Hrs[17] = -0.25*(1 - t*t);
	Hrs[18] =  0.25*(1 - t*t);
	Hrs[19] = -0.25*(1 - t*t);

	Hrs[0] =  0.125*(1 - t) - 0.5*(Hrs[ 8] + Hrs[11] + Hrs[16]);
	Hrs[1] = -0.125*(1 - t) - 0.5*(Hrs[ 8] + Hrs[ 9] + Hrs[17]);
	Hrs[2] =  0.125*(1 - t) - 0.5*(Hrs[ 9] + Hrs[10] + Hrs[18]);
	Hrs[3] = -0.125*(1 - t) - 0.5*(Hrs[10] + Hrs[11] + Hrs[19]);
	Hrs[4] =  0.125*(1 + t) - 0.5*(Hrs[12] + Hrs[15] + Hrs[16]);
	Hrs[5] = -0.125*(1 + t) - 0.5*(Hrs[12] + Hrs[13] + Hrs[17]);
	Hrs[6] =  0.125*(1 + t) - 0.5*(Hrs[13] + Hrs[14] + Hrs[18]);
	Hrs[7] = -0.125*(1 + t) - 0.5*(Hrs[14] + Hrs[15] + Hrs[19]);

	// Hst
	Hst[ 8] =  0.25*(1 - r*r);
	Hst[ 9] =  0.5*s*(1 + r);
	Hst[10] = -0.25*(1 - r*r);
	Hst[11] =  0.5*s*(1 - r);
	Hst[12] = -0.25*(1 - r*r);
	Hst[13] = -0.5*s*(1 + r);
	Hst[14] =  0.25*(1 - r*r);
	Hst[15] = -0.5*s*(1 - r);
	Hst[16] =  0.5*t*(1 - r);
	Hst[17] =  0.5*t*(1 + r);
	Hst[18] = -0.5*t*(1 + r);
	Hst[19] = -0.5*t*(1 - r);

	Hst[0] =  0.125*(1 - r) - 0.5*(Hst[ 8] + Hst[11] + Hst[16]);
	Hst[1] =  0.125*(1 + r) - 0.5*(Hst[ 8] + Hst[ 9] + Hst[17]);
	Hst[2] = -0.125*(1 + r) - 0.5*(Hst[ 9] + Hst[10] + Hst[18]);
	Hst[3] = -0.125*(1 - r) - 0.5*(Hst[10] + Hst[11] + Hst[19]);
	Hst[4] = -0.125*(1 - r) - 0.5*(Hst[12] + Hst[15] + Hst[16]);
	Hst[5] = -0.125*(1 + r) - 0.5*(Hst[12] + Hst[13] + Hst[17]);
	Hst[6] =  0.125*(1 + r) - 0.5*(Hst[13] + Hst[14] + Hst[18]);
	Hst[7] =  0.125*(1 - r) - 0.5*(Hst[14] + Hst[15] + Hst[19]);

	// Hrt
	Hrt[ 8] =  0.5*r*(1 - s);
	Hrt[ 9] = -0.25*(1 - s*s);
	Hrt[10] =  0.5*r*(1 + s);
	Hrt[11] =  0.25*(1 - s*s);
	Hrt[12] = -0.5*r*(1 - s);
	Hrt[13] =  0.25*(1 - s*s);
	Hrt[14] = -0.5*r*(1 + s);
	Hrt[15] = -0.25*(1 - s*s);
	Hrt[16] =  0.5*t*(1 - s);
	Hrt[17] = -0.5*t*(1 - s);
	Hrt[18] = -0.5*t*(1 + s);
	Hrt[19] =  0.5*t*(1 + s);	
	
	Hrt[0] =  0.125*(1 - s) - 0.5*(Hrt[ 8] + Hrt[11] + Hrt[16]);
	Hrt[1] = -0.125*(1 - s) - 0.5*(Hrt[ 8] + Hrt[ 9] + Hrt[17]);
	Hrt[2] = -0.125*(1 + s) - 0.5*(Hrt[ 9] + Hrt[10] + Hrt[18]);
	Hrt[3] =  0.125*(1 + s) - 0.5*(Hrt[10] + Hrt[11] + Hrt[19]);
	Hrt[4] = -0.125*(1 - s) - 0.5*(Hrt[12] + Hrt[15] + Hrt[16]);
	Hrt[5] =  0.125*(1 - s) - 0.5*(Hrt[12] + Hrt[13] + Hrt[17]);
	Hrt[6] =  0.125*(1 + s) - 0.5*(Hrt[13] + Hrt[14] + Hrt[18]);
	Hrt[7] = -0.125*(1 + s) - 0.5*(Hrt[14] + Hrt[15] + Hrt[19]);


}

//=============================================================================
//              H E X 2 0 G 2 7
//=============================================================================

FEHex20G27::FEHex20G27() : FEHex20_(NINT, FE_HEX20G27)
{
	// integration point coordinates
	const double a = 0.774596669241483;
	const double w1 = 5.0 / 9.0;
	const double w2 = 8.0 / 9.0;
	gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -a; gw[ 0] = w1*w1*w1;
	gr[ 1] =  0; gs[ 1] = -a; gt[ 1] = -a; gw[ 1] = w2*w1*w1;
	gr[ 2] =  a; gs[ 2] = -a; gt[ 2] = -a; gw[ 2] = w1*w1*w1;
	gr[ 3] = -a; gs[ 3] =  0; gt[ 3] = -a; gw[ 3] = w1*w2*w1;
	gr[ 4] =  0; gs[ 4] =  0; gt[ 4] = -a; gw[ 4] = w2*w2*w1;
	gr[ 5] =  a; gs[ 5] =  0; gt[ 5] = -a; gw[ 5] = w1*w2*w1;
	gr[ 6] = -a; gs[ 6] =  a; gt[ 6] = -a; gw[ 6] = w1*w1*w1;
	gr[ 7] =  0; gs[ 7] =  a; gt[ 7] = -a; gw[ 7] = w2*w1*w1;
	gr[ 8] =  a; gs[ 8] =  a; gt[ 8] = -a; gw[ 8] = w1*w1*w1;
	gr[ 9] = -a; gs[ 9] = -a; gt[ 9] =  0; gw[ 9] = w1*w1*w2;
	gr[10] =  0; gs[10] = -a; gt[10] =  0; gw[10] = w2*w1*w2;
	gr[11] =  a; gs[11] = -a; gt[11] =  0; gw[11] = w1*w1*w2;
	gr[12] = -a; gs[12] =  0; gt[12] =  0; gw[12] = w1*w2*w2;
	gr[13] =  0; gs[13] =  0; gt[13] =  0; gw[13] = w2*w2*w2;
	gr[14] =  a; gs[14] =  0; gt[14] =  0; gw[14] = w1*w2*w2;
	gr[15] = -a; gs[15] =  a; gt[15] =  0; gw[15] = w1*w1*w2;
	gr[16] =  0; gs[16] =  a; gt[16] =  0; gw[16] = w2*w1*w2;
	gr[17] =  a; gs[17] =  a; gt[17] =  0; gw[17] = w1*w1*w2;
	gr[18] = -a; gs[18] = -a; gt[18] =  a; gw[18] = w1*w1*w1;
	gr[19] =  0; gs[19] = -a; gt[19] =  a; gw[19] = w2*w1*w1;
	gr[20] =  a; gs[20] = -a; gt[20] =  a; gw[20] = w1*w1*w1;
	gr[21] = -a; gs[21] =  0; gt[21] =  a; gw[21] = w1*w2*w1;
	gr[22] =  0; gs[22] =  0; gt[22] =  a; gw[22] = w2*w2*w1;
	gr[23] =  a; gs[23] =  0; gt[23] =  a; gw[23] = w1*w2*w1;
	gr[24] = -a; gs[24] =  a; gt[24] =  a; gw[24] = w1*w1*w1;
	gr[25] =  0; gs[25] =  a; gt[25] =  a; gw[25] = w2*w1*w1;
	gr[26] =  a; gs[26] =  a; gt[26] =  a; gw[26] = w1*w1*w1;

	init();
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FEHex20G27::project_to_nodes(double* ai, double* ao)
{
	
}

//=============================================================================
//              H E X 2 7
//=============================================================================

//-----------------------------------------------------------------------------
//! values of shape functions
void FEHex27_::shape_fnc(double* H, double r, double s, double t)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};
	double T[3] = {0.5*t*(t-1.0), 0.5*t*(t+1.0), 1.0 - t*t};

	H[ 0] = R[0]*S[0]*T[0];
	H[ 1] = R[1]*S[0]*T[0];
	H[ 2] = R[1]*S[1]*T[0];
	H[ 3] = R[0]*S[1]*T[0];
	H[ 4] = R[0]*S[0]*T[1];
	H[ 5] = R[1]*S[0]*T[1];
	H[ 6] = R[1]*S[1]*T[1];
	H[ 7] = R[0]*S[1]*T[1];
	H[ 8] = R[2]*S[0]*T[0];
	H[ 9] = R[1]*S[2]*T[0];
	H[10] = R[2]*S[1]*T[0];
	H[11] = R[0]*S[2]*T[0];
	H[12] = R[2]*S[0]*T[1];
	H[13] = R[1]*S[2]*T[1];
	H[14] = R[2]*S[1]*T[1];
	H[15] = R[0]*S[2]*T[1];
	H[16] = R[0]*S[0]*T[2];
	H[17] = R[1]*S[0]*T[2];
	H[18] = R[1]*S[1]*T[2];
	H[19] = R[0]*S[1]*T[2];
	H[20] = R[2]*S[0]*T[2];
	H[21] = R[1]*S[2]*T[2];
	H[22] = R[2]*S[1]*T[2];
	H[23] = R[0]*S[2]*T[2];
	H[24] = R[2]*S[2]*T[0];
	H[25] = R[2]*S[2]*T[1];
	H[26] = R[2]*S[2]*T[2];
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FEHex27_::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};
	double T[3] = {0.5*t*(t-1.0), 0.5*t*(t+1.0), 1.0 - t*t};

	double DR[3] = {r - 0.5, r  + 0.5, -2.0*r};
	double DS[3] = {s - 0.5, s  + 0.5, -2.0*s};
	double DT[3] = {t - 0.5, t  + 0.5, -2.0*t};

	Hr[ 0] = DR[0]*S[0]*T[0];
	Hr[ 1] = DR[1]*S[0]*T[0];
	Hr[ 2] = DR[1]*S[1]*T[0];
	Hr[ 3] = DR[0]*S[1]*T[0];
	Hr[ 4] = DR[0]*S[0]*T[1];
	Hr[ 5] = DR[1]*S[0]*T[1];
	Hr[ 6] = DR[1]*S[1]*T[1];
	Hr[ 7] = DR[0]*S[1]*T[1];
	Hr[ 8] = DR[2]*S[0]*T[0];
	Hr[ 9] = DR[1]*S[2]*T[0];
	Hr[10] = DR[2]*S[1]*T[0];
	Hr[11] = DR[0]*S[2]*T[0];
	Hr[12] = DR[2]*S[0]*T[1];
	Hr[13] = DR[1]*S[2]*T[1];
	Hr[14] = DR[2]*S[1]*T[1];
	Hr[15] = DR[0]*S[2]*T[1];
	Hr[16] = DR[0]*S[0]*T[2];
	Hr[17] = DR[1]*S[0]*T[2];
	Hr[18] = DR[1]*S[1]*T[2];
	Hr[19] = DR[0]*S[1]*T[2];
	Hr[20] = DR[2]*S[0]*T[2];
	Hr[21] = DR[1]*S[2]*T[2];
	Hr[22] = DR[2]*S[1]*T[2];
	Hr[23] = DR[0]*S[2]*T[2];
	Hr[24] = DR[2]*S[2]*T[0];
	Hr[25] = DR[2]*S[2]*T[1];
	Hr[26] = DR[2]*S[2]*T[2];

	Hs[ 0] = R[0]*DS[0]*T[0];
	Hs[ 1] = R[1]*DS[0]*T[0];
	Hs[ 2] = R[1]*DS[1]*T[0];
	Hs[ 3] = R[0]*DS[1]*T[0];
	Hs[ 4] = R[0]*DS[0]*T[1];
	Hs[ 5] = R[1]*DS[0]*T[1];
	Hs[ 6] = R[1]*DS[1]*T[1];
	Hs[ 7] = R[0]*DS[1]*T[1];
	Hs[ 8] = R[2]*DS[0]*T[0];
	Hs[ 9] = R[1]*DS[2]*T[0];
	Hs[10] = R[2]*DS[1]*T[0];
	Hs[11] = R[0]*DS[2]*T[0];
	Hs[12] = R[2]*DS[0]*T[1];
	Hs[13] = R[1]*DS[2]*T[1];
	Hs[14] = R[2]*DS[1]*T[1];
	Hs[15] = R[0]*DS[2]*T[1];
	Hs[16] = R[0]*DS[0]*T[2];
	Hs[17] = R[1]*DS[0]*T[2];
	Hs[18] = R[1]*DS[1]*T[2];
	Hs[19] = R[0]*DS[1]*T[2];
	Hs[20] = R[2]*DS[0]*T[2];
	Hs[21] = R[1]*DS[2]*T[2];
	Hs[22] = R[2]*DS[1]*T[2];
	Hs[23] = R[0]*DS[2]*T[2];
	Hs[24] = R[2]*DS[2]*T[0];
	Hs[25] = R[2]*DS[2]*T[1];
	Hs[26] = R[2]*DS[2]*T[2];

	Ht[ 0] = R[0]*S[0]*DT[0];
	Ht[ 1] = R[1]*S[0]*DT[0];
	Ht[ 2] = R[1]*S[1]*DT[0];
	Ht[ 3] = R[0]*S[1]*DT[0];
	Ht[ 4] = R[0]*S[0]*DT[1];
	Ht[ 5] = R[1]*S[0]*DT[1];
	Ht[ 6] = R[1]*S[1]*DT[1];
	Ht[ 7] = R[0]*S[1]*DT[1];
	Ht[ 8] = R[2]*S[0]*DT[0];
	Ht[ 9] = R[1]*S[2]*DT[0];
	Ht[10] = R[2]*S[1]*DT[0];
	Ht[11] = R[0]*S[2]*DT[0];
	Ht[12] = R[2]*S[0]*DT[1];
	Ht[13] = R[1]*S[2]*DT[1];
	Ht[14] = R[2]*S[1]*DT[1];
	Ht[15] = R[0]*S[2]*DT[1];
	Ht[16] = R[0]*S[0]*DT[2];
	Ht[17] = R[1]*S[0]*DT[2];
	Ht[18] = R[1]*S[1]*DT[2];
	Ht[19] = R[0]*S[1]*DT[2];
	Ht[20] = R[2]*S[0]*DT[2];
	Ht[21] = R[1]*S[2]*DT[2];
	Ht[22] = R[2]*S[1]*DT[2];
	Ht[23] = R[0]*S[2]*DT[2];
	Ht[24] = R[2]*S[2]*DT[0];
	Ht[25] = R[2]*S[2]*DT[1];
	Ht[26] = R[2]*S[2]*DT[2];
}

//-----------------------------------------------------------------------------
//! values of shape function second derivatives
//! \todo implement this (needed for biphasic problems)
void FEHex27_::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	
}

//=============================================================================
//              H E X 2 7 G 2 7
//=============================================================================

FEHex27G27::FEHex27G27() : FEHex27_(NINT, FE_HEX27G27)
{
	// integration point coordinates
	const double a = 0.774596669241483;
	const double w1 = 5.0 / 9.0;
	const double w2 = 8.0 / 9.0;
	gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -a; gw[ 0] = w1*w1*w1;
	gr[ 1] =  0; gs[ 1] = -a; gt[ 1] = -a; gw[ 1] = w2*w1*w1;
	gr[ 2] =  a; gs[ 2] = -a; gt[ 2] = -a; gw[ 2] = w1*w1*w1;
	gr[ 3] = -a; gs[ 3] =  0; gt[ 3] = -a; gw[ 3] = w1*w2*w1;
	gr[ 4] =  0; gs[ 4] =  0; gt[ 4] = -a; gw[ 4] = w2*w2*w1;
	gr[ 5] =  a; gs[ 5] =  0; gt[ 5] = -a; gw[ 5] = w1*w2*w1;
	gr[ 6] = -a; gs[ 6] =  a; gt[ 6] = -a; gw[ 6] = w1*w1*w1;
	gr[ 7] =  0; gs[ 7] =  a; gt[ 7] = -a; gw[ 7] = w2*w1*w1;
	gr[ 8] =  a; gs[ 8] =  a; gt[ 8] = -a; gw[ 8] = w1*w1*w1;
	gr[ 9] = -a; gs[ 9] = -a; gt[ 9] =  0; gw[ 9] = w1*w1*w2;
	gr[10] =  0; gs[10] = -a; gt[10] =  0; gw[10] = w2*w1*w2;
	gr[11] =  a; gs[11] = -a; gt[11] =  0; gw[11] = w1*w1*w2;
	gr[12] = -a; gs[12] =  0; gt[12] =  0; gw[12] = w1*w2*w2;
	gr[13] =  0; gs[13] =  0; gt[13] =  0; gw[13] = w2*w2*w2;
	gr[14] =  a; gs[14] =  0; gt[14] =  0; gw[14] = w1*w2*w2;
	gr[15] = -a; gs[15] =  a; gt[15] =  0; gw[15] = w1*w1*w2;
	gr[16] =  0; gs[16] =  a; gt[16] =  0; gw[16] = w2*w1*w2;
	gr[17] =  a; gs[17] =  a; gt[17] =  0; gw[17] = w1*w1*w2;
	gr[18] = -a; gs[18] = -a; gt[18] =  a; gw[18] = w1*w1*w1;
	gr[19] =  0; gs[19] = -a; gt[19] =  a; gw[19] = w2*w1*w1;
	gr[20] =  a; gs[20] = -a; gt[20] =  a; gw[20] = w1*w1*w1;
	gr[21] = -a; gs[21] =  0; gt[21] =  a; gw[21] = w1*w2*w1;
	gr[22] =  0; gs[22] =  0; gt[22] =  a; gw[22] = w2*w2*w1;
	gr[23] =  a; gs[23] =  0; gt[23] =  a; gw[23] = w1*w2*w1;
	gr[24] = -a; gs[24] =  a; gt[24] =  a; gw[24] = w1*w1*w1;
	gr[25] =  0; gs[25] =  a; gt[25] =  a; gw[25] = w2*w1*w1;
	gr[26] =  a; gs[26] =  a; gt[26] =  a; gw[26] = w1*w1*w1;

	init();
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FEHex27G27::project_to_nodes(double* ai, double* ao)
{
	
}

//=============================================================================
//
//                  S U R F A C E   E L E M E N T S
//
//=============================================================================

FESurfaceElementTraits::FESurfaceElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et) : FEElementTraits(ni, ne, FE_ELEM_SURFACE, es, et)
{
	gr.resize(ni);
	gs.resize(ni);
	gw.resize(ni);

	Gr.resize(ni, ne);
	Gs.resize(ni, ne);
}

//-----------------------------------------------------------------------------
//! Initialize the surface element traits data variables.
//
void FESurfaceElementTraits::init()
{
	assert(nint > 0);
	assert(neln > 0);

	// evaluate shape functions
	const int NE = FEElement::MAX_NODES;
	double N[NE];
	for (int n=0; n<nint; ++n)
	{
		shape(N, gr[n], gs[n]);
		for (int i=0; i<neln; ++i) H[n][i] = N[i]; 
	}
	
	// evaluate shape function derivatives
	double Nr[NE], Ns[NE];
	for (int n=0; n<nint; ++n)
	{
		shape_deriv(Nr, Ns, gr[n], gs[n]);
		for (int i=0; i<neln; ++i)
		{
			Gr[n][i] = Nr[i];
			Gs[n][i] = Ns[i];
		}
	}
}

//=============================================================================
//                          F E Q U A D 4
//=============================================================================

//-----------------------------------------------------------------------------
void FEQuad4_::shape(double* H, double r, double s)
{
	H[0] = 0.25*(1-r)*(1-s);
	H[1] = 0.25*(1+r)*(1-s);
	H[2] = 0.25*(1+r)*(1+s);
	H[3] = 0.25*(1-r)*(1+s);
}

//-----------------------------------------------------------------------------
void FEQuad4_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
	Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
	Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
	Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
}

//-----------------------------------------------------------------------------
void FEQuad4_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] =  0.25; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = -0.25; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] =  0.25; Hss[2] = 0;
	Hrr[3] = 0; Hrs[3] = -0.25; Hss[3] = 0;
}

//=============================================================================
//                          F E Q U A D G 4 
//=============================================================================

FEQuad4G4::FEQuad4G4() : FEQuad4_(NINT, FE_QUAD4G4) 
{
	const double a = 1.0 / sqrt(3.0);
	gr[0] = -a; gs[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gw[3] = 1;
	init(); 
	Hi = H.inverse();
}

//-----------------------------------------------------------------------------
void FEQuad4G4::project_to_nodes(double* ai, double* ao)
{
	int ni = NINT;
	int ne = NELN;
	assert(ni == ne);
	for (int i=0; i<ne; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<ni; ++j) ao[i] += Hi[i][j]*ai[j];
	}
}

//=============================================================================
//                          F E Q U A D N I
//=============================================================================

FEQuad4NI::FEQuad4NI() : FEQuad4_(NINT, FE_QUAD4NI) 
{
	gr[0] = -1; gs[0] = -1; gw[0] = 1;
	gr[1] =  1; gs[1] = -1; gw[1] = 1;
	gr[2] =  1; gs[2] =  1; gw[2] = 1;
	gr[3] = -1; gs[3] =  1; gw[3] = 1;
	init(); 
}

//-----------------------------------------------------------------------------
void FEQuad4NI::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
	ao[3] = ai[3];
}

//=============================================================================
//                          F E T R I 
//=============================================================================

//-----------------------------------------------------------------------------
void FETri3_::shape(double* H, double r, double s)
{
	H[0] = 1.0 - r - s;
	H[1] = r;
	H[2] = s;
}

//-----------------------------------------------------------------------------
void FETri3_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -1; Hs[0] = -1;
	Hr[1] =  1; Hs[1] =  0;
	Hr[2] =  0; Hs[2] =  1;
}

//-----------------------------------------------------------------------------
void FETri3_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] = 0; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = 0; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] = 0; Hss[2] = 0;
}

//=============================================================================
//                          F E T R I G 1 
//=============================================================================

//-----------------------------------------------------------------------------
FETri3G1::FETri3G1() : FETri3_(NINT, FE_TRI3G1)
{
	const double a = 1.0/3.0;
	gr[0] = a; gs[0] = a; gw[0] = 0.5;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri3G1::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
}

//=============================================================================
//                          F E T R I G 3 
//=============================================================================

//-----------------------------------------------------------------------------
FETri3G3::FETri3G3() : FETri3_(NINT, FE_TRI3G3)
{
	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;
	init(); 
	Hi = H.inverse();
}

//-----------------------------------------------------------------------------
void FETri3G3::project_to_nodes(double* ai, double* ao)
{
	assert(NINT == NELN);
	for (int i=0; i<NELN; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<NINT; ++j) ao[i] += Hi[i][j]*ai[j];
	}
}

//=============================================================================
//                          F E T R I G 7
//=============================================================================

//-----------------------------------------------------------------------------
FETri3G7::FETri3G7() : FETri3_(NINT, FE_TRI3G7) 
{ 
	const double w = 1.0/2.0;
	gr[0] = 0.333333333333333; gs[0] = 0.333333333333333; gw[0] = w*0.225000000000000;
	gr[1] = 0.797426985353087; gs[1] = 0.101286507323456; gw[1] = w*0.125939180544827;
	gr[2] = 0.101286507323456; gs[2] = 0.797426985353087; gw[2] = w*0.125939180544827;
	gr[3] = 0.101286507323456; gs[3] = 0.101286507323456; gw[3] = w*0.125939180544827;
	gr[4] = 0.470142064105115; gs[4] = 0.470142064105115; gw[4] = w*0.132394152788506;
	gr[5] = 0.470142064105115; gs[5] = 0.059715871789770; gw[5] = w*0.132394152788506;
	gr[6] = 0.059715871789770; gs[6] = 0.470142064105115; gw[6] = w*0.132394152788506;
	init(); 

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri3G7::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}

//=============================================================================
//                          F E T R I N I
//=============================================================================

//-----------------------------------------------------------------------------
FETri3NI::FETri3NI() : FETri3_(NINT, FE_TRI3NI)
{ 
	const double a = 1.0 / 6.0;
	gr[0] = 0; gs[0] = 0; gw[0] = a;
	gr[1] = 1; gs[1] = 0; gw[1] = a;
	gr[2] = 0; gs[2] = 1; gw[2] = a;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri3NI::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
}

//============================================================================
//                             F E T R I 6
//============================================================================

//-----------------------------------------------------------------------------
void FETri6_::shape(double* H, double r, double s)
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
void FETri6_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -3.0 + 4.0*r + 4.0*s;
	Hr[1] =  4.0*r - 1.0;
	Hr[2] =  0.0;
	Hr[3] =  4.0 - 8.0*r - 4.0*s;
	Hr[4] =  4.0*s;
	Hr[5] = -4.0*s;

	Hs[0] = -3.0 + 4.0*s + 4.0*r;
	Hs[1] =  0.0;
	Hs[2] =  4.0*s - 1.0;
	Hs[3] = -4.0*r;
	Hs[4] =  4.0*r;
	Hs[5] =  4.0 - 8.0*s - 4.0*r;
}

//-----------------------------------------------------------------------------
void FETri6_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] =  4.0; Hrs[0] =  4.0; Hss[0] =  4.0;
	Hrr[1] =  4.0; Hrs[1] =  0.0; Hss[1] =  0.0;
	Hrr[2] =  0.0; Hrs[2] =  0.0; Hss[2] =  4.0;
	Hrr[3] = -8.0; Hrs[3] = -4.0; Hss[3] =  0.0;
	Hrr[4] =  0.0; Hrs[4] =  4.0; Hss[4] =  0.0;
	Hrr[5] =  0.0; Hrs[5] = -4.0; Hss[5] = -8.0;
}

//=============================================================================
//                          F E T R I 6 G 3
//=============================================================================

FETri6G3::FETri6G3() : FETri6_(NINT, FE_TRI6G3) 
{ 
	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri6G3::project_to_nodes(double* ai, double* ao)
{
	matrix H(3, 3);
	for (int n=0; n<3; ++n)
	{
		H[n][0] = 1.0 - gr[n] - gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}
	H.inverse();

	for (int i=0; i<3; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<3; ++j) ao[i] += H[i][j]*ai[j];
	}

	ao[3] = 0.5*(ao[0] + ao[1]);
	ao[4] = 0.5*(ao[1] + ao[2]);
	ao[5] = 0.5*(ao[2] + ao[0]);
}

//=============================================================================
//                          F E T R I 6 G 4
//=============================================================================

FETri6G4::FETri6G4() : FETri6_(NINT, FE_TRI6G4) 
{ 
	const double a = 1.0/3.0;
	const double b = 1.0/5.0;
	const double c = 3.0/5.0;
	gr[0] = a; gs[0] = a; gw[0] = -27.0/96.0;
	gr[1] = c; gs[1] = b; gw[1] =  25.0/96.0;
	gr[2] = b; gs[2] = b; gw[2] =  25.0/96.0;
	gr[3] = b; gs[3] = c; gw[3] =  25.0/96.0;
	init(); 
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FETri6G4::project_to_nodes(double* ai, double* ao)
{
	
}

//=============================================================================
//                          F E T R I 6 G 7
//=============================================================================

FETri6G7::FETri6G7() : FETri6_(NINT, FE_TRI6G7) 
{ 
	const double w = 1.0/2.0;
	gr[0] = 0.333333333333333; gs[0] = 0.333333333333333; gw[0] = w*0.225000000000000;
	gr[1] = 0.797426985353087; gs[1] = 0.101286507323456; gw[1] = w*0.125939180544827;
	gr[2] = 0.101286507323456; gs[2] = 0.797426985353087; gw[2] = w*0.125939180544827;
	gr[3] = 0.101286507323456; gs[3] = 0.101286507323456; gw[3] = w*0.125939180544827;
	gr[4] = 0.470142064105115; gs[4] = 0.470142064105115; gw[4] = w*0.132394152788506;
	gr[5] = 0.470142064105115; gs[5] = 0.059715871789770; gw[5] = w*0.132394152788506;
	gr[6] = 0.059715871789770; gs[6] = 0.470142064105115; gw[6] = w*0.132394152788506;
	init(); 

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri6G7::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}

//=============================================================================
//                          T R I 6 G L 7
//=============================================================================

FETri6GL7::FETri6GL7() : FETri6_(NINT, FE_TRI6GL7) 
{ 
	const double a = 1.0/40.0;
	const double b = 1.0/15.0;
	gr[0] = 0.0; gs[0] = 0.0; gw[0] = a;
	gr[1] = 1.0; gs[1] = 0.0; gw[1] = a;
	gr[2] = 0.0; gs[2] = 1.0; gw[2] = a;
	gr[3] = 0.5; gs[3] = 0.0; gw[3] = b;
	gr[4] = 0.5; gs[4] = 0.5; gw[4] = b;
	gr[5] = 0.0; gs[5] = 0.5; gw[5] = b;
	gr[6] = 1.0/3.0; gs[6] = 1.0/3.0; gw[6] = 9.0*a;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri6GL7::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0]; ao[1] = ai[1]; ao[2] = ai[2];
	ao[3] = ai[3]; ao[4] = ai[4]; ao[5] = ai[5];
}

//=============================================================================
//                          F E T R I 6 N I
//=============================================================================

FETri6NI::FETri6NI() : FETri6_(NINT, FE_TRI6NI)
{ 
	const double a = 0.0;
	const double b = 1.0/6.0;
	gr[0] = 0.0; gs[0] = 0.0; gw[0] = a;
	gr[1] = 1.0; gs[1] = 0.0; gw[1] = a;
	gr[2] = 0.0; gs[2] = 1.0; gw[2] = a;
	gr[3] = 0.5; gs[3] = 0.0; gw[3] = b;
	gr[4] = 0.5; gs[4] = 0.5; gw[4] = b;
	gr[5] = 0.0; gs[5] = 0.5; gw[5] = b;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri6NI::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
	ao[3] = ai[3];
	ao[4] = ai[4];
	ao[5] = ai[5];
}

//============================================================================
//                             F E T R I 6 M
//============================================================================

// parameter used in the tri6m (6-node triangle with modified shape functions)
const double fetri6m_alpha = 0.2;

//-----------------------------------------------------------------------------
void FETri6m_::shape(double* H, double r, double s)
{
	double r1 = 1.0 - r - s;
	double r2 = r;
	double r3 = s;

	double N[6];
	N[0] = r1*(2.0*r1 - 1.0);
	N[1] = r2*(2.0*r2 - 1.0);
	N[2] = r3*(2.0*r3 - 1.0);
	N[3] = 4.0*r1*r2;
	N[4] = 4.0*r2*r3;
	N[5] = 4.0*r3*r1;

	const double a = fetri6m_alpha;
	const double b = 1.0 - 2.0*a;
	H[0] = N[0] + a*(N[3] + N[5]);
	H[1] = N[1] + a*(N[3] + N[4]);
	H[2] = N[2] + a*(N[4] + N[5]);
	H[3] = b*N[3];
	H[4] = b*N[4];
	H[5] = b*N[5];
}

//-----------------------------------------------------------------------------
void FETri6m_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	double Nr[6], Ns[6];
	Nr[0] = -3.0 + 4.0*r + 4.0*s;
	Nr[1] =  4.0*r - 1.0;
	Nr[2] =  0.0;
	Nr[3] =  4.0 - 8.0*r - 4.0*s;
	Nr[4] =  4.0*s;
	Nr[5] = -4.0*s;

	Ns[0] = -3.0 + 4.0*s + 4.0*r;
	Ns[1] =  0.0;
	Ns[2] =  4.0*s - 1.0;
	Ns[3] = -4.0*r;
	Ns[4] =  4.0*r;
	Ns[5] =  4.0 - 8.0*s - 4.0*r;

	const double a = fetri6m_alpha;
	const double b = 1.0 - 2.0*a;
	Hr[0] = Nr[0] + a*(Nr[3] + Nr[5]);
	Hr[1] = Nr[1] + a*(Nr[3] + Nr[4]);
	Hr[2] = Nr[2] + a*(Nr[4] + Nr[5]);
	Hr[3] = b*Nr[3];
	Hr[4] = b*Nr[4];
	Hr[5] = b*Nr[5];

	Hs[0] = Ns[0] + a*(Ns[3] + Ns[5]);
	Hs[1] = Ns[1] + a*(Ns[3] + Ns[4]);
	Hs[2] = Ns[2] + a*(Ns[4] + Ns[5]);
	Hs[3] = b*Ns[3];
	Hs[4] = b*Ns[4];
	Hs[5] = b*Ns[5];
}

//-----------------------------------------------------------------------------
void FETri6m_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	double Nrr[6], Nrs[6], Nss[6];
	Nrr[0] =  4.0; Nrs[0] =  4.0; Nss[0] =  4.0;
	Nrr[1] =  4.0; Nrs[1] =  0.0; Nss[1] =  0.0;
	Nrr[2] =  0.0; Nrs[2] =  0.0; Nss[2] =  4.0;
	Nrr[3] = -8.0; Nrs[3] = -4.0; Nss[3] =  0.0;
	Nrr[4] =  0.0; Nrs[4] =  4.0; Nss[4] =  0.0;
	Nrr[5] =  0.0; Nrs[5] = -4.0; Nss[5] = -8.0;

	const double a = fetri6m_alpha;
	const double b = 1.0 - 2.0*a;
	Hrr[0] = Nrr[0] + a*(Nrr[3] + Nrr[5]);
	Hrr[1] = Nrr[1] + a*(Nrr[3] + Nrr[4]);
	Hrr[2] = Nrr[2] + a*(Nrr[4] + Nrr[5]);
	Hrr[3] = b*Nrr[3];
	Hrr[4] = b*Nrr[4];
	Hrr[5] = b*Nrr[5];

	Hrs[0] = Nrs[0] + a*(Nrs[3] + Nrs[5]);
	Hrs[1] = Nrs[1] + a*(Nrs[3] + Nrs[4]);
	Hrs[2] = Nrs[2] + a*(Nrs[4] + Nrs[5]);
	Hrs[3] = b*Nrs[3];
	Hrs[4] = b*Nrs[4];
	Hrs[5] = b*Nrs[5];

	Hss[0] = Nss[0] + a*(Nss[3] + Nss[5]);
	Hss[1] = Nss[1] + a*(Nss[3] + Nss[4]);
	Hss[2] = Nss[2] + a*(Nss[4] + Nss[5]);
	Hss[3] = b*Nss[3];
	Hss[4] = b*Nss[4];
	Hss[5] = b*Nss[5];
}

//=============================================================================
//                          F E T R I 6 M G 7
//=============================================================================

FETri6mG7::FETri6mG7() : FETri6m_(NINT, FE_TRI6MG7) 
{ 
	const double w = 1.0/2.0;
	gr[0] = 0.333333333333333; gs[0] = 0.333333333333333; gw[0] = w*0.225000000000000;
	gr[1] = 0.797426985353087; gs[1] = 0.101286507323456; gw[1] = w*0.125939180544827;
	gr[2] = 0.101286507323456; gs[2] = 0.797426985353087; gw[2] = w*0.125939180544827;
	gr[3] = 0.101286507323456; gs[3] = 0.101286507323456; gw[3] = w*0.125939180544827;
	gr[4] = 0.470142064105115; gs[4] = 0.470142064105115; gw[4] = w*0.132394152788506;
	gr[5] = 0.470142064105115; gs[5] = 0.059715871789770; gw[5] = w*0.132394152788506;
	gr[6] = 0.059715871789770; gs[6] = 0.470142064105115; gw[6] = w*0.132394152788506;
	init(); 

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri6mG7::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}

//=============================================================================
//                          F E T R I 7 G 3
//=============================================================================

FETri7G3::FETri7G3() : FETri7_(NINT, FE_TRI7G3) 
{ 
	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri7G3::project_to_nodes(double* ai, double* ao)
{
	matrix H(3, 3);
	for (int n=0; n<3; ++n)
	{
		H[n][0] = 1.0 - gr[n] - gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}
	H.inverse();

	for (int i=0; i<3; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<3; ++j) ao[i] += H[i][j]*ai[j];
	}

	ao[3] = 0.5*(ao[0] + ao[1]);
	ao[4] = 0.5*(ao[1] + ao[2]);
	ao[5] = 0.5*(ao[2] + ao[0]);
	ao[6] = (ao[0]+ao[1]+ao[2])/3.0;
}

//=============================================================================
//                          F E T R I 6 G 4
//=============================================================================

FETri7G4::FETri7G4() : FETri7_(NINT, FE_TRI7G4) 
{ 
	const double a = 1.0/3.0;
	const double b = 1.0/5.0;
	const double c = 3.0/5.0;
	gr[0] = a; gs[0] = a; gw[0] = -27.0/96.0;
	gr[1] = c; gs[1] = b; gw[1] =  25.0/96.0;
	gr[2] = b; gs[2] = b; gw[2] =  25.0/96.0;
	gr[3] = b; gs[3] = c; gw[3] =  25.0/96.0;
	init(); 
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FETri7G4::project_to_nodes(double* ai, double* ao)
{
	
}

//============================================================================
//                             F E T R I 7
//============================================================================

//-----------------------------------------------------------------------------
void FETri7_::shape(double* H, double r, double s)
{
	double r1 = 1.0 - r - s;
	double r2 = r;
	double r3 = s;

	H[6] = 27.0*r1*r2*r3;
	H[0] = r1*(2.0*r1 - 1.0) + H[6]/9.0;
	H[1] = r2*(2.0*r2 - 1.0) + H[6]/9.0;
	H[2] = r3*(2.0*r3 - 1.0) + H[6]/9.0;
	H[3] = 4.0*r1*r2 - 4.0*H[6]/9.0;
	H[4] = 4.0*r2*r3 - 4.0*H[6]/9.0;
	H[5] = 4.0*r3*r1 - 4.0*H[6]/9.0;
}

//-----------------------------------------------------------------------------
void FETri7_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[6] = 27.0*s*(1.0 - 2.0*r - s);
	Hr[0] = -3.0 + 4.0*r + 4.0*s     + Hr[6]/9.0;
	Hr[1] =  4.0*r - 1.0             + Hr[6]/9.0;
	Hr[2] =  0.0                     + Hr[6]/9.0;
	Hr[3] =  4.0 - 8.0*r - 4.0*s - 4.0*Hr[6]/9.0;
	Hr[4] =  4.0*s               - 4.0*Hr[6]/9.0;
	Hr[5] = -4.0*s               - 4.0*Hr[6]/9.0;

	Hs[6] = 27.0*r*(1.0 - r - 2.0*s);
	Hs[0] = -3.0 + 4.0*s + 4.0*r     + Hs[6]/9.0;
	Hs[1] =  0.0                     + Hs[6]/9.0;
	Hs[2] =  4.0*s - 1.0             + Hs[6]/9.0;
	Hs[3] = -4.0*r               - 4.0*Hs[6]/9.0;
	Hs[4] =  4.0*r               - 4.0*Hs[6]/9.0;
	Hs[5] =  4.0 - 8.0*s - 4.0*r - 4.0*Hs[6]/9.0;
}

//-----------------------------------------------------------------------------
void FETri7_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[6] = -54.0*s;
	Hss[6] = -54.0*r;
	Hrs[6] = 27.0*(1.0 - 2.0*r - 2.0*s);

	Hrr[0] =  4.0 +     Hrr[6]/9.0; Hrs[0] =  4.0 +     Hrs[6]/9.0; Hss[0] =  4.0 +     Hss[6]/9.0;
	Hrr[1] =  4.0 +     Hrr[6]/9.0; Hrs[1] =  0.0 +     Hrs[6]/9.0; Hss[1] =  0.0 +     Hss[6]/9.0;
	Hrr[2] =  0.0 +     Hrr[6]/9.0; Hrs[2] =  0.0 +     Hrs[6]/9.0; Hss[2] =  4.0 +     Hss[6]/9.0;
	Hrr[3] = -8.0 - 4.0*Hrr[6]/9.0; Hrs[3] = -4.0 - 4.0*Hrs[6]/9.0; Hss[3] =  0.0 - 4.0*Hss[6]/9.0;
	Hrr[4] =  0.0 - 4.0*Hrr[6]/9.0; Hrs[4] =  4.0 - 4.0*Hrs[6]/9.0; Hss[4] =  0.0 - 4.0*Hss[6]/9.0;
	Hrr[5] =  0.0 - 4.0*Hrr[6]/9.0; Hrs[5] = -4.0 - 4.0*Hrs[6]/9.0; Hss[5] = -8.0 - 4.0*Hss[6]/9.0;
}


//=============================================================================
//                          F E T R I 7 G 7
//=============================================================================

FETri7G7::FETri7G7() : FETri7_(NINT, FE_TRI7G7) 
{ 
	const double w = 1.0/2.0;
	gr[0] = 0.333333333333333; gs[0] = 0.333333333333333; gw[0] = w*0.225000000000000;
	gr[1] = 0.797426985353087; gs[1] = 0.101286507323456; gw[1] = w*0.125939180544827;
	gr[2] = 0.101286507323456; gs[2] = 0.797426985353087; gw[2] = w*0.125939180544827;
	gr[3] = 0.101286507323456; gs[3] = 0.101286507323456; gw[3] = w*0.125939180544827;
	gr[4] = 0.470142064105115; gs[4] = 0.470142064105115; gw[4] = w*0.132394152788506;
	gr[5] = 0.470142064105115; gs[5] = 0.059715871789770; gw[5] = w*0.132394152788506;
	gr[6] = 0.059715871789770; gs[6] = 0.470142064105115; gw[6] = w*0.132394152788506;
	init(); 

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri7G7::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}


//=============================================================================
//                          F E T R I 7 G L 7
//=============================================================================

FETri7GL7::FETri7GL7() : FETri7_(NINT, FE_TRI7GL7) 
{ 
	const double a = 1.0/40.0;
	const double b = 1.0/15.0;
	gr[0] = 0.0; gs[0] = 0.0; gw[0] = a;
	gr[1] = 1.0; gs[1] = 0.0; gw[1] = a;
	gr[2] = 0.0; gs[2] = 1.0; gw[2] = a;
	gr[3] = 0.5; gs[3] = 0.0; gw[3] = b;
	gr[4] = 0.5; gs[4] = 0.5; gw[4] = b;
	gr[5] = 0.0; gs[5] = 0.5; gw[5] = b;
	gr[6] = 1.0/3.0; gs[6] = 1.0/3.0; gw[6] = 9.0*a;
	init(); 
}

//-----------------------------------------------------------------------------
void FETri7GL7::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0]; ao[1] = ai[1]; ao[2] = ai[2];
	ao[3] = ai[3]; ao[4] = ai[4]; ao[5] = ai[5];
	ao[6] = ai[6];
}

//=============================================================================
//          F E Q U A D 8 
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
void FEQuad8_::shape(double* H, double r, double s)
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
void FEQuad8_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[4] = -r*(1 - s);
	Hr[5] = 0.5*(1 - s*s);
	Hr[6] = -r*(1 + s);
	Hr[7] = -0.5*(1 - s*s);

	Hr[0] = -0.25*(1 - s) - 0.5*(Hr[4] + Hr[7]);
	Hr[1] =  0.25*(1 - s) - 0.5*(Hr[4] + Hr[5]);
	Hr[2] =  0.25*(1 + s) - 0.5*(Hr[5] + Hr[6]);
	Hr[3] = -0.25*(1 + s) - 0.5*(Hr[6] + Hr[7]);

	Hs[4] = -0.5*(1 - r*r);
	Hs[5] = -s*(1 + r);
	Hs[6] = 0.5*(1 - r*r);
	Hs[7] = -s*(1 - r);

	Hs[0] = -0.25*(1 - r) - 0.5*(Hs[4] + Hs[7]);
	Hs[1] = -0.25*(1 + r) - 0.5*(Hs[4] + Hs[5]);
	Hs[2] =  0.25*(1 + r) - 0.5*(Hs[5] + Hs[6]);
	Hs[3] =  0.25*(1 - r) - 0.5*(Hs[6] + Hs[7]);
}

//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
//! \todo implement this
void FEQuad8_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
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

	Hrr[0] = - 0.5*(Hrr[4] + Hrr[7]);
	Hrr[1] = - 0.5*(Hrr[4] + Hrr[5]);
	Hrr[2] = - 0.5*(Hrr[5] + Hrr[6]);
	Hrr[3] = - 0.5*(Hrr[6] + Hrr[7]);

	Hrs[0] =  0.25 - 0.5*(Hrs[4] + Hrs[7]);
	Hrs[1] = -0.25 - 0.5*(Hrs[4] + Hrs[5]);
	Hrs[2] =  0.25 - 0.5*(Hrs[5] + Hrs[6]);
	Hrs[3] = -0.25 - 0.5*(Hrs[6] + Hrs[7]);

	Hss[0] = - 0.5*(Hss[4] + Hss[7]);
	Hss[1] = - 0.5*(Hss[4] + Hss[5]);
	Hss[2] = - 0.5*(Hss[5] + Hss[6]);
	Hss[3] = - 0.5*(Hss[6] + Hss[7]);
}

//=============================================================================
//       F E Q U A D 8 G 9
//=============================================================================

FEQuad8G9::FEQuad8G9() : FEQuad8_(NINT, FE_QUAD8G9)
{
	// integration point coordinates
	const double a = sqrt(0.6);
	const double w1 = 25.0/81.0;
	const double w2 = 40.0/81.0;
	const double w3 = 64.0/81.0;
	gr[ 0] = -a; gs[ 0] = -a;  gw[ 0] = w1;
	gr[ 1] =  0; gs[ 1] = -a;  gw[ 1] = w2;
	gr[ 2] =  a; gs[ 2] = -a;  gw[ 2] = w1;
	gr[ 3] = -a; gs[ 3] =  0;  gw[ 3] = w2;
	gr[ 4] =  0; gs[ 4] =  0;  gw[ 4] = w3;
	gr[ 5] =  a; gs[ 5] =  0;  gw[ 5] = w2;
	gr[ 6] = -a; gs[ 6] =  a;  gw[ 6] = w1;
	gr[ 7] =  0; gs[ 7] =  a;  gw[ 7] = w2;
	gr[ 8] =  a; gs[ 8] =  a;  gw[ 8] = w1;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FEQuad8G9::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}


//=============================================================================
//          F E Q U A D 9 
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
void FEQuad9_::shape(double* H, double r, double s)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};

	H[0] = R[0]*S[0];
	H[1] = R[1]*S[0];
	H[2] = R[1]*S[1];
	H[3] = R[0]*S[1];
	H[4] = R[2]*S[0];
	H[5] = R[1]*S[2];
	H[6] = R[2]*S[1];
	H[7] = R[0]*S[2];
	H[8] = R[2]*S[2];
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
void FEQuad9_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};
	double DR[3] = {r-0.5, r+0.5, -2.0*r};
	double DS[3] = {s-0.5, s+0.5, -2.0*s};

	Hr[0] = DR[0]*S[0];
	Hr[1] = DR[1]*S[0];
	Hr[2] = DR[1]*S[1];
	Hr[3] = DR[0]*S[1];
	Hr[4] = DR[2]*S[0];
	Hr[5] = DR[1]*S[2];
	Hr[6] = DR[2]*S[1];
	Hr[7] = DR[0]*S[2];
	Hr[8] = DR[2]*S[2];

	Hs[0] = R[0]*DS[0];
	Hs[1] = R[1]*DS[0];
	Hs[2] = R[1]*DS[1];
	Hs[3] = R[0]*DS[1];
	Hs[4] = R[2]*DS[0];
	Hs[5] = R[1]*DS[2];
	Hs[6] = R[2]*DS[1];
	Hs[7] = R[0]*DS[2];
	Hs[8] = R[2]*DS[2];
}

//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
void FEQuad9_::shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};
	double DR[3] = {r-0.5, r+0.5, -2.0*r};
	double DS[3] = {s-0.5, s+0.5, -2.0*s};
	double DDR[3] = {1.0, 1.0, -2.0};
	double DDS[3] = {1.0, 1.0, -2.0};

	Grr[0] = DDR[0]*S[0]; Grs[0] = DR[0]*DS[0]; Gss[0] = R[0]*DDS[0];
	Grr[1] = DDR[1]*S[0]; Grs[1] = DR[1]*DS[0]; Gss[1] = R[1]*DDS[0];
	Grr[2] = DDR[1]*S[1]; Grs[2] = DR[1]*DS[1]; Gss[2] = R[1]*DDS[1];
	Grr[3] = DDR[0]*S[1]; Grs[3] = DR[0]*DS[1]; Gss[3] = R[0]*DDS[1];
	Grr[4] = DDR[2]*S[0]; Grs[4] = DR[2]*DS[0]; Gss[4] = R[2]*DDS[0];
	Grr[5] = DDR[1]*S[2]; Grs[5] = DR[1]*DS[2]; Gss[5] = R[1]*DDS[2];
	Grr[6] = DDR[2]*S[1]; Grs[6] = DR[2]*DS[1]; Gss[6] = R[2]*DDS[1];
	Grr[7] = DDR[0]*S[2]; Grs[7] = DR[0]*DS[2]; Gss[7] = R[0]*DDS[2];
	Grr[8] = DDR[2]*S[2]; Grs[8] = DR[2]*DS[2]; Gss[8] = R[2]*DDS[2];		
}

//=============================================================================
//       F E Q U A D 9 G 9
//=============================================================================

FEQuad9G9::FEQuad9G9() : FEQuad9_(NINT, FE_QUAD9G9)
{
	// integration point coordinates
	const double a = sqrt(0.6);
	const double w1 = 25.0/81.0;
	const double w2 = 40.0/81.0;
	const double w3 = 64.0/81.0;
	gr[ 0] = -a; gs[ 0] = -a;  gw[ 0] = w1;
	gr[ 1] =  0; gs[ 1] = -a;  gw[ 1] = w2;
	gr[ 2] =  a; gs[ 2] = -a;  gw[ 2] = w1;
	gr[ 3] = -a; gs[ 3] =  0;  gw[ 3] = w2;
	gr[ 4] =  0; gs[ 4] =  0;  gw[ 4] = w3;
	gr[ 5] =  a; gs[ 5] =  0;  gw[ 5] = w2;
	gr[ 6] = -a; gs[ 6] =  a;  gw[ 6] = w1;
	gr[ 7] =  0; gs[ 7] =  a;  gw[ 7] = w2;
	gr[ 8] =  a; gs[ 8] =  a;  gw[ 8] = w1;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FEQuad9G9::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}

//=============================================================================
//
//                      S H E L L   E L E M E N T S
//
//=============================================================================

FEShellElementTraits::FEShellElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et) : FEElementTraits(ni, ne, FE_ELEM_SHELL, es, et)
{
	gr.resize(ni);
	gs.resize(ni);
	gt.resize(ni);
	gw.resize(ni);

	Hr.resize(ni, ne);
	Hs.resize(ni, ne);

	D0.resize(ne);
	Dt.resize(ne);
}


//=============================================================================
//                          F E S H E L L Q U A D E L E M E N T
//=============================================================================

void FEShellQuadElementTraits::init()
{
	int n;
	
	const double a = 1.0 / sqrt(3.0);
	const double b = sqrt(3.0/5.0);
	const double w = 5.0 / 9.0;
	
	gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -b; gw[ 0] = w;
	gr[ 1] =  a; gs[ 1] = -a; gt[ 1] = -b; gw[ 1] = w;
	gr[ 2] =  a; gs[ 2] =  a; gt[ 2] = -b; gw[ 2] = w;
	gr[ 3] = -a; gs[ 3] =  a; gt[ 3] = -b; gw[ 3] = w;
	
	gr[ 4] = -a; gs[ 4] = -a; gt[ 4] =  0; gw[ 4] = 8.0/9.0;
	gr[ 5] =  a; gs[ 5] = -a; gt[ 5] =  0; gw[ 5] = 8.0/9.0;
	gr[ 6] =  a; gs[ 6] =  a; gt[ 6] =  0; gw[ 6] = 8.0/9.0;
	gr[ 7] = -a; gs[ 7] =  a; gt[ 7] =  0; gw[ 7] = 8.0/9.0;
	
	gr[ 8] = -a; gs[ 8] = -a; gt[ 8] =  b; gw[ 8] = w;
	gr[ 9] =  a; gs[ 9] = -a; gt[ 9] =  b; gw[ 9] = w;
	gr[10] =  a; gs[10] =  a; gt[10] =  b; gw[10] = w;
	gr[11] = -a; gs[11] =  a; gt[11] =  b; gw[11] = w;
	
	
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}

//	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Hr[n][0] = -0.25*(1-gs[n]);
		Hr[n][1] =  0.25*(1-gs[n]);
		Hr[n][2] =  0.25*(1+gs[n]);
		Hr[n][3] = -0.25*(1+gs[n]);
		
		Hs[n][0] = -0.25*(1-gr[n]);
		Hs[n][1] = -0.25*(1+gr[n]);
		Hs[n][2] =  0.25*(1+gr[n]);
		Hs[n][3] =  0.25*(1-gr[n]);
	}
}

//=============================================================================
//                          F E S H E L L T R I E L E M E N T
//=============================================================================

void FEShellTriElementTraits::init()
{
	int n;

	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	const double w = 5.0 / 9.0;

	gr[0] = a; gs[0] = a; gt[0] = -b; gw[0] = a*w;
	gr[1] = b; gs[1] = a; gt[1] = -b; gw[1] = a*w;
	gr[2] = a; gs[2] = b; gt[2] = -b; gw[2] = a*w;

	gr[3] = a; gs[3] = a; gt[3] =  0; gw[3] = a*w;
	gr[4] = b; gs[4] = a; gt[4] =  0; gw[4] = a*w;
	gr[5] = a; gs[5] = b; gt[5] =  0; gw[5] = a*w;

	gr[6] = a; gs[6] = a; gt[6] =  b; gw[6] = a*w;
	gr[7] = b; gs[7] = a; gt[7] =  b; gw[7] = a*w;
	gr[8] = a; gs[8] = b; gt[8] =  b; gw[8] = a*w;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}

//	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Hr[n][0] = -1;
		Hr[n][1] =  1;
		Hr[n][2] =  0;

		Hs[n][0] = -1;
		Hs[n][1] =  0;
		Hs[n][2] =  1;
	}
}

//=============================================================================
//               F E F E R G U S O N S H E L L E L E M E N T T R A I T S
//=============================================================================
void FEFergusonShellElementTraits::shape(double *H, double r)
{
    H[0] = SQR(r-1)*(2+r)/4;
    H[1] = SQR(r+1)*(2-r)/4;
    H[2] = SQR(r-1)*(r+1)/4;
    H[3] = SQR(r+1)*(r-1)/4;
}

void FEFergusonShellElementTraits::shape_deriv(double* Gr, double r)
{
    Gr[0] = (SQR(r)-1)*0.75;
    Gr[1] = -(SQR(r)-1)*0.75;
    Gr[2] = (r-1)*(3*r+1)/4;
    Gr[3] = (r+1)*(3*r-1)/4;
}

void FEFergusonShellElementTraits::shape_deriv2(double* Grr, double r)
{
    Grr[0] = 1.5*r;
    Grr[1] = -1.5*r;
    Grr[2] = (3*r-1)/2;
    Grr[3] = (3*r+1)/2;
}

//=============================================================================
//               F E F E R G U S O N S H E L L Q U A D E L E M E N T
//=============================================================================

void FEFergusonShellQuadElementTraits::init()
{
    int n;
    
    const double a = 1.0 / sqrt(3.0);
    const double b = sqrt(3.0/5.0);
    const double w = 5.0 / 9.0;
    double Fr[4], Fs[4];
    double Gr[4], Gs[4];
    
    gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -b; gw[ 0] = w;
    gr[ 1] =  a; gs[ 1] = -a; gt[ 1] = -b; gw[ 1] = w;
    gr[ 2] =  a; gs[ 2] =  a; gt[ 2] = -b; gw[ 2] = w;
    gr[ 3] = -a; gs[ 3] =  a; gt[ 3] = -b; gw[ 3] = w;
    
    gr[ 4] = -a; gs[ 4] = -a; gt[ 4] =  0; gw[ 4] = 8.0/9.0;
    gr[ 5] =  a; gs[ 5] = -a; gt[ 5] =  0; gw[ 5] = 8.0/9.0;
    gr[ 6] =  a; gs[ 6] =  a; gt[ 6] =  0; gw[ 6] = 8.0/9.0;
    gr[ 7] = -a; gs[ 7] =  a; gt[ 7] =  0; gw[ 7] = 8.0/9.0;
    
    gr[ 8] = -a; gs[ 8] = -a; gt[ 8] =  b; gw[ 8] = w;
    gr[ 9] =  a; gs[ 9] = -a; gt[ 9] =  b; gw[ 9] = w;
    gr[10] =  a; gs[10] =  a; gt[10] =  b; gw[10] = w;
    gr[11] = -a; gs[11] =  a; gt[11] =  b; gw[11] = w;
    
    
    for (n=0; n<NINT; ++n)
    {
        shape(Fr, gr[n]);
        shape(Fs, gs[n]);

        H[n][0] = Fr[0]*Fs[0];
        H[n][1] = Fr[1]*Fs[0];
        H[n][2] = Fr[1]*Fs[1];
        H[n][3] = Fr[0]*Fs[1];

        P[n][0] = Fr[2]*Fs[0];
        P[n][1] = Fr[3]*Fs[0];
        P[n][2] = Fr[3]*Fs[1];
        P[n][3] = Fr[2]*Fs[1];
        
        Q[n][0] = Fr[0]*Fs[2];
        Q[n][1] = Fr[1]*Fs[2];
        Q[n][2] = Fr[1]*Fs[3];
        Q[n][3] = Fr[0]*Fs[3];
        
        shape_deriv(Gr, gr[n]);
        shape_deriv(Gs, gs[n]);
        
        Hr[n][0] = Gr[0]*Fs[0];
        Hr[n][1] = Gr[1]*Fs[0];
        Hr[n][2] = Gr[1]*Fs[1];
        Hr[n][3] = Gr[0]*Fs[1];
        
        Hs[n][0] = Fr[0]*Gs[0];
        Hs[n][1] = Fr[1]*Gs[0];
        Hs[n][2] = Fr[1]*Gs[1];
        Hs[n][3] = Fr[0]*Gs[1];
        
        Pr[n][0] = Gr[2]*Fs[0];
        Pr[n][1] = Gr[3]*Fs[0];
        Pr[n][2] = Gr[3]*Fs[1];
        Pr[n][3] = Gr[2]*Fs[1];
        
        Ps[n][0] = Fr[2]*Gs[0];
        Ps[n][1] = Fr[3]*Gs[0];
        Ps[n][2] = Fr[3]*Gs[1];
        Ps[n][3] = Fr[2]*Gs[1];
        
        Qr[n][0] = Gr[0]*Fs[2];
        Qr[n][1] = Gr[1]*Fs[2];
        Qr[n][2] = Gr[1]*Fs[3];
        Qr[n][3] = Gr[0]*Fs[3];
        
        Qs[n][0] = Fr[0]*Gs[2];
        Qs[n][1] = Fr[1]*Gs[2];
        Qs[n][2] = Fr[1]*Gs[3];
        Qs[n][3] = Fr[0]*Gs[3];
    }
}

//=============================================================================
//                          F E T R U S S E L E M E N T
//=============================================================================

void FETrussElementTraits::init()
{

}

//=============================================================================
//
//                  2 D   E L E M E N T S
//
//=============================================================================

FE2DElementTraits::FE2DElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et) : FEElementTraits(ni, ne, FE_ELEM_2D, es, et)
{
	gr.resize(ni);
	gs.resize(ni);
	gw.resize(ni);

	Gr.resize(ni, ne);
	Gs.resize(ni, ne);
    
    Grr.resize(ni, ne);
    Gsr.resize(ni, ne);
    
    Grs.resize(ni, ne);
    Gss.resize(ni, ne);
}

//-----------------------------------------------------------------------------
//! Initialize the 2D element traits data variables.
//
void FE2DElementTraits::init()
{
	assert(nint > 0);
	assert(neln > 0);

	// evaluate shape functions
	const int NE = FEElement::MAX_NODES;
	double N[NE];
	for (int n=0; n<nint; ++n)
	{
		shape(N, gr[n], gs[n]);
		for (int i=0; i<neln; ++i) H[n][i] = N[i]; 
	}
	
	// evaluate shape function derivatives
	double Nr[NE], Ns[NE];
	for (int n=0; n<nint; ++n)
	{
		shape_deriv(Nr, Ns, gr[n], gs[n]);
		for (int i=0; i<neln; ++i)
		{
			Gr[n][i] = Nr[i];
			Gs[n][i] = Ns[i];
		}
	}
}

//=============================================================================
//                          F E 2 D T R I 
//=============================================================================

//-----------------------------------------------------------------------------
void FE2DTri3_::shape(double* H, double r, double s)
{
	H[0] = 1.0 - r - s;
	H[1] = r;
	H[2] = s;
}

//-----------------------------------------------------------------------------
void FE2DTri3_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -1; Hs[0] = -1;
	Hr[1] =  1; Hs[1] =  0;
	Hr[2] =  0; Hs[2] =  1;
}

//-----------------------------------------------------------------------------
void FE2DTri3_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] = 0; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = 0; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] = 0; Hss[2] = 0;
}

//=============================================================================
//                          F E 2 D T R I G 1 
//=============================================================================

//-----------------------------------------------------------------------------
FE2DTri3G1::FE2DTri3G1() : FE2DTri3_(NINT, FE2D_TRI3G1)
{
	const double a = 1.0/3.0;
	gr[0] = a; gs[0] = a; gw[0] = 0.5;
	init(); 
}

//-----------------------------------------------------------------------------
void FE2DTri3G1::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
}

//============================================================================
//                             F E 2 D T R I 6
//============================================================================

//-----------------------------------------------------------------------------
void FE2DTri6_::shape(double* H, double r, double s)
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
void FE2DTri6_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -3.0 + 4.0*r + 4.0*s;
	Hr[1] =  4.0*r - 1.0;
	Hr[2] =  0.0;
	Hr[3] =  4.0 - 8.0*r - 4.0*s;
	Hr[4] =  4.0*s;
	Hr[5] = -4.0*s;

	Hs[0] = -3.0 + 4.0*s + 4.0*r;
	Hs[1] =  0.0;
	Hs[2] =  4.0*s - 1.0;
	Hs[3] = -4.0*r;
	Hs[4] =  4.0*r;
	Hs[5] =  4.0 - 8.0*s - 4.0*r;
}

//-----------------------------------------------------------------------------
void FE2DTri6_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] =  4.0; Hrs[0] =  4.0; Hss[0] =  4.0;
	Hrr[1] =  4.0; Hrs[1] =  0.0; Hss[1] =  0.0;
	Hrr[2] =  0.0; Hrs[2] =  0.0; Hss[2] =  4.0;
	Hrr[3] = -8.0; Hrs[3] = -4.0; Hss[3] =  0.0;
	Hrr[4] =  0.0; Hrs[4] =  4.0; Hss[4] =  0.0;
	Hrr[5] =  0.0; Hrs[5] = -4.0; Hss[5] = -8.0;
}

//=============================================================================
//                          F E 2 D T R I 6 G 3
//=============================================================================

FE2DTri6G3::FE2DTri6G3() : FE2DTri6_(NINT, FE2D_TRI6G3)
{ 
	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;
	init(); 
}

//-----------------------------------------------------------------------------
void FE2DTri6G3::project_to_nodes(double* ai, double* ao)
{
	matrix H(3, 3);
	for (int n=0; n<3; ++n)
	{
		H[n][0] = 1.0 - gr[n] - gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}
	H.inverse();

	for (int i=0; i<3; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<3; ++j) ao[i] += H[i][j]*ai[j];
	}

	ao[3] = 0.5*(ao[0] + ao[1]);
	ao[4] = 0.5*(ao[1] + ao[2]);
	ao[5] = 0.5*(ao[2] + ao[0]);
}

//=============================================================================
//                          F E 2 D Q U A D 4
//=============================================================================

//-----------------------------------------------------------------------------
void FE2DQuad4_::shape(double* H, double r, double s)
{
	H[0] = 0.25*(1-r)*(1-s);
	H[1] = 0.25*(1+r)*(1-s);
	H[2] = 0.25*(1+r)*(1+s);
	H[3] = 0.25*(1-r)*(1+s);
}

//-----------------------------------------------------------------------------
void FE2DQuad4_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
	Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
	Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
	Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
}

//-----------------------------------------------------------------------------
void FE2DQuad4_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] =  0.25; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = -0.25; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] =  0.25; Hss[2] = 0;
	Hrr[3] = 0; Hrs[3] = -0.25; Hss[3] = 0;
}

//=============================================================================
//                          F E 2 D Q U A D G 4 
//=============================================================================

FE2DQuad4G4::FE2DQuad4G4() : FE2DQuad4_(NINT, FE2D_QUAD4G4) 
{
	const double a = 1.0 / sqrt(3.0);
	gr[0] = -a; gs[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gw[3] = 1;
	init(); 
	Hi = H.inverse();
}

//-----------------------------------------------------------------------------
void FE2DQuad4G4::project_to_nodes(double* ai, double* ao)
{
	int ni = NINT;
	int ne = NELN;
	assert(ni == ne);
	for (int i=0; i<ne; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<ni; ++j) ao[i] += Hi[i][j]*ai[j];
	}
}

//=============================================================================
//          F E 2 D Q U A D 8 
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
void FE2DQuad8_::shape(double* H, double r, double s)
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
void FE2DQuad8_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[4] = -r*(1 - s);
	Hr[5] = 0.5*(1 - s*s);
	Hr[6] = -r*(1 + s);
	Hr[7] = -0.5*(1 - s*s);

	Hr[0] = -0.25*(1 - s) - 0.5*(Hr[4] + Hr[7]);
	Hr[1] =  0.25*(1 - s) - 0.5*(Hr[4] + Hr[5]);
	Hr[2] =  0.25*(1 + s) - 0.5*(Hr[5] + Hr[6]);
	Hr[3] = -0.25*(1 + s) - 0.5*(Hr[6] + Hr[7]);

	Hs[4] = -0.5*(1 - r*r);
	Hs[5] = -s*(1 + r);
	Hs[6] = 0.5*(1 - r*r);
	Hs[7] = -s*(1 - r);

	Hs[0] = -0.25*(1 - r) - 0.5*(Hs[4] + Hs[7]);
	Hs[1] = -0.25*(1 + r) - 0.5*(Hs[4] + Hs[5]);
	Hs[2] =  0.25*(1 + r) - 0.5*(Hs[5] + Hs[6]);
	Hs[3] =  0.25*(1 - r) - 0.5*(Hs[6] + Hs[7]);
}

//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
void FE2DQuad8_::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
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

	Hrr[0] = - 0.5*(Hrr[4] + Hrr[7]);
	Hrr[1] = - 0.5*(Hrr[4] + Hrr[5]);
	Hrr[2] = - 0.5*(Hrr[5] + Hrr[6]);
	Hrr[3] = - 0.5*(Hrr[6] + Hrr[7]);

	Hrs[0] =  0.25 - 0.5*(Hrs[4] + Hrs[7]);
	Hrs[1] = -0.25 - 0.5*(Hrs[4] + Hrs[5]);
	Hrs[2] =  0.25 - 0.5*(Hrs[5] + Hrs[6]);
	Hrs[3] = -0.25 - 0.5*(Hrs[6] + Hrs[7]);

	Hss[0] = - 0.5*(Hss[4] + Hss[7]);
	Hss[1] = - 0.5*(Hss[4] + Hss[5]);
	Hss[2] = - 0.5*(Hss[5] + Hss[6]);
	Hss[3] = - 0.5*(Hss[6] + Hss[7]);
}

//=============================================================================
//       F E 2 D Q U A D 8 G 9
//=============================================================================

FE2DQuad8G9::FE2DQuad8G9() : FE2DQuad8_(NINT, FE2D_QUAD8G9)
{
	// integration point coordinates
	const double a = sqrt(0.6);
	const double w1 = 25.0/81.0;
	const double w2 = 40.0/81.0;
	const double w3 = 64.0/81.0;
	gr[ 0] = -a; gs[ 0] = -a;  gw[ 0] = w1;
	gr[ 1] =  0; gs[ 1] = -a;  gw[ 1] = w2;
	gr[ 2] =  a; gs[ 2] = -a;  gw[ 2] = w1;
	gr[ 3] = -a; gs[ 3] =  0;  gw[ 3] = w2;
	gr[ 4] =  0; gs[ 4] =  0;  gw[ 4] = w3;
	gr[ 5] =  a; gs[ 5] =  0;  gw[ 5] = w2;
	gr[ 6] = -a; gs[ 6] =  a;  gw[ 6] = w1;
	gr[ 7] =  0; gs[ 7] =  a;  gw[ 7] = w2;
	gr[ 8] =  a; gs[ 8] =  a;  gw[ 8] = w1;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FE2DQuad8G9::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}

//=============================================================================
//          F E 2 D Q U A D 9 
//=============================================================================

//-----------------------------------------------------------------------------
// shape function at (r,s)
void FE2DQuad9_::shape(double* H, double r, double s)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};

	H[0] = R[0]*S[0];
	H[1] = R[1]*S[0];
	H[2] = R[1]*S[1];
	H[3] = R[0]*S[1];
	H[4] = R[2]*S[0];
	H[5] = R[1]*S[2];
	H[6] = R[2]*S[1];
	H[7] = R[0]*S[2];
	H[8] = R[2]*S[2];
}

//-----------------------------------------------------------------------------
// shape function derivatives at (r,s)
void FE2DQuad9_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};
	double DR[3] = {r-0.5, r+0.5, -2.0*r};
	double DS[3] = {s-0.5, s+0.5, -2.0*s};

	Hr[0] = DR[0]*S[0];
	Hr[1] = DR[1]*S[0];
	Hr[2] = DR[1]*S[1];
	Hr[3] = DR[0]*S[1];
	Hr[4] = DR[2]*S[0];
	Hr[5] = DR[1]*S[2];
	Hr[6] = DR[2]*S[1];
	Hr[7] = DR[0]*S[2];
	Hr[8] = DR[2]*S[2];

	Hs[0] = R[0]*DS[0];
	Hs[1] = R[1]*DS[0];
	Hs[2] = R[1]*DS[1];
	Hs[3] = R[0]*DS[1];
	Hs[4] = R[2]*DS[0];
	Hs[5] = R[1]*DS[2];
	Hs[6] = R[2]*DS[1];
	Hs[7] = R[0]*DS[2];
	Hs[8] = R[2]*DS[2];
}

//-----------------------------------------------------------------------------
//! shape function derivatives at (r,s)
void FE2DQuad9_::shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
{
	double R[3] = {0.5*r*(r-1.0), 0.5*r*(r+1.0), 1.0 - r*r};
	double S[3] = {0.5*s*(s-1.0), 0.5*s*(s+1.0), 1.0 - s*s};
	double DR[3] = {r-0.5, r+0.5, -2.0*r};
	double DS[3] = {s-0.5, s+0.5, -2.0*s};
	double DDR[3] = {1.0, 1.0, -2.0};
	double DDS[3] = {1.0, 1.0, -2.0};

	Grr[0] = DDR[0]*S[0]; Grs[0] = DR[0]*DS[0]; Gss[0] = R[0]*DDS[0];
	Grr[1] = DDR[1]*S[0]; Grs[1] = DR[1]*DS[0]; Gss[1] = R[1]*DDS[0];
	Grr[2] = DDR[1]*S[1]; Grs[2] = DR[1]*DS[1]; Gss[2] = R[1]*DDS[1];
	Grr[3] = DDR[0]*S[1]; Grs[3] = DR[0]*DS[1]; Gss[3] = R[0]*DDS[1];
	Grr[4] = DDR[2]*S[0]; Grs[4] = DR[2]*DS[0]; Gss[4] = R[2]*DDS[0];
	Grr[5] = DDR[1]*S[2]; Grs[5] = DR[1]*DS[2]; Gss[5] = R[1]*DDS[2];
	Grr[6] = DDR[2]*S[1]; Grs[6] = DR[2]*DS[1]; Gss[6] = R[2]*DDS[1];
	Grr[7] = DDR[0]*S[2]; Grs[7] = DR[0]*DS[2]; Gss[7] = R[0]*DDS[2];
	Grr[8] = DDR[2]*S[2]; Grs[8] = DR[2]*DS[2]; Gss[8] = R[2]*DDS[2];		
}

//=============================================================================
//       F E 2 D Q U A D 9 G 9
//=============================================================================

FE2DQuad9G9::FE2DQuad9G9() : FE2DQuad9_(NINT, FE2D_QUAD9G9)
{
	// integration point coordinates
	const double a = sqrt(0.6);
	const double w1 = 25.0/81.0;
	const double w2 = 40.0/81.0;
	const double w3 = 64.0/81.0;
	gr[ 0] = -a; gs[ 0] = -a;  gw[ 0] = w1;
	gr[ 1] =  0; gs[ 1] = -a;  gw[ 1] = w2;
	gr[ 2] =  a; gs[ 2] = -a;  gw[ 2] = w1;
	gr[ 3] = -a; gs[ 3] =  0;  gw[ 3] = w2;
	gr[ 4] =  0; gs[ 4] =  0;  gw[ 4] = w3;
	gr[ 5] =  a; gs[ 5] =  0;  gw[ 5] = w2;
	gr[ 6] = -a; gs[ 6] =  a;  gw[ 6] = w1;
	gr[ 7] =  0; gs[ 7] =  a;  gw[ 7] = w2;
	gr[ 8] =  a; gs[ 8] =  a;  gw[ 8] = w1;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN,NELN);
	Ai.resize(NELN,NELN);
	A = H.transpose()*H;
	Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FE2DQuad9G9::project_to_nodes(double* ai, double* ao)
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += Ai[i][j]*b[j];
	}
}

//=============================================================================
//
//                  L I N E    E L E M E N T S
//
//=============================================================================

FELineElementTraits::FELineElementTraits(int ni, int ne, FE_Element_Shape es, FE_Element_Type et) : FEElementTraits(ni, ne, FE_ELEM_EDGE, es, et)
{
	gr.resize(ni);
	gw.resize(ni);
	Gr.resize(ni, ne);
    Grr.resize(ni, ne);
}

//-----------------------------------------------------------------------------
void FELineElementTraits::init()
{
	assert(nint > 0);
	assert(neln > 0);

	// evaluate shape functions
	const int NE = FEElement::MAX_NODES;
	double N[NE];
	for (int n=0; n<nint; ++n)
	{
		shape(N, gr[n]);
		for (int i=0; i<neln; ++i) H[n][i] = N[i]; 
	}
	
	// evaluate shape function derivatives
	double Nr[NE];
	for (int n=0; n<nint; ++n)
	{
		shape_deriv(Nr, gr[n]);
		for (int i=0; i<neln; ++i)
		{
			Gr[n][i] = Nr[i];
		}
	}
}

//=============================================================================
//                         FELine2_
//=============================================================================

//-----------------------------------------------------------------------------
void FELine2_::shape(double* H, double r)
{
	H[0] = 0.5*(1.0 - r);
	H[1] = 0.5*(1.0 + r);
}

//-----------------------------------------------------------------------------
void FELine2_::shape_deriv(double* Hr, double r)
{
	Hr[0] = -0.5;
	Hr[1] =  0.5;
}

//-----------------------------------------------------------------------------
void FELine2_::shape_deriv2(double* Hrr, double r)
{
	Hrr[0] = 0;
	Hrr[1] = 0;
}

//=============================================================================
//                          FELine2G1 
//=============================================================================

//-----------------------------------------------------------------------------
FELine2G1::FELine2G1() : FELine2_(NINT, FE_LINE2G1)
{
	gr[0] = 0.0; gw[0] = 2.0;
	init(); 
}

//-----------------------------------------------------------------------------
void FELine2G1::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
}
