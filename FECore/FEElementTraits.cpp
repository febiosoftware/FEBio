// FEElementTraits.cpp: implementation of the FEElementTraits class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElementTraits.h"
#include "FEElement.h"
#include "FEException.h"
//=============================================================================
FESolidElementTraits::FESolidElementTraits(int ni, int ne, FE_Element_Type et) : FEElementTraits(ni, ne, et) 
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
		
	m_Jt.resize(ni);
	m_Jti.resize(ni);
	m_detJt.resize(ni);

	m_J0.resize(ni);
	m_J0i.resize(ni);
	m_detJ0.resize(ni);
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
			Grr[n][i] = Hrr[i]; Grs[n][i] = Hrs[i]; Grr[n][i] = Hrr[i]; 
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
	
}

//=============================================================================
//              H E X 2 0 G 2 7
//=============================================================================

FEHex20G27::FEHex20G27() : FEHex20_(NINT, FE_HEX20G27)
{
	// integration point coordinates
	const double a = 0.774596669241483;
	const double w1 = 5.0 / 9.0;
	const double w2 = 8.9 / 9.0;
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
//
//                  S U R F A C E   E L E M E N T S
//
//=============================================================================

FESurfaceElementTraits::FESurfaceElementTraits(int ni, int ne, FE_Element_Type et) : FEElementTraits(ni, ne, et)
{
	m_ntype = et;

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
//                          F E Q U A D E L E M E N T
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
//                          F E Q U A D G 4 E L E M E N T
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
//                          F E N I Q U A D E L E M E N T
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
//                          F E T R I E L E M E N T
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
//                          F E T R I G 1 E L E M E N T
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
//                          F E T R I G 3 E L E M E N T
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

//-----------------------------------------------------------------------------
//! \todo implement this
void FEQuad8G9::project_to_nodes(double* ai, double* ao)
{
	
}

//=============================================================================
//                          F E N I T R I E L E M E N T
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
//                          F E T R I 6 G 3 E L E M E N T
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
//                          F E T R I 6 G 4 E L E M E N T
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
//                          F E T R I 6 G 7 E L E M E N T
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
//                          F E N I T R I 6 E L E M E N T
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
void FEQuad8_::shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
{
	
}

//=============================================================================
//       F E Q U A D 8 G 7
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
}

//=============================================================================
//
//                      S H E L L   E L E M E N T S
//
//=============================================================================

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
//                          F E T R U S S E L E M E N T
//=============================================================================

void FETrussElementTraits::init()
{

}
