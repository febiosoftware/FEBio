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
#include "FEElementTraits.h"
#include "FEElement.h"
#include "FEException.h"
#include "FESolidElementShape.h"
#include "FESurfaceElementShape.h"
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

FEElementTraits::FEElementTraits(int ni, int ne, FE_Element_Class c, FE_Element_Shape s, FE_Element_Type t)
{
	m_neln = ne;
	m_nint = ni;
	m_faces = 0;
	m_spec.eclass = c;
	m_spec.eshape = s;
	m_spec.etype  = t;
	m_H.resize(ni, ne);
}

//-----------------------------------------------------------------------------
//! project mat3ds integration point data to nodes
void FEElementTraits::project_to_nodes(mat3ds* si, mat3ds* so) const
{
	double ai[FEElement::MAX_INTPOINTS];
	double ao[FEElement::MAX_NODES];
	for (int i = 0; i<3; ++i) {
		for (int j = i; j<3; ++j) {
			for (int n = 0; n<m_nint; ++n) ai[n] = si[n](i, j);
			project_to_nodes(ai, ao);
			for (int n = 0; n<m_neln; ++n) so[n](i, j) = ao[n];
		}
	}
}

//-----------------------------------------------------------------------------
//! project mat3d integration point data to nodes
void FEElementTraits::project_to_nodes(mat3d* si, mat3d* so) const
{
	double ai[FEElement::MAX_INTPOINTS];
	double ao[FEElement::MAX_NODES];
	for (int i = 0; i<3; ++i) {
		for (int j = 0; j<3; ++j) {
			for (int n = 0; n<m_nint; ++n) ai[n] = si[n](i, j);
			project_to_nodes(ai, ao);
			for (int n = 0; n<m_neln; ++n) so[n](i, j) = ao[n];
		}
	}
}

//-----------------------------------------------------------------------------
//! project vec3d integration point data to nodes
void FEElementTraits::project_to_nodes(vec3d* si, vec3d* so) const
{
	double ai[FEElement::MAX_INTPOINTS];
	double ao[FEElement::MAX_NODES];

	for (int n = 0; n<m_nint; ++n) ai[n] = si[n].x; project_to_nodes(ai, ao); for (int n = 0; n<m_neln; ++n) so[n].x = ao[n];
	for (int n = 0; n<m_nint; ++n) ai[n] = si[n].y; project_to_nodes(ai, ao); for (int n = 0; n<m_neln; ++n) so[n].y = ao[n];
	for (int n = 0; n<m_nint; ++n) ai[n] = si[n].z; project_to_nodes(ai, ao); for (int n = 0; n<m_neln; ++n) so[n].z = ao[n];
}

//=============================================================================
FESolidElementTraits::FESolidElementTraits(int ni, int ne, FE_Element_Shape eshape, FE_Element_Type etype) : FEElementTraits(ni, ne, FE_ELEM_SOLID, eshape, etype) 
{
	m_shape = nullptr;

	gr.resize(ni);
	gs.resize(ni);
	gt.resize(ni);
	gw.resize(ni);

	m_Gr.resize(ni, ne);
	m_Gs.resize(ni, ne);
	m_Gt.resize(ni, ne);

	Grr.resize(ni, ne);
	Gsr.resize(ni, ne);
	Gtr.resize(ni, ne);
		
	Grs.resize(ni, ne);
	Gss.resize(ni, ne);
	Gts.resize(ni, ne);
		
	Grt.resize(ni, ne);
	Gst.resize(ni, ne);
	Gtt.resize(ni, ne);

	// TODO: Move this to the individual classes?
	m_faces = 0;
	switch (eshape)
	{
	case ET_TET4:
	case ET_TET5:
	case ET_TET10:
	case ET_TET15:
	case ET_TET20: 
		m_faces = 4;
		break;
	case ET_PENTA6:
	case ET_PENTA15:
		m_faces = 5;
		break;
	case ET_HEX8:
	case ET_HEX20:
	case ET_HEX27:
		m_faces = 6;
		break;
	case ET_PYRA5:
    case ET_PYRA13:
		m_faces = 5;
		break;
	default:
		assert(false);
	}
}

//-----------------------------------------------------------------------------
int FESolidElementTraits::ShapeFunctions(int order)
{
	FESolidElementShape* shape = m_shapeP[order];
	return (shape ? shape->nodes(): 0);
}

//-----------------------------------------------------------------------------
//! initialize element traits data
void FESolidElementTraits::init()
{
	assert(m_nint > 0);
	assert(m_neln > 0);
	const int NELN = FEElement::MAX_NODES;

	// get shape class
	m_shape = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(m_spec.eshape));
	assert(m_shape && (m_shape->shape() == m_spec.eshape));

	// calculate shape function values at gauss points
	double N[NELN];
	for (int n=0; n<m_nint; ++n)
	{
		m_shape->shape_fnc(N, gr[n], gs[n], gt[n]);
		for (int i=0; i<m_neln; ++i) m_H[n][i] = N[i];
	}

	// calculate local derivatives of shape functions at gauss points
	double Hr[NELN], Hs[NELN], Ht[NELN];
	for (int n=0; n<m_nint; ++n)
	{
		m_shape->shape_deriv(Hr, Hs, Ht, gr[n], gs[n], gt[n]);
		for (int i=0; i<m_neln; ++i)
		{
			m_Gr[n][i] = Hr[i];
			m_Gs[n][i] = Hs[i];
			m_Gt[n][i] = Ht[i];
		}
	}
	
	// calculate local second derivatives of shape functions at gauss points
	double Hrr[NELN], Hss[NELN], Htt[NELN], Hrs[NELN], Hst[NELN], Hrt[NELN];
	for (int n=0; n<m_nint; ++n)
	{
		m_shape->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, gr[n], gs[n], gt[n]);
		for (int i=0; i<m_neln; ++i)
		{
			Grr[n][i] = Hrr[i]; Grs[n][i] = Hrs[i]; Grt[n][i] = Hrt[i]; 
			Gsr[n][i] = Hrs[i]; Gss[n][i] = Hss[i]; Gst[n][i] = Hst[i]; 
			Gtr[n][i] = Hrt[i]; Gts[n][i] = Hst[i]; Gtt[n][i] = Htt[i]; 
		}
	}

	// NOTE: Below, is a new interface for dealing with mixed element formulations.
	//       This is still a work in progress.

	// Get the max interpolation order
	const int maxOrder = (int) m_shapeP.size() - 1;
	m_Hp.resize(maxOrder + 1);
	m_Gr_p.resize(maxOrder + 1);
	m_Gs_p.resize(maxOrder + 1);
	m_Gt_p.resize(maxOrder + 1);
	for (int i = 0; i <= maxOrder; ++i)
	{
		FESolidElementShape* shape = m_shapeP[i];
		matrix& H = m_Hp[i];
		matrix& Gr = m_Gr_p[i];
		matrix& Gs = m_Gs_p[i];
		matrix& Gt = m_Gt_p[i];
		if (i == 0)
		{
			H.resize(m_nint, 1);
			Gr.resize(m_nint, 1);
			Gs.resize(m_nint, 1);
			Gt.resize(m_nint, 1);
			for (int n = 0; n < m_nint; ++n)
			{
				H[n][0] = 1.0;
				Gr[n][0] = Gs[n][0] = Gt[n][0] = 0.0;
			}
		}
		else if (m_shapeP[i])
		{
			// get the nodes
			int neln = shape->nodes();

			// shape function values
			H.resize(m_nint, neln);
			for (int n = 0; n<m_nint; ++n)
			{
				m_shapeP[i]->shape_fnc(N, gr[n], gs[n], gt[n]);
				for (int j = 0; j<neln; ++j) H[n][j] = N[j];
			}

			// calculate local derivatives of shape functions at gauss points
			Gr.resize(m_nint, neln);
			Gs.resize(m_nint, neln);
			Gt.resize(m_nint, neln);
			for (int n = 0; n<m_nint; ++n)
			{
				shape->shape_deriv(Hr, Hs, Ht, gr[n], gs[n], gt[n]);
				for (int j = 0; j<neln; ++j)
				{
					Gr[n][j] = Hr[j];
					Gs[n][j] = Hs[j];
					Gt[n][j] = Ht[j];
				}
			}
		}
	}
}

//! values of shape functions
void FESolidElementTraits::shape_fnc(double* H, double r, double s, double t)
{
	return m_shape->shape_fnc(H, r, s, t); 
}

//! values of shape function derivatives
void FESolidElementTraits::shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t)
{
	return m_shape->shape_deriv(Hr, Hs, Ht, r, s, t);
}

//! values of shape function second derivatives
void FESolidElementTraits::shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t)
{
	return m_shape->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t);
}

//=============================================================================
//                                F E H E X 8
//=============================================================================

void FEHex8_::init()
{
	// allocate shape classes
	m_shapeP.resize(2);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_HEX8));

	// initialize base class
	FESolidElementTraits::init();
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
	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void FEHex8G8::project_to_nodes(double* ai, double* ao) const
{
	for (int j=0; j<NELN; ++j)
	{
		ao[j] = 0;
		for (int k=0; k<NINT; ++k) 
		{
			ao[j] += m_Hi[j][k]*ai[k];
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
void FEHex8RI::project_to_nodes(double* ai, double* ao) const
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
void FEHex8G1::project_to_nodes(double* ai, double* ao) const
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

//! initialize element traits data
void FETet4_::init()
{
	// allocate shape classes
	m_shapeP.resize(2);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TET4));

	// initialize base class
	FESolidElementTraits::init();
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
	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void FETet4G4::project_to_nodes(double* ai, double* ao) const
{
	for (int j=0; j<NELN; ++j)
	{
		ao[j] = 0;
		for (int k=0; k<NINT; ++k) 
		{
			ao[j] += m_Hi[j][k]*ai[k];
		}
	}
}

//=============================================================================
//                              F E T E T 5
//=============================================================================

//! initialize element traits data
void FETet5_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TET4));
	m_shapeP[2] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TET5));

	// initialize base class
	FESolidElementTraits::init();
}

//=============================================================================
//                          T E T 5 G 4
//=============================================================================

FETet5G4::FETet5G4() : FETet5_(NINT, FE_TET5G4)
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
}

//-----------------------------------------------------------------------------
void FETet5G4::project_to_nodes(double* ai, double* ao) const
{
	// TODO: implement this
	assert(false);
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
void FETet4G1::project_to_nodes(double* ai, double* ao) const
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
	ao[3] = ai[0];
}

//=============================================================================
//                       P E N T A 6
//=============================================================================

void FEPenta6_::init()
{
	// allocate shape classes
	m_shapeP.resize(2);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_PENTA6));

	// initialize base class
	FESolidElementTraits::init();
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

	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void FEPenta6G6::project_to_nodes(double* ai, double* ao) const
{
	for (int j=0; j<NELN; ++j)
	{
		ao[j] = 0;
		for (int k=0; k<NINT; ++k) 
		{
			ao[j] += m_Hi[j][k]*ai[k];
		}
	}
}

//=============================================================================
//                       P E N T A 1 5
//=============================================================================

//=============================================================================
//                          F E P E N T A 1 5 G 8
//=============================================================================

int FEPenta15G8::ni[NELN] = {};

FEPenta15G8::FEPenta15G8() : FEPenta15_(NINT, FE_PENTA15G8)
{
    const double a = 1.0/3.0;
    const double b = 1.0/5.0;
    const double c = 3.0/5.0;
    const double d = sqrt(a);
    gr[0] = a; gs[0] = a; gt[0] = -d; gw[0] = -27.0/96.0;
    gr[1] = c; gs[1] = b; gt[1] = -d; gw[1] =  25.0/96.0;
    gr[2] = b; gs[2] = b; gt[2] = -d; gw[2] =  25.0/96.0;
    gr[3] = b; gs[3] = c; gt[3] = -d; gw[3] =  25.0/96.0;
    gr[4] = a; gs[4] = a; gt[4] =  d; gw[4] = -27.0/96.0;
    gr[5] = c; gs[5] = b; gt[5] =  d; gw[5] =  25.0/96.0;
    gr[6] = b; gs[6] = b; gt[6] =  d; gw[6] =  25.0/96.0;
    gr[7] = b; gs[7] = c; gt[7] =  d; gw[7] =  25.0/96.0;
    
    init();
    
	m_MT.resize(NELN, NINT);
    for (int i=0; i<NINT; ++i)
        for (int n=0; n<NELN; ++n)
			m_MT(n,i) = m_H(i,n);
    
	m_Hi.resize(NELN, NELN);
	m_Hi = m_MT*m_MT.transpose();
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
//! Use least-squares extrapolation
void FEPenta15G8::project_to_nodes(double* ai, double* ao) const
{
    double v[NELN];
    for (int n=0; n<NELN; ++n) {
        v[n] = 0;
        for (int i=0; i<NINT; ++i) {
            v[n] += m_MT(n,i)*ai[i];
        }
    }
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*v[k];
        }
    }
}

//=============================================================================
//                          F E P E N T A 1 5 G 2 1
//=============================================================================

int FEPenta15G21::ni[NELN] = { 1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 16, 17, 18, 19, 20 };

FEPenta15G21::FEPenta15G21() : FEPenta15_(NINT, FE_PENTA15G21)
{
    const double w = 1.0/2.0;
    const double a = 0.774596669241483;
    const double w1 = 5.0 / 9.0;
    const double w2 = 8.0 / 9.0;
    gr[ 0] = 0.333333333333333; gs[ 0] = 0.333333333333333; gt[ 0] = -a; gw[ 0] = w*w1*0.225000000000000;
    gr[ 1] = 0.797426985353087; gs[ 1] = 0.101286507323456; gt[ 1] = -a; gw[ 1] = w*w1*0.125939180544827;
    gr[ 2] = 0.101286507323456; gs[ 2] = 0.797426985353087; gt[ 2] = -a; gw[ 2] = w*w1*0.125939180544827;
    gr[ 3] = 0.101286507323456; gs[ 3] = 0.101286507323456; gt[ 3] = -a; gw[ 3] = w*w1*0.125939180544827;
    gr[ 4] = 0.470142064105115; gs[ 4] = 0.470142064105115; gt[ 4] = -a; gw[ 4] = w*w1*0.132394152788506;
    gr[ 5] = 0.470142064105115; gs[ 5] = 0.059715871789770; gt[ 5] = -a; gw[ 5] = w*w1*0.132394152788506;
    gr[ 6] = 0.059715871789770; gs[ 6] = 0.470142064105115; gt[ 6] = -a; gw[ 6] = w*w1*0.132394152788506;
    gr[ 7] = 0.333333333333333; gs[ 7] = 0.333333333333333; gt[ 7] =  0; gw[ 7] = w*w2*0.225000000000000;
    gr[ 8] = 0.797426985353087; gs[ 8] = 0.101286507323456; gt[ 8] =  0; gw[ 8] = w*w2*0.125939180544827;
    gr[ 9] = 0.101286507323456; gs[ 9] = 0.797426985353087; gt[ 9] =  0; gw[ 9] = w*w2*0.125939180544827;
    gr[10] = 0.101286507323456; gs[10] = 0.101286507323456; gt[10] =  0; gw[10] = w*w2*0.125939180544827;
    gr[11] = 0.470142064105115; gs[11] = 0.470142064105115; gt[11] =  0; gw[11] = w*w2*0.132394152788506;
    gr[12] = 0.470142064105115; gs[12] = 0.059715871789770; gt[12] =  0; gw[12] = w*w2*0.132394152788506;
    gr[13] = 0.059715871789770; gs[13] = 0.470142064105115; gt[13] =  0; gw[13] = w*w2*0.132394152788506;
    gr[14] = 0.333333333333333; gs[14] = 0.333333333333333; gt[14] =  a; gw[14] = w*w1*0.225000000000000;
    gr[15] = 0.797426985353087; gs[15] = 0.101286507323456; gt[15] =  a; gw[15] = w*w1*0.125939180544827;
    gr[16] = 0.101286507323456; gs[16] = 0.797426985353087; gt[16] =  a; gw[16] = w*w1*0.125939180544827;
    gr[17] = 0.101286507323456; gs[17] = 0.101286507323456; gt[17] =  a; gw[17] = w*w1*0.125939180544827;
    gr[18] = 0.470142064105115; gs[18] = 0.470142064105115; gt[18] =  a; gw[18] = w*w1*0.132394152788506;
    gr[19] = 0.470142064105115; gs[19] = 0.059715871789770; gt[19] =  a; gw[19] = w*w1*0.132394152788506;
    gr[20] = 0.059715871789770; gs[20] = 0.470142064105115; gt[20] =  a; gw[20] = w*w1*0.132394152788506;

    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEPenta15G21::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//=============================================================================
//                           T E T 1 0
//=============================================================================

//! initialize element traits data
void FETet10_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TET4));
	m_shapeP[2] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TET10));

	// initialize base class
	FESolidElementTraits::init();
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
void FETet10G1::project_to_nodes(double* ai, double* ao) const
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
void FETet10G4::project_to_nodes(double* ai, double* ao) const
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
void FETet10G8::project_to_nodes(double* ai, double* ao) const
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
void FETet10GL11::project_to_nodes(double* ai, double* ao) const
{
	ao[0] = ai[0]; ao[1] = ai[1]; ao[2] = ai[2]; ao[3] = ai[3];
	ao[4] = ai[4]; ao[5] = ai[5]; ao[6] = ai[6]; ao[7] = ai[7]; ao[8] = ai[8]; ao[9] = ai[9];
}

//=============================================================================
//                           T E T 1 5
//=============================================================================

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
void FETet15G4::project_to_nodes(double* ai, double* ao) const
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
void FETet15G8::project_to_nodes(double* ai, double* ao) const
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
void FETet15G11::project_to_nodes(double* ai, double* ao) const
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
void FETet15G15::project_to_nodes(double* ai, double* ao) const
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
//                           T E T 2 0
//=============================================================================

//=============================================================================
//                          T E T 2 0 G 1 5
//=============================================================================
FETet20G15::FETet20G15() : FETet20_(NINT, FE_TET20G15)
{
	gr[0] = 0.25; gs[0] = 0.25; gt[0] = 0.25; gw[0] = 0.030283678097089;

	gr[1] = 0.333333333333333; gs[1] = 0.333333333333333; gt[1] = 0.333333333333333; gw[1] = 0.006026785714286;
	gr[2] = 0.000000000000000; gs[2] = 0.333333333333333; gt[2] = 0.333333333333333; gw[2] = 0.006026785714286;
	gr[3] = 0.333333333333333; gs[3] = 0.000000000000000; gt[3] = 0.333333333333333; gw[3] = 0.006026785714286;
	gr[4] = 0.333333333333333; gs[4] = 0.333333333333333; gt[4] = 0.000000000000000; gw[4] = 0.006026785714286;

	gr[5] = 0.090909090909091; gs[5] = 0.090909090909091; gt[5] = 0.090909090909091; gw[5] = 0.011645249086029;
	gr[6] = 0.727272727272727; gs[6] = 0.090909090909091; gt[6] = 0.090909090909091; gw[6] = 0.011645249086029;
	gr[7] = 0.090909090909091; gs[7] = 0.727272727272727; gt[7] = 0.090909090909091; gw[7] = 0.011645249086029;
	gr[8] = 0.090909090909091; gs[8] = 0.090909090909091; gt[8] = 0.727272727272727; gw[8] = 0.011645249086029;

	gr[9] = 0.433449846426336; gs[9] = 0.066550153573664; gt[9] = 0.066550153573664; gw[9] = 0.010949141561386;
	gr[10] = 0.066550153573664; gs[10] = 0.433449846426336; gt[10] = 0.066550153573664; gw[10] = 0.010949141561386;
	gr[11] = 0.066550153573664; gs[11] = 0.066550153573664; gt[11] = 0.433449846426336; gw[11] = 0.010949141561386;
	gr[12] = 0.066550153573664; gs[12] = 0.433449846426336; gt[12] = 0.433449846426336; gw[12] = 0.010949141561386;
	gr[13] = 0.433449846426336; gs[13] = 0.066550153573664; gt[13] = 0.433449846426336; gw[13] = 0.010949141561386;
	gr[14] = 0.433449846426336; gs[14] = 0.433449846426336; gt[14] = 0.066550153573664; gw[14] = 0.010949141561386;

	init();

	// setup the shape function matrix
	N.resize(15, 4);
	for (int i = 0; i<15; ++i)
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

void FETet20G15::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(4);
	for (int i = 0; i<4; ++i)
	{
		b[i] = 0;
		for (int j = 0; j<NINT; ++j) b[i] += N[j][i] * ai[j];
	}

	for (int i = 0; i<4; ++i)
	{
		ao[i] = 0.0;
		for (int j = 0; j<4; ++j) ao[i] += Ai[i][j] * b[j];
	}

	const double w1 = 1.0 / 3.0;
	const double w2 = 2.0 / 3.0;

	ao[ 4] = ao[0] * w2 + ao[1] * w1;
	ao[ 5] = ao[0] * w1 + ao[1] * w2;
	ao[ 6] = ao[1] * w2 + ao[2] * w1;
	ao[ 7] = ao[1] * w1 + ao[2] * w2;
	ao[ 8] = ao[0] * w2 + ao[2] * w1;
	ao[ 9] = ao[0] * w1 + ao[2] * w2;
	ao[10] = ao[0] * w2 + ao[3] * w1;
	ao[11] = ao[0] * w1 + ao[3] * w2;
	ao[12] = ao[1] * w2 + ao[3] * w1;
	ao[13] = ao[1] * w1 + ao[3] * w2;
	ao[14] = ao[2] * w2 + ao[3] * w1;
	ao[15] = ao[2] * w1 + ao[3] * w2;

	ao[16] = (ao[0] + ao[1] + ao[3]) / 3.0;
	ao[17] = (ao[1] + ao[2] + ao[3]) / 3.0;
	ao[18] = (ao[2] + ao[0] + ao[3]) / 3.0;
	ao[19] = (ao[0] + ao[1] + ao[2]) / 3.0;
}

//=============================================================================
//              H E X 2 0
//=============================================================================

int FEHex20_::Nodes(int order)
{
	switch (order)
	{
	case 2: return 20; break;
	case 1: return 8; break;
	case 0: return 1; break;
	default:
		assert(false);
		return 20;
	}
}

void FEHex20_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_HEX8));
	m_shapeP[2] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_HEX20));

	// initialize base class
	FESolidElementTraits::init();
}

//=============================================================================
//              H E X 2 0 G 8
//=============================================================================

int FEHex20G8::ni[NELN] = {};

FEHex20G8::FEHex20G8() : FEHex20_(NINT, FE_HEX20G8)
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
    
    MT.resize(NELN, NINT);
    for (int i=0; i<NINT; ++i)
        for (int n=0; n<NELN; ++n)
            MT(n,i) = m_H(i,n);
    
    Hi.resize(NELN, NELN);
    Hi = MT*MT.transpose();
    Hi = Hi.inverse();
}

//-----------------------------------------------------------------------------
//! Use least-squares extrapolation
void FEHex20G8::project_to_nodes(double* ai, double* ao) const
{
    double v[NELN];
    for (int n=0; n<NELN; ++n) {
        v[n] = 0;
        for (int i=0; i<NINT; ++i) {
            v[n] += MT(n,i)*ai[i];
        }
    }
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += Hi[j][k]*v[k];
        }
    }
}

//=============================================================================
//              H E X 2 0 G 2 7
//=============================================================================

int FEHex20G27::ni[NELN] = { 0, 1, 2, 3, 5, 6, 7, 8, 9, 11, 15, 17, 18, 19, 20, 21, 23, 24, 25, 26 };

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
    
    Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
            Hi(i,n) = m_H(ni[i],n);
    Hi = Hi.inverse();
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FEHex20G27::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += Hi[j][k]*ai[ni[k]];
        }
    }
}

//=============================================================================
//              H E X 2 7
//=============================================================================

void FEHex27_::init()
{
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_HEX8));
	m_shapeP[2] = dynamic_cast<FESolidElementShape*>(FEElementLibrary::GetElementShapeClass(ET_HEX27));

	// initialize base class
	FESolidElementTraits::init();
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
	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
//! \todo implement this
void FEHex27G27::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NINT; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[k];
        }
    }
}

//=============================================================================
//              P Y R A 5
//=============================================================================

//=============================================================================
//              P Y R A 5 G 8
//=============================================================================

FEPyra5G8::FEPyra5G8() : FEPyra5_(NINT, FE_PYRA5G8)
{
	// integration point coordinates
	const double a = 1.0 / sqrt(3.0);
	gr[0] = -a; gs[0] = -a; gt[0] = -a; gw[0] = 1;
	gr[1] = a; gs[1] = -a; gt[1] = -a; gw[1] = 1;
	gr[2] = a; gs[2] = a; gt[2] = -a; gw[2] = 1;
	gr[3] = -a; gs[3] = a; gt[3] = -a; gw[3] = 1;
	gr[4] = -a; gs[4] = -a; gt[4] = a; gw[4] = 1;
	gr[5] = a; gs[5] = -a; gt[5] = a; gw[5] = 1;
	gr[6] = a; gs[6] = a; gt[6] = a; gw[6] = 1;
	gr[7] = -a; gs[7] = a; gt[7] = a; gw[7] = 1;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN, NELN);
	m_Ai.resize(NELN, NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

void FEPyra5G8::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i = 0; i<NELN; ++i)
	{
		b[i] = 0;
		for (int j = 0; j<NINT; ++j) b[i] += m_H[j][i] * ai[j];
	}

	for (int i = 0; i<NELN; ++i)
	{
		ao[i] = 0;
		for (int j = 0; j<NELN; ++j) ao[i] += m_Ai[i][j] * b[j];
	}
}


//=============================================================================
//              P Y R A 1 3
//=============================================================================

//=============================================================================
//              P Y R A 1 3 G 8
//=============================================================================

FEPyra13G8::FEPyra13G8() : FEPyra13_(NINT, FE_PYRA13G8)
{
    // integration point coordinates
    const double a = 1.0 / sqrt(3.0);
    gr[0] = -a; gs[0] = -a; gt[0] = -a; gw[0] = 1;
    gr[1] = a; gs[1] = -a; gt[1] = -a; gw[1] = 1;
    gr[2] = a; gs[2] = a; gt[2] = -a; gw[2] = 1;
    gr[3] = -a; gs[3] = a; gt[3] = -a; gw[3] = 1;
    gr[4] = -a; gs[4] = -a; gt[4] = a; gw[4] = 1;
    gr[5] = a; gs[5] = -a; gt[5] = a; gw[5] = 1;
    gr[6] = a; gs[6] = a; gt[6] = a; gw[6] = 1;
    gr[7] = -a; gs[7] = a; gt[7] = a; gw[7] = 1;
    init();
    
    // we need Ai to project integration point data to the nodes
    matrix A(NELN, NELN);
    m_Ai.resize(NELN, NELN);
    A = m_H.transpose()*m_H;
    m_Ai = A.inverse();
}

void FEPyra13G8::project_to_nodes(double* ai, double* ao) const
{
    vector<double> b(NELN);
    for (int i = 0; i<NELN; ++i)
    {
        b[i] = 0;
        for (int j = 0; j<NINT; ++j) b[i] += m_H[j][i] * ai[j];
    }
    
    for (int i = 0; i<NELN; ++i)
    {
        ao[i] = 0;
        for (int j = 0; j<NELN; ++j) ao[i] += m_Ai[i][j] * b[j];
    }
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
	assert(m_nint > 0);
	assert(m_neln > 0);

	// get shape class
	m_shape = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(m_spec.eshape));
	assert(m_shape && (m_shape->shape() == m_spec.eshape));

	// evaluate shape functions
	const int NE = FEElement::MAX_NODES;
	double N[NE];
	for (int n=0; n<m_nint; ++n)
	{
		shape_fnc(N, gr[n], gs[n]);
		for (int i=0; i<m_neln; ++i) m_H[n][i] = N[i];
	}
	
	// evaluate shape function derivatives
	double Nr[NE], Ns[NE];
	for (int n=0; n<m_nint; ++n)
	{
		shape_deriv(Nr, Ns, gr[n], gs[n]);
		for (int i=0; i<m_neln; ++i)
		{
			Gr[n][i] = Nr[i];
			Gs[n][i] = Ns[i];
		}
	}

	// NOTE: Below, is a new interface for dealing with mixed element formulations.
	//       This is still a work in progress.

	// Get the max interpolation order
	const int maxOrder = (int)m_shapeP.size() - 1;
	m_Hp.resize(maxOrder + 1);
	Gr_p.resize(maxOrder + 1);
	Gs_p.resize(maxOrder + 1);
	for (int i = 0; i <= maxOrder; ++i)
	{
		FESurfaceElementShape* shape = m_shapeP[i];
		matrix& H = m_Hp[i];
		matrix& Gr = Gr_p[i];
		matrix& Gs = Gs_p[i];
		if (i == 0)
		{
			H.resize(m_nint, 1);
			Gr.resize(m_nint, 1);
			Gs.resize(m_nint, 1);
			for (int n = 0; n < m_nint; ++n)
			{
				H[n][0] = 1.0;
				Gr[n][0] = Gs[n][0] = 0.0;
			}
		}
		else if (m_shapeP[i])
		{
			// get the nodes
			int neln = shape->nodes();

			// shape function values
			H.resize(m_nint, neln);
			for (int n = 0; n<m_nint; ++n)
			{
				m_shapeP[i]->shape_fnc(N, gr[n], gs[n]);
				for (int j = 0; j<neln; ++j) H[n][j] = N[j];
			}

			// calculate local derivatives of shape functions at gauss points
			Gr.resize(m_nint, neln);
			Gs.resize(m_nint, neln);
			for (int n = 0; n<m_nint; ++n)
			{
				shape->shape_deriv(Nr, Ns, gr[n], gs[n]);
				for (int j = 0; j<neln; ++j)
				{
					Gr[n][j] = Nr[j];
					Gs[n][j] = Ns[j];
				}
			}
		}
	}
}

// shape functions at (r,s)
void FESurfaceElementTraits::shape_fnc(double* H, double r, double s)
{
	m_shape->shape_fnc(H, r, s);
}

// shape function derivatives at (r,s)
void FESurfaceElementTraits::shape_deriv(double* Gr, double* Gs, double r, double s)
{
	m_shape->shape_deriv(Gr, Gs, r, s);
}

// shape function derivatives at (r,s)
void FESurfaceElementTraits::shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
{
	m_shape->shape_deriv2(Grr, Gss, Grs, r, s);
}

// shape functions at (r,s)
void FESurfaceElementTraits::shape_fnc(int order, double* H, double r, double s)
{
	if (order == -1) m_shape->shape_fnc(H, r, s);
	else m_shapeP[order]->shape_fnc(H, r, s);
}

// shape function derivatives at (r,s)
void FESurfaceElementTraits::shape_deriv(int order, double* Gr, double* Gs, double r, double s)
{
	if (order == -1) shape_deriv(Gr, Gs, r, s);
	else m_shapeP[order]->shape_deriv(Gr, Gs, r, s);
}

//=============================================================================
//                          F E Q U A D 4
//=============================================================================

void FEQuad4_::init()
{
	// allocate shape classes
	m_shapeP.resize(2);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_QUAD4));

    // centroid coordinates
    cr = cs = 0;

    // initialize base class
	FESurfaceElementTraits::init();
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
	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void FEQuad4G4::project_to_nodes(double* ai, double* ao) const
{
	int ni = NINT;
	int ne = NELN;
	assert(ni == ne);
	for (int i=0; i<ne; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<ni; ++j) ao[i] += m_Hi[i][j]*ai[j];
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
void FEQuad4NI::project_to_nodes(double* ai, double* ao) const
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
	ao[3] = ai[3];
}

//=============================================================================
//                          F E T R I 
//=============================================================================

void FETri3_::init()
{
	// allocate shape classes
	m_shapeP.resize(2);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI3));

    // centroid coordinates
    cr = cs = 1.0/3.0;
    
	// initialize base class
	FESurfaceElementTraits::init();
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
void FETri3G1::project_to_nodes(double* ai, double* ao) const
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
	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void FETri3G3::project_to_nodes(double* ai, double* ao) const
{
	assert(NINT == NELN);
	for (int i=0; i<NELN; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<NINT; ++j) ao[i] += m_Hi[i][j]*ai[j];
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri3G7::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
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
void FETri3NI::project_to_nodes(double* ai, double* ao) const
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
}

//============================================================================
//                             F E T R I 6
//============================================================================

void FETri6_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI3));
	m_shapeP[2] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI6));

    // centroid coordinates
    cr = cs = 1.0/3.0;
    
	// initialize base class
	FESurfaceElementTraits::init();
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
void FETri6G3::project_to_nodes(double* ai, double* ao) const
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
void FETri6G4::project_to_nodes(double* ai, double* ao) const
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri6G7::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
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
void FETri6GL7::project_to_nodes(double* ai, double* ao) const
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
void FETri6NI::project_to_nodes(double* ai, double* ao) const
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
/*
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri6mG7::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
	}
}
*/

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
void FETri7G3::project_to_nodes(double* ai, double* ao) const
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
void FETri7G4::project_to_nodes(double* ai, double* ao) const
{
	
}

//============================================================================
//                             F E T R I 7
//============================================================================

void FETri7_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI3));
	m_shapeP[2] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI7));

    // centroid coordinates
    cr = cs = 1.0/3.0;
    
	// initialize base class
	FESurfaceElementTraits::init();
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri7G7::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
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
void FETri7GL7::project_to_nodes(double* ai, double* ao) const
{
	ao[0] = ai[0]; ao[1] = ai[1]; ao[2] = ai[2];
	ao[3] = ai[3]; ao[4] = ai[4]; ao[5] = ai[5];
	ao[6] = ai[6];
}

//============================================================================
//                             F E T R I 1 0
//============================================================================

void FETri10_::init()
{
	// allocate shape classes
	m_shapeP.resize(4);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI3));
	m_shapeP[2] = 0; // this element cannot be used for quadratic interpolation!
	m_shapeP[3] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_TRI10));

    // centroid coordinates
    cr = cs = 1.0/3.0;
    
	// initialize base class
	FESurfaceElementTraits::init();
}

//=============================================================================
//                          F E T R I 1 0 G 7
//=============================================================================

FETri10G7::FETri10G7() : FETri10_(NINT, FE_TRI10G7)
{
	const double w = 1.0 / 2.0;
	gr[0] = 0.333333333333333; gs[0] = 0.333333333333333; gw[0] = w*0.225000000000000;
	gr[1] = 0.797426985353087; gs[1] = 0.101286507323456; gw[1] = w*0.125939180544827;
	gr[2] = 0.101286507323456; gs[2] = 0.797426985353087; gw[2] = w*0.125939180544827;
	gr[3] = 0.101286507323456; gs[3] = 0.101286507323456; gw[3] = w*0.125939180544827;
	gr[4] = 0.470142064105115; gs[4] = 0.470142064105115; gw[4] = w*0.132394152788506;
	gr[5] = 0.470142064105115; gs[5] = 0.059715871789770; gw[5] = w*0.132394152788506;
	gr[6] = 0.059715871789770; gs[6] = 0.470142064105115; gw[6] = w*0.132394152788506;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN, NELN);
	m_Ai.resize(NELN, NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri10G7::project_to_nodes(double* ai, double* ao) const
{
	// TODO: Implement this
}


//=============================================================================
//                          F E T R I 1 0 G 1 2
//=============================================================================

FETri10G12::FETri10G12() : FETri10_(NINT, FE_TRI10G12)
{
	gr[ 0] = 0.063089014; gs[ 0] = 0.873821971;	gw[ 0] = 0.025422453;
	gr[ 1] = 0.873821971; gs[ 1] = 0.063089014;	gw[ 1] = 0.025422453;
	gr[ 2] = 0.063089014; gs[ 2] = 0.063089014;	gw[ 2] = 0.025422453;
	gr[ 3] = 0.249286745; gs[ 3] = 0.501426510;	gw[ 3] = 0.058393138;
	gr[ 4] = 0.501426510; gs[ 4] = 0.249286745;	gw[ 4] = 0.058393138;
	gr[ 5] = 0.249286745; gs[ 5] = 0.249286745;	gw[ 5] = 0.058393138;
	gr[ 6] = 0.053145050; gs[ 6] = 0.636502499;	gw[ 6] = 0.041425538;
	gr[ 7] = 0.636502499; gs[ 7] = 0.053145050;	gw[ 7] = 0.041425538;
	gr[ 8] = 0.310352451; gs[ 8] = 0.636502499;	gw[ 8] = 0.041425538;
	gr[ 9] = 0.636502499; gs[ 9] = 0.310352451;	gw[ 9] = 0.041425538;
	gr[10] = 0.310352451; gs[10] = 0.053145050;	gw[10] = 0.041425538;
	gr[11] = 0.053145050; gs[11] = 0.310352451;	gw[11] = 0.041425538;
	init();

	// we need Ai to project integration point data to the nodes
	matrix A(NELN, NELN);
	m_Ai.resize(NELN, NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FETri10G12::project_to_nodes(double* ai, double* ao) const
{
	// TODO: Implement this
}

//=============================================================================
//          F E Q U A D 8 
//=============================================================================

void FEQuad8_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_QUAD4));
	m_shapeP[2] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_QUAD8));

    // centroid coordinates
    cr = cs = 0;
    
	// initialize base class
	FESurfaceElementTraits::init();
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FEQuad8G9::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
	}
}


//=============================================================================
//                          F E Q U A D 8 N I
//=============================================================================

FEQuad8NI::FEQuad8NI() : FEQuad8_(NINT, FE_QUAD8NI)
{
    double w = 1./9.;
    gr[0] = -1; gs[0] = -1; gw[0] = w;
    gr[1] =  1; gs[1] = -1; gw[1] = w;
    gr[2] =  1; gs[2] =  1; gw[2] = w;
    gr[3] = -1; gs[3] =  1; gw[3] = w;
    gr[4] =  0; gs[4] = -1; gw[0] = 4*w;
    gr[5] =  1; gs[5] =  0; gw[1] = 4*w;
    gr[6] =  0; gs[6] =  1; gw[2] = 4*w;
    gr[7] = -1; gs[7] =  0; gw[3] = 4*w;
    init();
}

//-----------------------------------------------------------------------------
void FEQuad8NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
    ao[3] = ai[3];
    ao[4] = ai[4];
    ao[5] = ai[5];
    ao[6] = ai[6];
    ao[7] = ai[7];
}

//=============================================================================
//          F E Q U A D 9 
//=============================================================================

void FEQuad9_::init()
{
	// allocate shape classes
	m_shapeP.resize(3);
	m_shapeP[0] = 0;
	m_shapeP[1] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_QUAD4));
	m_shapeP[2] = dynamic_cast<FESurfaceElementShape*>(FEElementLibrary::GetElementShapeClass(ET_QUAD9));

    // centroid coordinates
    cr = cs = 0;
    
	// initialize base class
	FESurfaceElementTraits::init();
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FEQuad9G9::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
	}
}

//=============================================================================
//                          F E Q U A D 9 N I
//=============================================================================

FEQuad9NI::FEQuad9NI() : FEQuad9_(NINT, FE_QUAD9NI)
{
    double w = 1./9.;
    gr[0] = -1; gs[0] = -1; gw[0] = w;
    gr[1] =  1; gs[1] = -1; gw[1] = w;
    gr[2] =  1; gs[2] =  1; gw[2] = w;
    gr[3] = -1; gs[3] =  1; gw[3] = w;
    gr[4] =  0; gs[4] = -1; gw[0] = 4*w;
    gr[5] =  1; gs[5] =  0; gw[1] = 4*w;
    gr[6] =  0; gs[6] =  1; gw[2] = 4*w;
    gr[7] = -1; gs[7] =  0; gw[3] = 4*w;
    gr[8] =  0; gs[8] =  0; gw[3] = 16*w;
    init();
}

//-----------------------------------------------------------------------------
void FEQuad9NI::project_to_nodes(double* ai, double* ao) const
{
    ao[0] = ai[0];
    ao[1] = ai[1];
    ao[2] = ai[2];
    ao[3] = ai[3];
    ao[4] = ai[4];
    ao[5] = ai[5];
    ao[6] = ai[6];
    ao[7] = ai[7];
    ao[8] = ai[8];
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

	m_faces = 1;
}

//-----------------------------------------------------------------------------
//! initialize element traits data
void FEShellElementTraits::init()
{
    assert(m_nint > 0);
    assert(m_neln > 0);
    const int NELN = FEElement::MAX_NODES;
    
    // calculate shape function values at gauss points
    double N[NELN];
    for (int n=0; n<m_nint; ++n)
    {
        shape_fnc(N, gr[n], gs[n]);
        for (int i=0; i<m_neln; ++i) m_H[n][i] = N[i];
    }
    
    // calculate local derivatives of shape functions at gauss points
    double Nr[NELN], Ns[NELN];
    for (int n=0; n<m_nint; ++n)
    {
        shape_deriv(Nr, Ns, gr[n], gs[n]);
        for (int i=0; i<m_neln; ++i)
        {
            Hr[n][i] = Nr[i];
            Hs[n][i] = Ns[i];
        }
    }
}

//=============================================================================
//                                F E S H E L L Q U A D 4
//=============================================================================

//-----------------------------------------------------------------------------
void FEShellQuad4_::shape_fnc(double* H, double r, double s)
{
    H[0] = 0.25*(1-r)*(1-s);
    H[1] = 0.25*(1+r)*(1-s);
    H[2] = 0.25*(1+r)*(1+s);
    H[3] = 0.25*(1-r)*(1+s);
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FEShellQuad4_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
    Hr[0] = -0.25*(1-s);
    Hr[1] =  0.25*(1-s);
    Hr[2] =  0.25*(1+s);
    Hr[3] = -0.25*(1+s);
    
    Hs[0] = -0.25*(1-r);
    Hs[1] = -0.25*(1+r);
    Hs[2] =  0.25*(1+r);
    Hs[3] =  0.25*(1-r);
}

//*****************************************************************************
//                          S H E L L Q U A D 4 G 8
//*****************************************************************************

int FEShellQuad4G8::ni[NELN] = { 4, 5, 6, 7 };

FEShellQuad4G8::FEShellQuad4G8() : FEShellQuad4_(NINT, FE_SHELL_QUAD4G8)
{
    const double a = 1.0 / sqrt(3.0);
    const double w = 1.0;
    
    gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -a; gw[ 0] = w;
    gr[ 1] =  a; gs[ 1] = -a; gt[ 1] = -a; gw[ 1] = w;
    gr[ 2] =  a; gs[ 2] =  a; gt[ 2] = -a; gw[ 2] = w;
    gr[ 3] = -a; gs[ 3] =  a; gt[ 3] = -a; gw[ 3] = w;
    
    gr[ 4] = -a; gs[ 4] = -a; gt[ 4] =  a; gw[ 4] = w;
    gr[ 5] =  a; gs[ 5] = -a; gt[ 5] =  a; gw[ 5] = w;
    gr[ 6] =  a; gs[ 6] =  a; gt[ 6] =  a; gw[ 6] = w;
    gr[ 7] = -a; gs[ 7] =  a; gt[ 7] =  a; gw[ 7] = w;
    
    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
    m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellQuad4G8::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//*****************************************************************************
//                          S H E L L Q U A D 4 G 1 2
//*****************************************************************************

int FEShellQuad4G12::ni[NELN] = { 8, 9, 10, 11 };

FEShellQuad4G12::FEShellQuad4G12() : FEShellQuad4_(NINT, FE_SHELL_QUAD4G12)
{
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
    
    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellQuad4G12::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//=============================================================================
//                                F E S H E L L T R I 3
//=============================================================================

//-----------------------------------------------------------------------------
void FEShellTri3_::shape_fnc(double* H, double r, double s)
{
    H[0] = 1-r-s;
    H[1] = r;
    H[2] = s;
}

//-----------------------------------------------------------------------------
//! values of shape function derivatives
void FEShellTri3_::shape_deriv(double* Hr, double* Hs, double r, double s)
{
    Hr[0] = -1;
    Hr[1] =  1;
    Hr[2] =  0;
    
    Hs[0] = -1;
    Hs[1] =  0;
    Hs[2] =  1;
}

//*****************************************************************************
//                          S H E L L T R I 3 G 6
//*****************************************************************************

int FEShellTri3G6::ni[NELN] = { 3, 4, 5 };

FEShellTri3G6::FEShellTri3G6() : FEShellTri3_(NINT, FE_SHELL_TRI3G6)
{
    //gauss intergration points
    const double a = 1.0/6.0;
    const double b = 2.0/3.0;
    const double c = 1.0 / sqrt(3.0);
    
    gr[0] = a; gs[0] = a; gt[0] = -c; gw[0] = a;
    gr[1] = b; gs[1] = a; gt[1] = -c; gw[1] = a;
    gr[2] = a; gs[2] = b; gt[2] = -c; gw[2] = a;
    gr[3] = a; gs[3] = a; gt[3] =  c; gw[3] = a;
    gr[4] = b; gs[4] = a; gt[4] =  c; gw[4] = a;
    gr[5] = a; gs[5] = b; gt[5] =  c; gw[5] = a;

    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellTri3G6::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//*****************************************************************************
//                          S H E L L T R I 3 G 9
//*****************************************************************************

int FEShellTri3G9::ni[NELN] = { 6, 7, 8 };

FEShellTri3G9::FEShellTri3G9() : FEShellTri3_(NINT, FE_SHELL_TRI3G9)
{
    const double a = 1.0 / 6.0;
    const double b = 2.0 / 3.0;
    const double w1 = 5.0 / 9.0;
    const double w2 = 8.0 / 9.0;
    
    gr[0] = a; gs[0] = a; gt[0] = -b; gw[0] = a*w1;
    gr[1] = b; gs[1] = a; gt[1] = -b; gw[1] = a*w1;
    gr[2] = a; gs[2] = b; gt[2] = -b; gw[2] = a*w1;
    
    gr[3] = a; gs[3] = a; gt[3] =  0; gw[3] = a*w2;
    gr[4] = b; gs[4] = a; gt[4] =  0; gw[4] = a*w2;
    gr[5] = a; gs[5] = b; gt[5] =  0; gw[5] = a*w2;
    
    gr[6] = a; gs[6] = a; gt[6] =  b; gw[6] = a*w1;
    gr[7] = b; gs[7] = a; gt[7] =  b; gw[7] = a*w1;
    gr[8] = a; gs[8] = b; gt[8] =  b; gw[8] = a*w1;
    
    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellTri3G9::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//=============================================================================
//                                F E S H E L L Q U A D 8
//=============================================================================

//-----------------------------------------------------------------------------
void FEShellQuad8_::shape_fnc(double* H, double r, double s)
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
//! values of shape function derivatives
void FEShellQuad8_::shape_deriv(double* Hr, double* Hs, double r, double s)
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

//*****************************************************************************
//                          S H E L L Q U A D 8 G 1 8
//*****************************************************************************

int FEShellQuad8G18::ni[NELN] = { 9, 10, 11, 12, 14, 15, 16, 17 };

FEShellQuad8G18::FEShellQuad8G18() : FEShellQuad8_(NINT, FE_SHELL_QUAD8G18)
{
    // integration point coordinates
    const double a = 0.774596669241483;
    const double c = 0.577350269189626;
    const double w1 = 5.0 / 9.0;
    const double w2 = 8.0 / 9.0;
    gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -c; gw[ 0] = w1*w1;
    gr[ 1] =  0; gs[ 1] = -a; gt[ 1] = -c; gw[ 1] = w2*w1;
    gr[ 2] =  a; gs[ 2] = -a; gt[ 2] = -c; gw[ 2] = w1*w1;
    gr[ 3] = -a; gs[ 3] =  0; gt[ 3] = -c; gw[ 3] = w1*w2;
    gr[ 4] =  0; gs[ 4] =  0; gt[ 4] = -c; gw[ 4] = w2*w2;
    gr[ 5] =  a; gs[ 5] =  0; gt[ 5] = -c; gw[ 5] = w1*w2;
    gr[ 6] = -a; gs[ 6] =  a; gt[ 6] = -c; gw[ 6] = w1*w1;
    gr[ 7] =  0; gs[ 7] =  a; gt[ 7] = -c; gw[ 7] = w2*w1;
    gr[ 8] =  a; gs[ 8] =  a; gt[ 8] = -c; gw[ 8] = w1*w1;
    gr[ 9] = -a; gs[ 9] = -a; gt[ 9] =  c; gw[ 9] = w1*w1;
    gr[10] =  0; gs[10] = -a; gt[10] =  c; gw[10] = w2*w1;
    gr[11] =  a; gs[11] = -a; gt[11] =  c; gw[11] = w1*w1;
    gr[12] = -a; gs[12] =  0; gt[12] =  c; gw[12] = w1*w2;
    gr[13] =  0; gs[13] =  0; gt[13] =  c; gw[13] = w2*w2;
    gr[14] =  a; gs[14] =  0; gt[14] =  c; gw[14] = w1*w2;
    gr[15] = -a; gs[15] =  a; gt[15] =  c; gw[15] = w1*w1;
    gr[16] =  0; gs[16] =  a; gt[16] =  c; gw[16] = w2*w1;
    gr[17] =  a; gs[17] =  a; gt[17] =  c; gw[17] = w1*w1;
    
    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellQuad8G18::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//*****************************************************************************
//                          S H E L L Q U A D 8 G 2 7
//*****************************************************************************

int FEShellQuad8G27::ni[NELN] = { 18, 19, 20, 21, 23, 24, 25, 26 };

FEShellQuad8G27::FEShellQuad8G27() : FEShellQuad8_(NINT, FE_SHELL_QUAD8G27)
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
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellQuad8G27::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//=============================================================================
//                                F E S H E L L T R I 6
//=============================================================================

//-----------------------------------------------------------------------------
void FEShellTri6_::shape_fnc(double* H, double r, double s)
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
//! values of shape function derivatives
void FEShellTri6_::shape_deriv(double* Hr, double* Hs, double r, double s)
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

//*****************************************************************************
//                          S H E L L T R I 6 G 1 4
//*****************************************************************************

int FEShellTri6G14::ni[NELN] = { 8, 9, 10, 11, 12, 13 };

FEShellTri6G14::FEShellTri6G14() : FEShellTri6_(NINT, FE_SHELL_TRI6G14)
{
    const double a = 0.774596669241483;
    const double c = 0.577350269189626;
    const double w = 1.0/2.0;
    
    gr[ 0] = 0.333333333333333; gs[ 0] = 0.333333333333333; gt[ 0] = -c; gw[ 0] = w*0.225000000000000;
    gr[ 1] = 0.797426985353087; gs[ 1] = 0.101286507323456; gt[ 1] = -c; gw[ 1] = w*0.125939180544827;
    gr[ 2] = 0.101286507323456; gs[ 2] = 0.797426985353087; gt[ 2] = -c; gw[ 2] = w*0.125939180544827;
    gr[ 3] = 0.101286507323456; gs[ 3] = 0.101286507323456; gt[ 3] = -c; gw[ 3] = w*0.125939180544827;
    gr[ 4] = 0.470142064105115; gs[ 4] = 0.470142064105115; gt[ 4] = -c; gw[ 4] = w*0.132394152788506;
    gr[ 5] = 0.470142064105115; gs[ 5] = 0.059715871789770; gt[ 5] = -c; gw[ 5] = w*0.132394152788506;
    gr[ 6] = 0.059715871789770; gs[ 6] = 0.470142064105115; gt[ 6] = -c; gw[ 6] = w*0.132394152788506;
    
    gr[ 7] = 0.333333333333333; gs[ 7] = 0.333333333333333; gt[ 7] =  c; gw[ 7] = w*0.225000000000000;
    gr[ 8] = 0.797426985353087; gs[ 8] = 0.101286507323456; gt[ 8] =  c; gw[ 8] = w*0.125939180544827;
    gr[ 9] = 0.101286507323456; gs[ 9] = 0.797426985353087; gt[ 9] =  c; gw[ 9] = w*0.125939180544827;
    gr[10] = 0.101286507323456; gs[10] = 0.101286507323456; gt[10] =  c; gw[10] = w*0.125939180544827;
    gr[11] = 0.470142064105115; gs[11] = 0.470142064105115; gt[11] =  c; gw[11] = w*0.132394152788506;
    gr[12] = 0.470142064105115; gs[12] = 0.059715871789770; gt[12] =  c; gw[12] = w*0.132394152788506;
    gr[13] = 0.059715871789770; gs[13] = 0.470142064105115; gt[13] =  c; gw[13] = w*0.132394152788506;
    
    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellTri6G14::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
    }
}

//*****************************************************************************
//                          S H E L L T R I 6 G 2 1
//*****************************************************************************

int FEShellTri6G21::ni[NELN] = { 15, 16, 17, 18, 19, 20 };

FEShellTri6G21::FEShellTri6G21() : FEShellTri6_(NINT, FE_SHELL_TRI6G21)
{
    const double a = 0.774596669241483;
    const double w = 1.0/2.0;
    const double w1 = 5.0 / 9.0;
    const double w2 = 8.0 / 9.0;
    
    gr[ 0] = 0.333333333333333; gs[ 0] = 0.333333333333333; gt[ 0] = -a; gw[ 0] = w*w1*0.225000000000000;
    gr[ 1] = 0.797426985353087; gs[ 1] = 0.101286507323456; gt[ 1] = -a; gw[ 1] = w*w1*0.125939180544827;
    gr[ 2] = 0.101286507323456; gs[ 2] = 0.797426985353087; gt[ 2] = -a; gw[ 2] = w*w1*0.125939180544827;
    gr[ 3] = 0.101286507323456; gs[ 3] = 0.101286507323456; gt[ 3] = -a; gw[ 3] = w*w1*0.125939180544827;
    gr[ 4] = 0.470142064105115; gs[ 4] = 0.470142064105115; gt[ 4] = -a; gw[ 4] = w*w1*0.132394152788506;
    gr[ 5] = 0.470142064105115; gs[ 5] = 0.059715871789770; gt[ 5] = -a; gw[ 5] = w*w1*0.132394152788506;
    gr[ 6] = 0.059715871789770; gs[ 6] = 0.470142064105115; gt[ 6] = -a; gw[ 6] = w*w1*0.132394152788506;
    
    gr[ 7] = 0.333333333333333; gs[ 7] = 0.333333333333333; gt[ 7] =  0; gw[ 7] = w*w2*0.225000000000000;
    gr[ 8] = 0.797426985353087; gs[ 8] = 0.101286507323456; gt[ 8] =  0; gw[ 8] = w*w2*0.125939180544827;
    gr[ 9] = 0.101286507323456; gs[ 9] = 0.797426985353087; gt[ 9] =  0; gw[ 9] = w*w2*0.125939180544827;
    gr[10] = 0.101286507323456; gs[10] = 0.101286507323456; gt[10] =  0; gw[10] = w*w2*0.125939180544827;
    gr[11] = 0.470142064105115; gs[11] = 0.470142064105115; gt[11] =  0; gw[11] = w*w2*0.132394152788506;
    gr[12] = 0.470142064105115; gs[12] = 0.059715871789770; gt[12] =  0; gw[12] = w*w2*0.132394152788506;
    gr[13] = 0.059715871789770; gs[13] = 0.470142064105115; gt[13] =  0; gw[13] = w*w2*0.132394152788506;
    
    gr[14] = 0.333333333333333; gs[14] = 0.333333333333333; gt[14] =  a; gw[14] = w*w1*0.225000000000000;
    gr[15] = 0.797426985353087; gs[15] = 0.101286507323456; gt[15] =  a; gw[15] = w*w1*0.125939180544827;
    gr[16] = 0.101286507323456; gs[16] = 0.797426985353087; gt[16] =  a; gw[16] = w*w1*0.125939180544827;
    gr[17] = 0.101286507323456; gs[17] = 0.101286507323456; gt[17] =  a; gw[17] = w*w1*0.125939180544827;
    gr[18] = 0.470142064105115; gs[18] = 0.470142064105115; gt[18] =  a; gw[18] = w*w1*0.132394152788506;
    gr[19] = 0.470142064105115; gs[19] = 0.059715871789770; gt[19] =  a; gw[19] = w*w1*0.132394152788506;
    gr[20] = 0.059715871789770; gs[20] = 0.470142064105115; gt[20] =  a; gw[20] = w*w1*0.132394152788506;
    
    init();
    
	m_Hi.resize(NELN, NELN);
    for (int i=0; i<NELN; ++i)
        for (int n=0; n<NELN; ++n)
			m_Hi(i,n) = m_H(ni[i],n);
	m_Hi = m_Hi.inverse();
}

//-----------------------------------------------------------------------------
//! project to nodes
void FEShellTri6G21::project_to_nodes(double* ai, double* ao) const
{
    for (int j=0; j<NELN; ++j)
    {
        ao[j] = 0;
        for (int k=0; k<NELN; ++k)
        {
            ao[j] += m_Hi[j][k]*ai[ni[k]];
        }
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
	assert(m_nint > 0);
	assert(m_neln > 0);

	// evaluate shape functions
	const int NE = FEElement::MAX_NODES;
	double N[NE];
	for (int n=0; n<m_nint; ++n)
	{
		shape(N, gr[n], gs[n]);
		for (int i=0; i<m_neln; ++i) m_H[n][i] = N[i];
	}
	
	// evaluate shape function derivatives
	double Nr[NE], Ns[NE];
	for (int n=0; n<m_nint; ++n)
	{
		shape_deriv(Nr, Ns, gr[n], gs[n]);
		for (int i=0; i<m_neln; ++i)
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
void FE2DTri3G1::project_to_nodes(double* ai, double* ao) const
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
void FE2DTri6G3::project_to_nodes(double* ai, double* ao) const
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
	m_Hi = m_H.inverse();
}

//-----------------------------------------------------------------------------
void FE2DQuad4G4::project_to_nodes(double* ai, double* ao) const
{
	int ni = NINT;
	int ne = NELN;
	assert(ni == ne);
	for (int i=0; i<ne; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<ni; ++j) ao[i] += m_Hi[i][j]*ai[j];
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FE2DQuad8G9::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
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
	m_Ai.resize(NELN,NELN);
	A = m_H.transpose()*m_H;
	m_Ai = A.inverse();
}

//-----------------------------------------------------------------------------
void FE2DQuad9G9::project_to_nodes(double* ai, double* ao) const
{
	vector<double> b(NELN);
	for (int i=0; i<NELN; ++i) 
	{
		b[i] = 0;
		for (int j=0; j<NINT; ++j) b[i] += m_H[j][i]*ai[j];
	}

	for (int i=0; i<NELN; ++i) 
	{
		ao[i] = 0;
		for (int j=0; j<NELN; ++j) ao[i] += m_Ai[i][j]*b[j];
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
	assert(m_nint > 0);
	assert(m_neln > 0);

	// evaluate shape functions
	const int NE = FEElement::MAX_NODES;
	double N[NE];
	for (int n=0; n<m_nint; ++n)
	{
		shape(N, gr[n]);
		for (int i=0; i<m_neln; ++i) m_H[n][i] = N[i];
	}
	
	// evaluate shape function derivatives
	double Nr[NE];
	for (int n=0; n<m_nint; ++n)
	{
		shape_deriv(Nr, gr[n]);
		for (int i=0; i<m_neln; ++i)
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
void FELine2G1::project_to_nodes(double* ai, double* ao) const
{
	ao[0] = ai[0];
	ao[1] = ai[0];
}
