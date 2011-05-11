// FEElement.cpp: implementation of the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElement.h"
#include <math.h>

//-----------------------------------------------------------------------------
FEElementState::FEElementState(const FEElementState& s)
{
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) m_data[i] = s.m_data[i]->Copy();
}

FEElementState& FEElementState::operator = (const FEElementState& s)
{
	Clear();
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) m_data[i] = s.m_data[i]->Copy();
	return (*this);
}

//-----------------------------------------------------------------------------
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy solid element data
	m_eJ = el.m_eJ;
	m_ep = el.m_ep;
	m_Lk = el.m_Lk;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy solid element data
	m_eJ = el.m_eJ;
	m_ep = el.m_ep;
	m_Lk = el.m_Lk;

	return (*this);
}

double FESolidElement::defgrad(mat3d& F, int n)
{
	// make sure this element is unpacked
	assert(m_pT->m_pel == this);

	double *Grn = Gr(n);
	double *Gsn = Gs(n);
	double *Gtn = Gt(n);

	vec3d* r = rt();

	double GX, GY, GZ;
	double x, y, z;
	double Gri, Gsi, Gti;
	double Ji[3][3];
	invjac0(Ji, n);

	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = Nodes();
	for (int i=0; i<neln; ++i)
	{
		Gri = Grn[i];
		Gsi = Gsn[i];
		Gti = Gtn[i];

		x = r[i].x;
		y = r[i].y;
		z = r[i].z;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
		GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
		GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
	
		// calculate deformation gradient F
		F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
		F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
		F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
	}

	double D = F.det();
	if (D <= 0) throw NegativeJacobian(m_nID, n, D, this);

	return D;
}

//-----------------------------------------------------------------------------
FEShellElement::FEShellElement(const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy shell data
	m_eJ = el.m_eJ;
	m_ep = el.m_ep;
	m_Lk = el.m_Lk;
	m_h0 = el.m_h0;
}

//! assignment operator
FEShellElement& FEShellElement::operator = (const FEShellElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy shell data
	m_eJ = el.m_eJ;
	m_ep = el.m_ep;
	m_Lk = el.m_Lk;
	m_h0 = el.m_h0;

	return (*this);
}

//-----------------------------------------------------------------------------
FETrussElement::FETrussElement(const FETrussElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// truss data
	m_a0 = el.m_a0;
}

FETrussElement& FETrussElement::operator = (const FETrussElement& el) 
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy truss data
	m_a0 = el.m_a0;

	return (*this); 
}

//-----------------------------------------------------------------------------
FEDiscreteElement::FEDiscreteElement(const FEDiscreteElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;
}

FEDiscreteElement& FEDiscreteElement::operator =(const FEDiscreteElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	return (*this);
}
