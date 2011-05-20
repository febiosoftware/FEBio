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

double FEShellElement::defgrad(mat3d& F, int n)
{
	// make sure this element is unpacked
	assert(m_pT->m_pel == this);

	double* Hrn = Hr(n);
	double* Hsn = Hs(n);
	double* Hn  = H(n);
	double NX, NY, NZ, MX, MY, MZ;
	double za;

	vec3d* r = rt();
	vec3d* D = Dt();

	double g = gt(n);

	double Ji[3][3];
	invjac0(Ji, n);

	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = Nodes();
	for (int i=0; i<neln; ++i)
	{
		const double& Hri = Hrn[i];
		const double& Hsi = Hsn[i];
		const double& Hi  = Hn[i];

		const double& x = r[i].x;
		const double& y = r[i].y;
		const double& z = r[i].z;

		const double& dx = D[i].x;
		const double& dy = D[i].y;
		const double& dz = D[i].z;

		za = 0.5*g*m_h0[i];

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		NX = Ji[0][0]*Hri+Ji[1][0]*Hsi;
		NY = Ji[0][1]*Hri+Ji[1][1]*Hsi;
		NZ = Ji[0][2]*Hri+Ji[1][2]*Hsi;

		MX = za*Ji[0][0]*Hri + za*Ji[1][0]*Hsi + Ji[2][0]*0.5*m_h0[i]*Hi;
		MY = za*Ji[0][1]*Hri + za*Ji[1][1]*Hsi + Ji[2][1]*0.5*m_h0[i]*Hi;
		MZ = za*Ji[0][2]*Hri + za*Ji[1][2]*Hsi + Ji[2][2]*0.5*m_h0[i]*Hi;

		// calculate deformation gradient F
		F[0][0] += NX*x + MX*dx; F[0][1] += NY*x + MY*dx; F[0][2] += NZ*x + MZ*dx;
		F[1][0] += NX*y + MX*dy; F[1][1] += NY*y + MY*dy; F[1][2] += NZ*y + MZ*dy;
		F[2][0] += NX*z + MX*dz; F[2][1] += NY*z + MY*dz; F[2][2] += NZ*z + MZ*dz;
	}

	double V = F.det();
	if (V <= 0) throw NegativeJacobian(m_nID, n, V, this);

	return V;
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
