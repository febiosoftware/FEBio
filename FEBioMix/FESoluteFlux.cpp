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
#include "FESoluteFlux.h"
#include <FECore/FEAnalysis.h>
#include <FECore/FEFacetSet.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESoluteFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux   , "flux");
	ADD_PARAMETER(m_blinear, "linear");
    ADD_PARAMETER(m_bshellb, "shell_bottom");
	ADD_PARAMETER(m_isol   , "solute_id")->setEnums("$(solutes)");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteFlux::FESoluteFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofC(pfem), m_dofU(pfem)
{ 
	m_flux = 1.0;
	m_blinear = false; 
    m_bshellb = false;
	m_isol = -1;
}
	
//-----------------------------------------------------------------------------
//! allocate storage
void FESoluteFlux::SetSurface(FESurface* ps)
{ 
	FESurfaceLoad::SetSurface(ps);
	m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! serialization
void FESoluteFlux::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);

	if (ar.IsShallow() == false)
	{
		ar & m_dofC & m_dofU;
	}
}

//-----------------------------------------------------------------------------
bool FESoluteFlux::Init()
{
	if (m_isol == -1) return false;

	// set up the dof lists
	m_dofC.Clear();
	m_dofU.Clear();
	if (m_bshellb == false)
	{ 
		m_dofC.AddDof(GetDOFIndex("concentration", m_isol - 1));

		m_dofU.AddDof(GetDOFIndex("x"));
		m_dofU.AddDof(GetDOFIndex("y"));
		m_dofU.AddDof(GetDOFIndex("z"));
	}
	else
	{
		m_dofC.AddDof(GetDOFIndex("shell concentration", m_isol - 1));

		m_dofU.AddDof(GetDOFIndex("sx"));
		m_dofU.AddDof(GetDOFIndex("sy"));
		m_dofU.AddDof(GetDOFIndex("sz"));

	}
    m_dof.AddDofs(m_dofU);
    m_dof.AddDofs(m_dofC);
	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FESoluteFlux::LoadVector(FEGlobalVector& R)
{
	double dt = CurrentTimeIncrement();

	m_psurf->SetShellBottom(m_bshellb);

	FESoluteFlux* flux = this;
	m_psurf->LoadVector(R, m_dofC, m_blinear, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {

		double wr = flux->m_flux(mp);
		if (flux->m_bshellb) wr = -wr;

		vec3d dxt = mp.dxr ^ mp.dxs;

		// volumetric flow rate
		double f = dxt.norm()*wr* dt;

		double H_i = dof_a.shape;
		fa[0] = H_i * f;
	});
}

//-----------------------------------------------------------------------------
void FESoluteFlux::StiffnessMatrix(FELinearSystem& LS)
{
	// time increment
	double dt = CurrentTimeIncrement();

	m_psurf->SetShellBottom(m_bshellb);
	
	// evaluate the stiffness contribution
	FESoluteFlux* flux = this;
	m_psurf->LoadStiffness(LS, m_dofC, m_dofU, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

		// shape functions and derivatives
		double H_i  = dof_a.shape;
		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		double wr = flux->m_flux(mp);
		if (flux->m_bshellb) wr = -wr;

		// calculate surface normal
		vec3d dxt = mp.dxr ^ mp.dxs;

		// calculate stiffness component
		vec3d t1 = dxt / dxt.norm()*wr;
		vec3d t2 = mp.dxs*Gr_j - mp.dxr*Gs_j;
		vec3d kab = (t1 ^ t2)*(H_i)*dt;

		Kab[0][0] = kab.x;
		Kab[0][1] = kab.y;
		Kab[0][2] = kab.z;
	});
}
