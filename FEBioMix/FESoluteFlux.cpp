/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESoluteFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux   , "flux");
	ADD_PARAMETER(m_blinear, "linear");
    ADD_PARAMETER(m_bshellb, "shell_bottom");
	ADD_PARAMETER(m_isol   , "solute_id");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteFlux::FESoluteFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofC(pfem), m_dofU(pfem)
{ 
	m_flux = 1.0;
	m_blinear = false; 
    m_bshellb = false;
	m_isol = 0;
}
	
//-----------------------------------------------------------------------------
//! allocate storage
void FESoluteFlux::SetSurface(FESurface* ps)
{ 
	FESurfaceLoad::SetSurface(ps);
	m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
bool FESoluteFlux::Init()
{
	if (m_isol == -1) return false;

	// set up the dof lists
	FEModel* fem = GetFEModel();
	m_dofC.Clear();
	m_dofU.Clear();
	if (m_bshellb == false)
	{ 
		m_dofC.AddDof(fem->GetDOFIndex("concentration", m_isol - 1));

		m_dofU.AddDof(fem->GetDOFIndex("x"));
		m_dofU.AddDof(fem->GetDOFIndex("y"));
		m_dofU.AddDof(fem->GetDOFIndex("z"));
	}
	else
	{
		m_dofC.AddDof(fem->GetDOFIndex("shell concentration", m_isol - 1));

		m_dofU.AddDof(fem->GetDOFIndex("sx"));
		m_dofU.AddDof(fem->GetDOFIndex("sy"));
		m_dofU.AddDof(fem->GetDOFIndex("sz"));

	}
	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FESoluteFlux::UnpackLM(FEElement& el, vector<int>& lm)
{
	// TODO: remove this
}

//-----------------------------------------------------------------------------
void FESoluteFlux::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	double dt = tp.timeIncrement;

	FESoluteFlux* flux = this;
	m_psurf->LoadVector(R, m_dofC, m_blinear, [=](FESurfaceMaterialPoint& mp, int node_a, std::vector<double>& fa) {

		double wr = flux->m_flux(mp);
		if (flux->m_bshellb) wr = -wr;

		vec3d dxt = mp.dxr ^ mp.dxs;

		// volumetric flow rate
		double f = dxt.norm()*wr* dt;

		double* H = mp.m_shape;
		fa[0] = H[node_a] * f;
	});
}

//-----------------------------------------------------------------------------
void FESoluteFlux::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
	// time increment
	double dt = tp.timeIncrement;
	
	// evaluate the stiffness contribution
	FESoluteFlux* flux = this;
	m_psurf->LoadStiffness(psolver, m_dofC, m_dofU, [=](FESurfaceMaterialPoint& mp, int node_a, int node_b, matrix& Kab) {

		double* N = mp.m_shape;
		double* Gr = mp.m_shape_deriv_r;
		double* Gs = mp.m_shape_deriv_s;

		double wr = flux->m_flux(mp);
		if (flux->m_bshellb) wr = -wr;

		// calculate surface normal
		vec3d dxt = mp.dxr ^ mp.dxs;

		// calculate stiffness component
		int i = node_a;
		int j = node_b;

		vec3d t1 = dxt / dxt.norm()*wr;
		vec3d t2 = mp.dxs*Gr[j] - mp.dxr*Gs[j];
		vec3d kab = (t1 ^ t2)*(N[i])*dt;

		Kab[0][0] = kab.x;
		Kab[0][1] = kab.y;
		Kab[0][2] = kab.z;
	});
}
