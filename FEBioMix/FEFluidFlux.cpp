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
#include "FEFluidFlux.h"
#include "FECore/FESolver.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FEBioMix.h"
#include "FEBiphasicAnalysis.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux    , "flux"   );
	ADD_PARAMETER(m_blinear , "linear" );
    ADD_PARAMETER(m_bshellb , "shell_bottom");
	ADD_PARAMETER(m_bmixture, "mixture");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFluidFlux::FEFluidFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofP(pfem), m_dofU(pfem), m_dofV(pfem)
{ 
	m_blinear = false; 
	m_bmixture = false;
    m_bshellb = false;
	m_flux = 1.0;
}

//-----------------------------------------------------------------------------
void FEFluidFlux::SetSurface(FESurface* ps) 
{ 
	FESurfaceLoad::SetSurface(ps);
	m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
bool FEFluidFlux::Init()
{
	FEModel* fem = GetFEModel();

	// get the degrees of freedom
	if (m_bshellb == false)
	{
		m_dofU.AddDof(fem->GetDOFIndex("x"));
		m_dofU.AddDof(fem->GetDOFIndex("y"));
		m_dofU.AddDof(fem->GetDOFIndex("z"));

		m_dofV.AddDof(fem->GetDOFIndex("vx"));
		m_dofV.AddDof(fem->GetDOFIndex("vy"));
		m_dofV.AddDof(fem->GetDOFIndex("vz"));

		m_dofP.AddDof(fem->GetDOFIndex("p"));
	}
	else
	{
		m_dofU.AddDof(fem->GetDOFIndex("sx"));
		m_dofU.AddDof(fem->GetDOFIndex("sy"));
		m_dofU.AddDof(fem->GetDOFIndex("sz"));

		m_dofV.AddDof(fem->GetDOFIndex("svx"));
		m_dofV.AddDof(fem->GetDOFIndex("svy"));
		m_dofV.AddDof(fem->GetDOFIndex("svz"));

		m_dofP.AddDof(fem->GetDOFIndex("q"));
	}

	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FEFluidFlux::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		ar & m_dofU & m_dofV & m_dofP;
	}
}

//-----------------------------------------------------------------------------
vec3d FEFluidFlux::SolidVelocity(FESurfaceMaterialPoint& pt)
{
	FESurfaceElement& el = *pt.SurfaceElement();
	int n = pt.m_index;

	vec3d vt[FEElement::MAX_NODES];
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
	{
		FENode& nd = m_psurf->GetMesh()->Node(el.m_node[i]);
		vt[i] = nd.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
	}

	// shape functions at integration point n
	double* N = el.H(n);

	// solid velocity at integration point
	vec3d vr(0, 0, 0);
	for (int i = 0; i<neln; ++i)
	{
		vr += vt[i] * N[i];
	}

	return vr;
}

//-----------------------------------------------------------------------------
void FEFluidFlux::LoadVector(FEGlobalVector& R)
{
	FEModel* fem = GetFEModel();

	// only add the mixture term for transient analysis and when the m_bmixture flag is true
	bool bmixture = m_bmixture;
	if (fem->GetCurrentStep()->m_nanalysis == FEBiphasicAnalysis::STEADY_STATE) bmixture = false;

	// get time increment
	double dt = CurrentTimeIncrement();

	// integrate over surface
	FEFluidFlux* flux = this;
	m_psurf->LoadVector(R, m_dofP, m_blinear, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {

		// fluid flux
		double wr = flux->m_flux(mp);
		if (flux->m_bshellb) wr = -wr;

		// surface normal
		vec3d dxt = mp.dxr ^ mp.dxs;

		// volumetric flow rate
		double f = dxt.norm()*wr;

		// add optional mixture term
		if (bmixture)
		{
			// solid velocity at integration point
			vec3d vr = flux->SolidVelocity(mp);
			f -= vr*dxt;
		}

		// put it all together
		double H_p = dof_a.shape;
		fa[0] = H_p * f * dt;
	});
}

//-----------------------------------------------------------------------------
void FEFluidFlux::StiffnessMatrix(FELinearSystem& LS)
{
	FEModel& fem = *GetFEModel();
	double dt = CurrentTimeIncrement();

	FEFluidFlux* flux = this;
	bool btransient = (fem.GetCurrentStep()->m_nanalysis != FEBiphasicAnalysis::STEADY_STATE);

	if (!m_blinear || m_bmixture)
	{
		m_psurf->LoadStiffness(LS, m_dofP, m_dofU, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

			// fluid flux
			double wr = flux->m_flux(mp);
			if (m_bshellb) wr = -wr;

			// calculate surface normal
			vec3d dxt = mp.dxr ^ mp.dxs;

			// shape functions and derivatives at integration point
			double N_p = dof_a.shape;

			double N_u = dof_b.shape;
			double Gr_u = dof_b.shape_deriv_r;
			double Gs_u = dof_b.shape_deriv_s;

			// calculate stiffness component
			int i = dof_a.index;
			int j = dof_b.index;

			vec3d kab(0, 0, 0);
			if (flux->m_bmixture == false)
			{
				vec3d t1 = (dxt / dxt.norm())*wr;
				vec3d t2 = mp.dxs*Gr_u - mp.dxr*Gs_u;
				kab = (t1^t2)*dt*N_p;
			}
			else if (btransient)
			{
				kab = (dxt*N_u)*N_p;
			}

			Kab[0][0] += kab.x;
			Kab[0][1] += kab.y;
			Kab[0][2] += kab.z;
		});
	}
}
