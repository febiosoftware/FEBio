/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEPoroTraction.h"
#include "FECore/FEModel.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPoroNormalTraction, FESurfaceLoad)
	ADD_PARAMETER(m_traction  , "traction" );
	ADD_PARAMETER(m_blinear   , "linear"   );
    ADD_PARAMETER(m_bshellb   , "shell_bottom");
	ADD_PARAMETER(m_beffective, "effective");
	ADD_PARAMETER(m_PC        , "value");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEPoroNormalTraction::FEPoroNormalTraction(FEModel* pfem) : FESurfaceLoad(pfem), m_dofUP(pfem)
{ 
	m_traction = 1.0;
	m_blinear = false; 
    m_bshellb = false;
	m_beffective = false;

	m_PC = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPoroNormalTraction::SetSurface(FESurface* ps)
{ 
	FESurfaceLoad::SetSurface(ps);
	m_PC.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
bool FEPoroNormalTraction::Init()
{
	FEModel* fem = GetFEModel();
	m_dofUP.Clear();
	if (m_bshellb == false)
	{
		m_dofUP.AddDof("x");
		m_dofUP.AddDof("y");
		m_dofUP.AddDof("z");
		if (m_dofUP.AddDof("p") == false) m_dofUP.AddDof(-1);
	}
	else
	{
		m_dofUP.AddDof(fem->GetDOFIndex("sx"));
		m_dofUP.AddDof(fem->GetDOFIndex("sy"));
		m_dofUP.AddDof(fem->GetDOFIndex("sz"));
		if (m_dofUP.AddDof("q") == false) m_dofUP.AddDof(-1);
	}
	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FEPoroNormalTraction::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	if (m_blinear) return;

	m_psurf->SetShellBottom(m_bshellb);

	bool bsymm = LS.IsSymmetric();

	FEPoroNormalTraction* traction = this;
	m_psurf->LoadStiffness(LS, m_dofUP, m_dofUP, [=](FESurfaceMaterialPoint& mp, int node_a, int node_b, matrix& Kab) {

			double* N = mp.m_shape;
			double* Gr = mp.m_shape_deriv_r;
			double* Gs = mp.m_shape_deriv_s;

			// traction at integration point
			double tr = traction->Traction(mp);
			if (traction->m_bshellb) tr = -tr;

			// calculate stiffness component
			Kab.zero();
			int i = node_a;
			int j = node_b;
			if (!bsymm) {
				// non-symmetric
				vec3d kab = (mp.dxs*Gr[j] - mp.dxr*Gs[j])*N[i] * tr;

				Kab[0][0] = 0;
				Kab[0][1] = -kab.z;
				Kab[0][2] = kab.y;

				Kab[1][0] += kab.z;
				Kab[1][1] += 0;
				Kab[1][2] += -kab.x;

				Kab[2][0] += -kab.y;
				Kab[2][1] += kab.x;
				Kab[2][2] += 0;

				// if prescribed traction is effective, add stiffness component
				if (traction->m_beffective)
				{
					vec3d kab = (mp.dxr ^ mp.dxs)* N[i] * N[j];

					Kab[0][3] = kab.x;
					Kab[1][3] = kab.y;
					Kab[2][3] = kab.z;
				}
			}
			else {
				// symmetric

				vec3d kab = ((mp.dxs*Gr[j] - mp.dxr*Gs[j])*N[i] - (mp.dxs*Gr[i] - mp.dxr*Gs[i])*N[j])*0.5 * tr;

				Kab[0][0] = 0;
				Kab[0][1] = -kab.z;
				Kab[0][2] = kab.y;

				Kab[1][0] = kab.z;
				Kab[1][1] = 0;
				Kab[1][2] = -kab.x;

				Kab[2][0] = -kab.y;
				Kab[2][1] = kab.x;
				Kab[2][2] = 0;

				// if prescribed traction is effective, add stiffness component
				if (traction->m_beffective)
				{
					vec3d kab = (mp.dxr ^ mp.dxs) * 0.5*N[i] * N[j];

					Kab[0][3] = kab.x;
					Kab[1][3] = kab.y;
					Kab[2][3] = kab.z;

					// TODO: This is not symmetric!
					Kab[0][3] = kab.x;
					Kab[1][3] = kab.y;
					Kab[2][3] = kab.z;
				}
		}
	});
}

//-----------------------------------------------------------------------------
double FEPoroNormalTraction::Traction(FESurfaceMaterialPoint& mp)
{
	FESurfaceElement& el = *mp.SurfaceElement();
	int neln = el.Nodes();

	// calculate nodal normal tractions
	double tr = m_traction*m_PC(mp);

	// if the prescribed traction is effective, evaluate the total traction
	if (m_beffective)
	{
		FEMesh& mesh = *m_psurf->GetMesh();

		// fluid pressure
		tr -= m_psurf->Evaluate(mp, m_dofUP[3]);
	}

	return tr;
}

//-----------------------------------------------------------------------------
void FEPoroNormalTraction::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	m_psurf->SetShellBottom(m_bshellb);

	FEPoroNormalTraction* traction = this;
	m_psurf->LoadVector(R, m_dofUP, m_blinear, [=](FESurfaceMaterialPoint& mp, int node_a, vector<double>& fa) {

		// traction at integration points
		double tr = traction->Traction(mp);
		if (traction->m_bshellb) tr = -tr;

		// force vector
		vec3d f = (mp.dxr ^ mp.dxs)*tr;

		double* H = mp.m_shape;
		fa[0] = H[node_a] * f.x;
		fa[1] = H[node_a] * f.y;
		fa[2] = H[node_a] * f.z;
		fa[3] = 0.0;
	});
}
