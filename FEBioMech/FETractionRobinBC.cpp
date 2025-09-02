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
#include "FETractionRobinBC.h"
#include "FEBioMech.h"
#include <FECore/FEFacetSet.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FETractionRobinBC, FESurfaceLoad)
    ADD_PARAMETER(m_epsk   , "spring_eps")->setUnits("P/L");
    ADD_PARAMETER(m_epsc   , "dashpot_eps")->setUnits("P.t/L");
	ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FETractionRobinBC::FETractionRobinBC(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_epsk = 0.0;
    m_epsc = 0.0;
    m_bshellb = false;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETractionRobinBC::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
// initialization
bool FETractionRobinBC::Init()
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// get the degrees of freedom
	m_dof.Clear();
	if (m_bshellb == false)
	{
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
	else
	{
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
	}
	if (m_dof.IsEmpty()) return false;

	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FETractionRobinBC::Update()
{
    FETimeInfo tp = GetFEModel()->GetTime();
    FESurface& surf = GetSurface();
    for (int i=0; i<surf.Elements(); ++i) {
        FESurfaceElement& el = surf.Element(i);
        int nint = el.GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint* mp = el.GetMaterialPoint(n);
            mp->Update(tp);
        }
    }
    surf.Update(tp);
}

//-----------------------------------------------------------------------------
void FETractionRobinBC::LoadVector(FEGlobalVector& R)
{
    if (GetFEModel()->GetTime().currentTime == 0) return;
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [=](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		// evaluate traction at this material point
        double epsk = m_epsk(pt);
        double epsc = m_epsc(pt);
        vec3d u = pt.m_rt - pt.m_r0;
        vec3d udot = (dt > 0) ? (pt.m_rt - pt.m_rp)/dt : vec3d(0,0,0);

        vec3d t = -u*epsk - udot*epsc;
		if (m_bshellb) t = -t;

		double J = (pt.dxr ^ pt.dxs).norm();

		double Na = dof_a.shape;

		val[0] = Na * t.x*J;
		val[1] = Na * t.y*J;
		val[2] = Na * t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FETractionRobinBC::StiffnessMatrix(FELinearSystem& LS)
{
    if (GetFEModel()->GetTime().currentTime == 0) return;
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // evaluate the integral
    FESurface& surf = GetSurface();
    surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {
        
        // evaluate pressure at this material point
        vec3d n = (pt.dxr ^ pt.dxs);
        double J = n.unit();
        
        mat3dd I(1);
        
        double epsk = m_epsk(pt);
        double epsc = m_epsc(pt);
        vec3d u = pt.m_rt - pt.m_r0;
        vec3d udot = (dt > 0) ? (pt.m_rt - pt.m_rp)/dt : vec3d(0,0,0);
        
        vec3d t = -u*epsk - udot*epsc;
        if (m_bshellb) t = -t;

        double Na  = dof_a.shape;
        double dNar = dof_a.shape_deriv_r;
        double dNas = dof_a.shape_deriv_s;
        
        double Nb  = dof_b.shape;
        double dNbr = dof_b.shape_deriv_r;
        double dNbs = dof_b.shape_deriv_s;
        
        mat3d Kab = (J*Na*Nb)*(epsk+epsc/dt)*mat3dd(1)
        +((t & n)*mat3da(pt.dxr*dNbs - pt.dxs*dNbr)*(Na*J));
        kab.set(0, 0, Kab);
    });
}
