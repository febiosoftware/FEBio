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
#include "FEMechModel.h"
#include <FECore/FEFacetSet.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FETractionRobinBC, FESurfaceLoad)
    ADD_PARAMETER(m_epsk   , "spring_eps")->setUnits("P/L");
	ADD_PARAMETER(m_bshellb , "shell_bottom");
    ADD_PARAMETER(m_nRB, "body")->setEnums("$(rigid_materials)");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FETractionRobinBC::FETractionRobinBC(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_epsk = 0.0;
    m_bshellb = false;
    m_nRB = -1;
    m_rb = nullptr;
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

    FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    if (m_nRB != -1) m_rb = fem.GetRigidBody(m_nRB);
    
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
    
    quatd q;
    if (m_rb) q = m_rb->GetRotation();
    
	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [=](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		// evaluate traction at this material point
        FESurfaceElement& el = *pt.SurfaceElement();
        vec3d G[2];
        surf.CoBaseVectors0(el, pt.m_index, G);
        
        // evaluate referential normal and area at this point
        vec3d nr = G[0] ^ G[1];
        double J = nr.unit();
        mat3ds Nr = dyad(nr);
        
        double epsk = m_epsk(pt);
        vec3d u = pt.m_rt - pt.m_r0;
        q.RotateVector(u);

        vec3d t = -u*epsk;
		if (m_bshellb) t = -t;

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
    quatd q;
    if (m_rb) q = m_rb->GetRotation();
    
    // evaluate the integral
    FESurface& surf = GetSurface();
    surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {
        
        FESurfaceElement& el = *pt.SurfaceElement();
        vec3d G[2];
        surf.CoBaseVectors0(el, pt.m_index, G);
        
        // evaluate referential normal and area at this point
        vec3d nr = G[0] ^ G[1];
        
        double J = nr.unit();
        
        double epsk = m_epsk(pt);

        double Na  = dof_a.shape;
        double Nb  = dof_b.shape;
        
        mat3d Kab = q.RotationMatrix()*(epsk*Na*Nb*J);
        kab.set(0, 0, Kab);
    });
}
