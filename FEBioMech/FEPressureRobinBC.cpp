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
#include "FEPressureRobinBC.h"
#include "FEBioMech.h"
#include "FEMechModel.h"
#include <FECore/FEFacetSet.h>
#include <FECore/FEModel.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEPressureRobinBC, FESurfaceLoad)
	ADD_PARAMETER(m_epsk   , "spring_eps")->setUnits("P/L");
	ADD_PARAMETER(m_bshellb , "shell_bottom");
    ADD_PARAMETER(m_nRB, "body")->setEnums("$(rigid_materials)");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEPressureRobinBC::FEPressureRobinBC(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_epsk = 0.0;
	m_bshellb = false;
    m_nRB = -1;
    m_rb = nullptr;
}

//-----------------------------------------------------------------------------
bool FEPressureRobinBC::Init()
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
void FEPressureRobinBC::Update()
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
void FEPressureRobinBC::LoadVector(FEGlobalVector& R)
{
    if (GetFEModel()->GetTime().currentTime == 0) return;
    
    quatd q;
    if (m_rb) q = m_rb->GetRotation();
    
	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
		
        FESurfaceElement& el = *pt.SurfaceElement();
        vec3d G[2];
        surf.CoBaseVectors0(el, pt.m_index, G);
        
		// evaluate referential normal and area at this point
        vec3d nr = G[0] ^ G[1];
        double J = nr.unit();
        
        double epsk = m_epsk(pt);
        vec3d u = pt.m_rt - pt.m_r0;
        double un = u*nr;

        // prescribe a pressure only when the boundary is moving in the direction of the normal
		double p = (un > 0) ? epsk*un : 0;
		if (m_bshellb) p = -p;

		// force vector
        q.RotateVector(nr);
		vec3d t = -nr*p;

		double Na = dof_a.shape;

		val[0] = Na*t.x*J;
		val[1] = Na*t.y*J;
		val[2] = Na*t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FEPressureRobinBC::StiffnessMatrix(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    if (fem.GetTime().currentTime == 0) return;
    
    quatd q;
    if (m_rb) q = m_rb->GetRotation().Conjugate();
    
    // evaluate the integral
	FESurface& surf = GetSurface();
    
    for (int k=0; k<surf.Elements(); ++k) {
        FESurfaceElement& el = surf.Element(k);
        int nint = el.GaussPoints();
        int neln = el.Nodes();
        FEElementMatrix ke(el);
        if (m_rb) ke.resize(3*neln, 3*(neln+1));
        else ke.resize(3*neln, 3*neln);
        ke.zero();
        vector<int> lmr(neln*3);
        vector<int> lmc((neln+1)*3);
        for (int n=0; n<nint; ++n) {
            vec3d G[2];
            surf.CoBaseVectors0(el, n, G);
            
            // evaluate referential normal and area at this point
            vec3d nr = G[0] ^ G[1];
            double J = nr.unit();
            mat3ds Nr = dyad(nr);
            
            FEMaterialPoint* mp = el.GetMaterialPoint(n);
            FESurfaceMaterialPoint& pt = *(mp->ExtractData<FESurfaceMaterialPoint>());
            double epsk = m_epsk(pt);
            vec3d u = mp->m_rt - mp->m_r0;
            double un = u*nr;
            q.RotateVector(nr);
            double* H = el.H(n);
            
            for (int i=0, i3=0; i<neln; ++i, i3+=3) {
                FENode& node = mesh.Node(el.m_node[i]);
                vector<int>& id = node.m_ID;
                lmr[i3  ] = id[m_dof[0]];
                lmr[i3+1] = id[m_dof[1]];
                lmr[i3+2] = id[m_dof[2]];
                double Ni = H[i];
                for (int j=0, j3=0; j<neln; ++j, j3+=3) {
                    double Nj = H[j];
                    mat3d Kij = q.RotationMatrix()*Nr*(epsk*Ni*Nj*J);
                    
                    ke[i3  ][j3  ] += Kij(0,0);
                    ke[i3  ][j3+1] += Kij(0,1);
                    ke[i3  ][j3+2] += Kij(0,2);
                    
                    ke[i3+1][j3  ] += Kij(1,0);
                    ke[i3+1][j3+1] += Kij(1,1);
                    ke[i3+1][j3+2] += Kij(1,2);
                    
                    ke[i3+2][j3  ] += Kij(2,0);
                    ke[i3+2][j3+1] += Kij(2,1);
                    ke[i3+2][j3+2] += Kij(2,2);
                }
                if (m_rb) {
                    lmc = lmr;
                    lmc.push_back(m_rb->m_LM[3]);
                    lmc.push_back(m_rb->m_LM[4]);
                    lmc.push_back(m_rb->m_LM[5]);
                    vec3d tmp = nr*(Ni*epsk*un*J);
                    mat3da Ki(tmp);
                    int j3 = 3*neln;
                    ke[i3  ][j3  ] += Ki(0,0);
                    ke[i3  ][j3+1] += Ki(0,1);
                    ke[i3  ][j3+2] += Ki(0,2);
                    
                    ke[i3+1][j3  ] += Ki(1,0);
                    ke[i3+1][j3+1] += Ki(1,1);
                    ke[i3+1][j3+2] += Ki(1,2);
                    
                    ke[i3+2][j3  ] += Ki(2,0);
                    ke[i3+2][j3+1] += Ki(2,1);
                    ke[i3+2][j3+2] += Ki(2,2);
                }
            }
        }
        if (m_rb)
            ke.SetIndices(lmr,lmc);
        else
            ke.SetIndices(lmr);
        LS.Assemble(ke);
    }
}
