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
#include "FEPipetteAspiration.h"
#include "FEBioMech.h"
#include <FECore/FEFacetSet.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEContactInterface.h"
#include "FESlidingElasticInterface.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEPipetteAspiration, FESurfaceLoad)
	ADD_PARAMETER(m_pressure, "pressure")->setUnits(UNIT_PRESSURE)->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
    ADD_PARAMETER(m_radius, "radius")->setUnits(UNIT_LENGTH)->setLongName("Pipette radius");
    ADD_PARAMETER(m_center, "center")->setUnits(UNIT_LENGTH)->setLongName("Pipette center");
    ADD_PARAMETER(m_normal, "normal")->setLongName("Pipette normal");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEPipetteAspiration::FEPipetteAspiration(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_pressure = 0.0;
    m_radius = 0.0;
    m_center = vec3d(0,0,0);
    m_normal = vec3d(0,0,1);
    m_contact = -1;
}

//-----------------------------------------------------------------------------
bool FEPipetteAspiration::Init()
{
	FESurface& surf = GetSurface();
    // check if this surface is a contact surface (primary)
    FEModel& fem = *GetFEModel();
    if (fem.SurfacePairConstraints() > 0)
    {
        // loop over all contact interfaces
        for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
        {
            FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
            FESlidingElasticInterface* pbw = dynamic_cast<FESlidingElasticInterface*>(pci);
            if (pbw) {
                FESurface& ps = *pci->GetPrimarySurface();
                // for now, use quick and dirty way to figure out if these surfaces are the same
                if (ps.Elements() == surf.Elements()) {
                    m_contact = i;
                    surf = ps;
                    break;
                }
            }
        }
    }
    else if (m_contact == -1) {
        feLogWarning("Pipette surface does not match any contact surface!");
        return false;
    }

	// get the degrees of freedom
	m_dof.Clear();
	if (surf.IsShellBottom() == false)
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
void FEPipetteAspiration::PrepStep()
{
    FEModel& fem = *GetFEModel();
    FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(m_contact));
    double psf = pci->GetPenaltyScaleFactor();
    FESlidingElasticInterface* sei = dynamic_cast<FESlidingElasticInterface*>(pci);
    FESurface& surf = *pci->GetPrimarySurface();
    m_tag.assign(surf.Elements(), std::vector<bool>(FEElement::MAX_INTPOINTS,false));

    // loop over all faces of the surface
    for (int i=0; i<surf.Elements(); ++i) {
        FESurfaceElement& el = surf.Element(i);
        // loop over all integration points of the surface
        for (int j=0; j<el.GaussPoints(); ++j) {
            FEMaterialPoint& pt = *el.GetMaterialPoint(j);
            // get the sliding elastic interface data at this material point
            FESlidingElasticSurface::Data& data = static_cast<FESlidingElasticSurface::Data&>(pt);
            // penalty
            double eps = sei->m_epsn*data.m_epsn*psf;
            // normal gap
            double g = data.m_gap;
            // normal traction Lagrange multiplier
            double Lm = data.m_Lmd;
            double pn = sei->m_btension ? (Lm + eps*g) : MBRACKET(Lm + eps*g);
            vec3d dx = pt.m_rt - m_center;
            double r = (dx - m_normal*(dx*m_normal)).Length();
            if ((r <= m_radius) && (pn == 0))
                m_tag[i][j] = true;
        }
    }
}

//-----------------------------------------------------------------------------
void FEPipetteAspiration::LoadVector(FEGlobalVector& R)
{
    FEModel& fem = *GetFEModel();
    FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(m_contact));
    FESurface& surf = *pci->GetPrimarySurface();

	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
		
        FESurfaceElement& sel = *pt.SurfaceElement();
		// evaluate pressure at this material point
		double P = (m_tag[sel.m_lid][pt.m_index]) ? -m_pressure(pt) : 0;
		if (surf.IsShellBottom()) P = -P;
        
		double J = (pt.dxr ^ pt.dxs).norm();

		// force vector
		vec3d N = (pt.dxr ^ pt.dxs); N.unit();
		vec3d t = N*P;

		double H_u = dof_a.shape;

		val[0] = H_u*t.x*J;
		val[1] = H_u*t.y*J;
		val[2] = H_u*t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FEPipetteAspiration::StiffnessMatrix(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(m_contact));
    FESurface& surf = *pci->GetPrimarySurface();

	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {

        FESurfaceElement& sel = *mp.SurfaceElement();
        // evaluate pressure at this material point
        double P = (m_tag[sel.m_lid][mp.m_index]) ? -m_pressure(mp) : 0;
       if (surf.IsShellBottom()) P = -P;
        
		double H_i  = dof_a.shape;
		double Gr_i = dof_a.shape_deriv_r;
		double Gs_i = dof_a.shape_deriv_s;

		double H_j  = dof_b.shape;
		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		vec3d vab(0,0,0);
        vab = (mp.dxs*Gr_j - mp.dxr*Gs_j)*(P*H_i);

		mat3da K(vab);
		kab.set(0, 0, K);
	});
}
