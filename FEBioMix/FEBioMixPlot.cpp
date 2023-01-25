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
#include "FEBioMixPlot.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicShellDomain.h"
#include "FEBiphasicSoluteDomain.h"
#include "FEBiphasicSoluteSolidDomain.h"
#include "FEBiphasicSoluteShellDomain.h"
#include "FETriphasicDomain.h"
#include "FEMultiphasicSolidDomain.h"
#include "FEMultiphasicShellDomain.h"
#include <FEBioMech/FEElasticSolidDomain.h>
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FEBiphasicContactSurface.h"
#include "FEBioMech/FEDonnanEquilibrium.h"
#include "FEBioMech/FEElasticMixture.h"
#include <FECore/FEModel.h>
#include <FECore/writeplot.h>
#include <FECore/FEEdgeList.h>

//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
// Plot local fluid load support
bool FEPlotLocalFluidLoadSupport::Save(FESurface& surf, FEDataStream& a)
{
    FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
    if (pcs == 0) return false;
    
    writeElementValue<double>(surf, a, [=](int nface) {
        double gn;
        pcs->GetLocalFLS(nface, gn);
        return gn;
    });
    return true;
}

//-----------------------------------------------------------------------------
// Plot effective friction coefficient
bool FEPlotEffectiveFrictionCoeff::Save(FESurface& surf, FEDataStream& a)
{
    FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
    if (pcs == 0) return false;
    
    writeElementValue<double>(surf, a, [=](int nface) {
        double gn;
        pcs->GetMuEffective(nface, gn);
        return gn;
    });
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMixtureFluidFlowRate::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        for (int j = 0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el, 0, 0)*pcs->FaceArea(el);
        }
        m_binit = false;
    }
    
    // calculate net flow rate normal to this surface
    for (int j = 0; j<NF; ++j)
    {
        FESurfaceElement& el = pcs->Element(j);
        
        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // evaluate the average fluid flux in this element
            int nint = pe->GaussPoints();
            vec3d w(0, 0, 0);
            for (int n = 0; n<nint; ++n)
            {
                FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                FEBiphasicMaterialPoint* ptf = mp.ExtractData<FEBiphasicMaterialPoint>();
                if (ptf) w += ptf->m_w;
            }
            w /= nint;
            
            // Evaluate contribution to net flow rate across surface.
            fn += w*m_area[j];
        }
    }
    
    // save results
    a << fn;
    
    return true;
}

//-----------------------------------------------------------------------------
// Plot contact gap
bool FEPlotPressureGap::Save(FESurface& surf, FEDataStream& a)
{
	FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	writeNodalProjectedElementValues<double>(surf, a, [](const FEMaterialPoint& mp) {
		const FEBiphasicContactPoint* pt = mp.ExtractData<FEBiphasicContactPoint>();
		return (pt ? pt->m_pg : 0);
	});

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidForce::Save(FESurface &surf, FEDataStream &a)
{
	FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	vec3d fn = pcs->GetFluidForce();
	a << fn;
    
	return true;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidForce2
{
public:
	FEFluidForce2(FESurface& surf, vector<double>& nodalPressures) : m_surf(surf), m_pe(nullptr), m_nodalPressures(nodalPressures) {}

	FEFluidForce2(const FEFluidForce2& fl) : m_surf(fl.m_surf), m_nodalPressures(fl.m_nodalPressures), m_pe(0) {}

	vec3d operator ()(const FEMaterialPoint& mp)
	{
		if (m_pe != mp.m_elem)
		{ 
			m_pe = mp.m_elem;
			int neln = m_pe->Nodes();
			for (int j = 0; j<neln; ++j) pn[j] = m_nodalPressures[m_pe->m_node[j]];
		}

		FESurfaceElement& face = static_cast<FESurfaceElement&>(*m_pe);

		// get the base vectors
		vec3d g[2];
		m_surf.CoBaseVectors(face, mp.m_index, g);

		// normal (magnitude = area)
		vec3d n = g[0] ^ g[1];

		// gauss weight
		double w = face.GaussWeights()[mp.m_index];

		// fluid pressure
		double p = face.eval(pn, mp.m_index);

		// contact force
		return n*(w*p);
	}

private:
	FESurface&	m_surf;
	FEElement* m_pe;
	double pn[FEElement::MAX_NODES];
	vector<double>& m_nodalPressures;
};

bool FEPlotFluidForce2::Save(FESurface &surf, FEDataStream &a)
{
	// get number of facets
	int NF = surf.Elements();
	if (NF == 0) return false;

	// this assumes that the surface sits on top of a single domain
	// so that we can figure out the domain from a single element
	FESurfaceElement& ref = surf.Element(0);
	if (ref.m_elem[0] == nullptr) return false;

	// get the element
	FEMesh& mesh = *surf.GetMesh();
	FEElement* el = ref.m_elem[0];
	if (el == 0) return false;

	// get the domain this element belongs to
	FEMeshPartition* dom = el->GetMeshPartition();
	if (dom == 0) return false;

	// see if this is a biphasic domain
	FEBiphasicSolidDomain* biphasicDomain = dynamic_cast<FEBiphasicSolidDomain*>(dom);
	if (biphasicDomain == 0) return false;

	// The biphasic solid domain contains actual nodal pressures.
	// we want to evaluate this over the surface. 
	vector<double> nodalPressures;
	biphasicDomain->GetNodalPressures(nodalPressures);

	// calculate element values
	writeIntegratedElementValue<vec3d>(surf, a, FEFluidForce2(surf, nodalPressures));

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidLoadSupport::Save(FESurface &surf, FEDataStream &a)
{
    FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
    if (pcs == 0) return false;
    
    double fn = pcs->GetFluidLoadSupport();
    a << fn;
    
    return true;
}

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FEMPSpecificStrainEnergy
{
public:
    FEMPSpecificStrainEnergy(FEMultiphasic* pm) : m_mat(pm) {}
    double operator()(const FEMaterialPoint& mp)
    {
        return m_mat->GetElasticMaterial()->StrainEnergyDensity(const_cast<FEMaterialPoint&>(mp))/m_mat->SolidReferentialApparentDensity(const_cast<FEMaterialPoint&>(mp));
    }
private:
    FEMultiphasic*    m_mat;
};

bool FEPlotMPSpecificStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEMultiphasic* pme = dom.GetMaterial()->ExtractProperty<FEMultiphasic>();
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FEMPSpecificStrainEnergy psi(pme);
        writeAverageElementValue<double>(dom, a, psi);
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotActualFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
    FEShellDomain& bsd = static_cast<FEShellDomain&>(dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&bd)) ||
		(dynamic_cast<FEBiphasicSoluteSolidDomain*>(&bd)) ||
		(dynamic_cast<FETriphasicDomain*     >(&bd)) ||
		(dynamic_cast<FEMultiphasicSolidDomain*   >(&bd)))
	{
		writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
			const FEBiphasicMaterialPoint* pt = mp.ExtractData<FEBiphasicMaterialPoint>();
			return (pt ? pt->m_pa : 0.0);
		});
		return true;
	}
    else if (dynamic_cast<FEElasticSolidDomain* >(&bd))
    {
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            FEElasticMixture* pem  = dynamic_cast<FEElasticMixture*> (dom.GetMaterial());
            
            if (pem == nullptr) return false;
            
            // extract fixed-charge density
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
                for (int k=0; k<pem->Materials(); ++k) {
                    FEDonnanEquilibriumMaterialPoint* pd = pt.GetPointData(k)->ExtractData<FEDonnanEquilibriumMaterialPoint>();
                    if (pd) ew += pd->m_p;
                }
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
 
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSolidStress::Save(FEDomain &dom, FEDataStream& a)
{
    FESolidDomain* bd = dynamic_cast<FESolidDomain*>(&dom);
    FEShellDomain* bsd = dynamic_cast<FEShellDomain*>(&dom);
    if (bd && (
               (dynamic_cast<FEBiphasicSolidDomain* >(bd)) ||
               (dynamic_cast<FEBiphasicSoluteSolidDomain*>(bd)) ||
               (dynamic_cast<FETriphasicDomain*     >(bd)) ||
               (dynamic_cast<FEMultiphasicSolidDomain*>(bd))))
    {
        writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
            const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
            return (pt ? pt->m_ss : mat3ds(0.0));
        });
        return true;
    }
    else if (bsd && (
                     (dynamic_cast<FEBiphasicShellDomain*>(bsd)) ||
                     (dynamic_cast<FEBiphasicSoluteShellDomain*>(bsd)) ||
                     (dynamic_cast<FEMultiphasicShellDomain*>(bsd))
                     ))
    {
        writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
            const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
            return (pt ? pt->m_ss : mat3ds(0.0));
        });
        return true;
    }
    
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidFlux::Save(FEDomain &dom, FEDataStream& a)
{
	FESolidDomain* bd = dynamic_cast<FESolidDomain*>(&dom);
    FEShellDomain* bsd = dynamic_cast<FEShellDomain*>(&dom);
	if (bd && (
        (dynamic_cast<FEBiphasicSolidDomain* >(bd)) ||
		(dynamic_cast<FEBiphasicSoluteSolidDomain*>(bd)) ||
		(dynamic_cast<FETriphasicDomain*     >(bd)) ||
		(dynamic_cast<FEMultiphasicSolidDomain*>(bd))))
	{
		writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
			const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
			return (pt ? pt->m_w : vec3d(0.0));
		});
		return true;
	}
    else if (bsd && (
                     (dynamic_cast<FEBiphasicShellDomain*>(bsd)) ||
                     (dynamic_cast<FEBiphasicSoluteShellDomain*>(bsd)) ||
                     (dynamic_cast<FEMultiphasicShellDomain*>(bsd))
                     ))
    {
		writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
			const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
			return (pt ? pt->m_w : vec3d(0.0));
		});
		return true;
    }

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotNodalFluidFlux::Save(FEDomain &dom, FEDataStream& a)
{
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&bd)) ||
		(dynamic_cast<FEBiphasicSoluteSolidDomain*>(&bd)) ||
		(dynamic_cast<FETriphasicDomain*     >(&bd)) ||
		(dynamic_cast<FEMultiphasicSolidDomain*>(&bd)))
	{
		writeNodalProjectedElementValues<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
			const FEBiphasicMaterialPoint* pt = mp.ExtractData<FEBiphasicMaterialPoint>();
			return pt->m_w;
		});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Finds the solute id of a solute with given ID nsol.
// This currently returns either nsol if a solute was found or -1 if not
int GetSoluteID(FEModel& fem, int nsol)
{
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd)
		{
			if (psd->GetID()-1 == nsol) return psd->GetID();
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
// Finds the id of a sbm with given ID nsbm.
// This currently returns either nsbm if a solute was found or -1 if not
int GetSBMID(FEModel& fem, int nsol)
{
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESBMData* psd = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
		if (psd)
		{
			if (psd->GetID()-1 == nsol) return psd->GetID();
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
// Finds the solute ID given the name of the solute
int GetSoluteID(FEModel& fem, const char* sz)
{
	string soluteName(sz);
	// find the solute with that name
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd)
		{
			if (psd->GetName() == soluteName) return psd->GetID();
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
// Finds the sbm ID given the name of the sbm
int GetSBMID(FEModel& fem, const char* sz)
{
	string sbmName(sz);
	// find the sbm with that name
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESBMData* psd = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
		if (psd)
		{
			if (psd->GetName() == sbmName) return psd->GetID();
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
// find the local SBM ID, given a global ID. If the material is not a 
// multiphasic material, this returns -1.
int GetLocalSBMID(FEMultiphasic* pmm, int nsbm)
{
	// figure out the SBM ID to export. This depends on the material type.
	int nsid = -1;

	// Check if this solute is present in this specific multiphasic mixture
	for (int i=0; i<pmm->SBMs(); ++i)
		if (pmm->GetSBM(i)->GetSBMID() == nsbm) {nsid = i; break;}
	
	return nsid;
}

//=================================================================================================
//-----------------------------------------------------------------------------
FEPlotActualSoluteConcentration::FEPlotActualSoluteConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i<ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
			m_sol.push_back(ps->GetID());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);
    SetUnits(UNIT_CONCENTRATION);
}

//-----------------------------------------------------------------------------
bool FEPlotActualSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int negs = 0;
	for (int i = 0; i<(int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
		if (lid[i] < 0) negs++;
	}
	if (negs == nsols) return false;

	// loop over all elements
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k=0; k<nsols; ++k)
		{
			int nsid = lid[k];
			if (nsid == -1) a << 0.f;
			else
			{
				// calculate average concentration
				double ew = 0;
				for (int j = 0; j<el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					ew += pm->GetActualSoluteConcentration(mp, nsid);
				}
				ew /= el.GaussPoints();
				a << ew;
			}
		}

	}
	return true;
}

//=================================================================================================
//-----------------------------------------------------------------------------
FEPlotPartitionCoefficient::FEPlotPartitionCoefficient(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
    DOFS& dofs = pfem->GetDOFS();
    int nsol = dofs.GetVariableSize("concentration");
    SetArraySize(nsol);
    
    // collect the names
    int ndata = pfem->GlobalDataItems();
    vector<string> s;
    for (int i = 0; i<ndata; ++i)
    {
        FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
        if (ps)
        {
            s.push_back(ps->GetName());
            m_sol.push_back(ps->GetID());
        }
    }
    assert(nsol == (int)s.size());
    SetArrayNames(s);
}

//-----------------------------------------------------------------------------
bool FEPlotPartitionCoefficient::Save(FEDomain &dom, FEDataStream& a)
{
    FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    if (pm == 0) return false;
    
    // figure out the local solute IDs. This depends on the material
    int nsols = (int)m_sol.size();
    vector<int> lid(nsols, -1);
    int negs = 0;
    for (int i = 0; i<(int)m_sol.size(); ++i)
    {
        lid[i] = pm->FindLocalSoluteID(m_sol[i]);
        if (lid[i] < 0) negs++;
    }
    if (negs == nsols) return false;
    
    // loop over all elements
    int N = dom.Elements();
    for (int i = 0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        for (int k=0; k<nsols; ++k)
        {
            int nsid = lid[k];
            if (nsid == -1) a << 0.f;
            else
            {
                // calculate average concentration
                double ew = 0;
                for (int j = 0; j<el.GaussPoints(); ++j)
                {
                    FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					ew += pm->GetPartitionCoefficient(mp, nsid);
                }
                ew /= el.GaussPoints();
                a << ew;
            }
        }
        
    }
    return true;
}

//-----------------------------------------------------------------------------
FEPlotSoluteFlux::FEPlotSoluteFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY_VEC3F, FMT_ITEM)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i<ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
			m_sol.push_back(ps->GetID());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);
    SetUnits(UNIT_MOLAR_FLUX);
}

//-----------------------------------------------------------------------------
bool FEPlotSoluteFlux::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if ((pm == 0) || (pm->Solutes() == 0)) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int nsc = 0;
	for (int i = 0; i<(int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
		if (lid[i] != -1) nsc++;
	}
	if (nsc == 0) return false;

	for (int i = 0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k=0; k<nsols; ++k)
		{
			int nsid = lid[k];
			if (nsid == -1) a << vec3d(0, 0, 0);
			else
			{
				// calculate average flux
				vec3d ew = vec3d(0, 0, 0);
				for (int j = 0; j<el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					ew += pm->GetSoluteFlux(mp, nsid);
				}

				ew /= el.GaussPoints();

				a << ew;
			}
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotSoluteVolumetricFlux::FEPlotSoluteVolumetricFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY_VEC3F, FMT_ITEM)
{
    DOFS& dofs = pfem->GetDOFS();
    int nsol = dofs.GetVariableSize("concentration");
    SetArraySize(nsol);
    
    // collect the names
    int ndata = pfem->GlobalDataItems();
    vector<string> s;
    for (int i = 0; i<ndata; ++i)
    {
        FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
        if (ps)
        {
            s.push_back(ps->GetName());
            m_sol.push_back(ps->GetID());
        }
    }
    assert(nsol == (int)s.size());
    SetArrayNames(s);
    SetUnits(UNIT_VELOCITY);
}

//-----------------------------------------------------------------------------
bool FEPlotSoluteVolumetricFlux::Save(FEDomain &dom, FEDataStream& a)
{
    FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    if ((pm == 0) || (pm->Solutes() == 0)) return false;
    
    // figure out the local solute IDs. This depends on the material
    int nsols = (int)m_sol.size();
    vector<int> lid(nsols, -1);
    int nsc = 0;
    for (int i = 0; i<(int)m_sol.size(); ++i)
    {
        lid[i] = pm->FindLocalSoluteID(m_sol[i]);
        if (lid[i] != -1) nsc++;
    }
    if (nsc == 0) return false;
    
    for (int i = 0; i<dom.Elements(); ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        for (int k=0; k<nsols; ++k)
        {
            int nsid = lid[k];
            if (nsid == -1) a << vec3d(0, 0, 0);
            else
            {
                // calculate average flux
                vec3d ew = vec3d(0, 0, 0);
                for (int j = 0; j<el.GaussPoints(); ++j)
                {
                    FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                    FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                    
                    if (pt && (pt->m_ca[nsid] > 0)) ew += pt->m_j[nsid]/pt->m_ca[nsid];
                }
                
                ew /= el.GaussPoints();
                
                a << ew;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotOsmolarity::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    FESolidDomain* sdom = dynamic_cast<FESolidDomain*>(&dom);
    FEShellDomain* ldom = dynamic_cast<FEShellDomain*>(&dom);
    if (sdom && psm)
	{
		writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
			double ew = psm->GetOsmolarity(mp);
			return ew;
		});

        return true;
    }
    else if (sdom && dynamic_cast<FEElasticSolidDomain*>(&dom)) {
        for (int i=0; i<sdom->Elements(); ++i)
        {
            FESolidElement& el = sdom->Element(i);
            FEElasticMixture* pem  = dynamic_cast<FEElasticMixture*> (dom.GetMaterial());
            
            if (pem == nullptr) return false;
            
            // extract fixed-charge density
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
                for (int k=0; k<pem->Materials(); ++k) {
                    FEDonnanEquilibriumMaterialPoint* pd = pt.GetPointData(k)->ExtractData<FEDonnanEquilibriumMaterialPoint>();
                    if (pd) ew += pd->m_osm;
                }
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
    else if (ldom && (
             dynamic_cast<FEBiphasicSoluteShellDomain*>(&dom) ||
             dynamic_cast<FEMultiphasicShellDomain*>(&dom))) {

		writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
			const FESolutesMaterialPoint* pt = mp.ExtractData<FESolutesMaterialPoint>();
			double ew = 0.0;
			for (int isol = 0; isol<(int)pt->m_ca.size(); ++isol) ew += pt->m_ca[isol];
			return ew;
		});

		return true;
    }
	return false;
}

//=================================================================================================
// FEPlotSBMConcentration
//=================================================================================================

//-----------------------------------------------------------------------------
FEPlotSBMConcentration::FEPlotSBMConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	// count SBMs
	int sbms = 0;
	int ndata = pfem->GlobalDataItems();
	vector<string> names;
	for (int i=0; i<ndata; ++i)
	{
		FESBMData* sbm = dynamic_cast<FESBMData*>(pfem->GetGlobalData(i));
		if (sbm) 
		{
			names.push_back(sbm->GetName());
			m_sbm.push_back(sbm->GetID());
			sbms++;
		}
	}

	SetArraySize(sbms);
	SetArrayNames(names);
    SetUnits(UNIT_CONCENTRATION);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local SBM IDs. This depend on the material
	int nsbm = (int)m_sbm.size();
	vector<int> lid(nsbm, -1);
	for (int i=0; i<(int)m_sbm.size(); ++i)
	{
		lid[i] = GetLocalSBMID(pm, m_sbm[i]);
	}

	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k=0; k<nsbm; ++k)
		{
			int nk = lid[k];
			if (nk == -1) a << 0.f;
			else
			{
				// calculate average concentration
				double ew = 0;
				for (int j = 0; j<el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());

					if (pt) ew += pm->SBMConcentration(mp, nk);
				}
				ew /= el.GaussPoints();

				a << ew;
			}
		}
	}
	return true;
}

//=================================================================================================
// FEPlotSBMArealConcentration
//=================================================================================================

//-----------------------------------------------------------------------------
FEPlotSBMArealConcentration::FEPlotSBMArealConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
    // count SBMs
    int sbms = 0;
    int ndata = pfem->GlobalDataItems();
    vector<string> names;
    for (int i=0; i<ndata; ++i)
    {
        FESBMData* sbm = dynamic_cast<FESBMData*>(pfem->GetGlobalData(i));
        if (sbm)
        {
            names.push_back(sbm->GetName());
            m_sbm.push_back(sbm->GetID());
            sbms++;
        }
    }
    
    SetArraySize(sbms);
    SetArrayNames(names);
    SetUnits(UNIT_MOLAR_AREAL_CONCENTRATION);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMArealConcentration::Save(FEDomain &dom, FEDataStream& a)
{
    FEShellDomain* bsd = static_cast<FEShellDomain*>(&dom);
    if (bsd == nullptr) return false;
    
    FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
    if (pm == 0) return false;
    
    // figure out the local SBM IDs. This depend on the material
    int nsbm = (int)m_sbm.size();
    vector<int> lid(nsbm, -1);
    for (int i=0; i<(int)m_sbm.size(); ++i)
    {
        lid[i] = GetLocalSBMID(pm, m_sbm[i]);
    }
    
    int N = dom.Elements();
    for (int i = 0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        for (int k=0; k<nsbm; ++k)
        {
            int nk = lid[k];
            if (nk == -1) a << 0.f;
            else
            {
                // calculate average concentration
                double ew = 0;
                for (int j = 0; j<el.GaussPoints(); ++j)
                {
                    FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                    FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                    
                    if (pt) ew += pm->SBMArealConcentration(mp, nk);
                }
                ew /= el.GaussPoints();
                
                a << ew;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotElectricPotential::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if (psm == nullptr) return false;
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		return psm->GetElectricPotential(mp);
		});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentDensity::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if (psm == nullptr) return false;
	writeAverageElementValue<vec3d>(dom, a, [=](const FEMaterialPoint& mp) {
		return psm->GetCurrentDensity(mp);
		});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialSolidVolumeFraction::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(dom.GetMaterial());
	if (pbm == nullptr) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		double phif0 = pbm->GetReferentialSolidVolumeFraction(mp);
		return phif0;
	});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPorosity::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasicSolidDomain* bmd = dynamic_cast<FEBiphasicSolidDomain*>(&dom);
    FEBiphasicShellDomain* bsd = dynamic_cast<FEBiphasicShellDomain*>(&dom);
    FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
    if (bmd || bsd || pmd || psd)
    {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
            const FEElasticMaterialPoint* et = (mp.ExtractData<FEElasticMaterialPoint>());
            const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
            return (pt ? (1 - pt->m_phi0t/et->m_J) : 0.0);
        });
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotPerm::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasic* bp = dom.GetMaterial()->ExtractProperty<FEBiphasic>();
    if (bp == 0) return false;

    writeAverageElementValue<mat3ds>(dom, a, [=](const FEMaterialPoint& mp) {
            return bp->Permeability(const_cast<FEMaterialPoint&>(mp));
        });
        return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFixedChargeDensity::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    FEElasticSolidDomain* ped = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (psm)
	{
		writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
            double cf = psm->GetFixedChargeDensity(mp);
			return cf;
		});
		return true;
	}
    else if (ped)
    {
        for (int i=0; i<ped->Elements(); ++i)
        {
            FESolidElement& el = ped->Element(i);
            FEElasticMixture* pem  = dynamic_cast<FEElasticMixture*> (dom.GetMaterial());
            
            if (pem == nullptr) return false;
            
            // extract fixed-charge density
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
                for (int k=0; k<pem->Materials(); ++k) {
                    FEDonnanEquilibriumMaterialPoint* pd = pt.GetPointData(k)->ExtractData<FEDonnanEquilibriumMaterialPoint>();
                    if (pd) ew += pd->m_cF;
                }
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialFixedChargeDensity::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    FEElasticSolidDomain* ped = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (psm)
	{
		writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
            double cf = psm->GetReferentialFixedChargeDensity(mp);
			return cf;
		});
		return true;
	}
    else if (ped)
    {
        for (int i=0; i<ped->Elements(); ++i)
        {
            FESolidElement& el = ped->Element(i);
            FEElasticMixture* pem  = dynamic_cast<FEElasticMixture*> (dom.GetMaterial());
            
            if (pem == nullptr) return false;
            
            // extract fixed-charge density
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
                for (int k=0; k<pem->Materials(); ++k) {
                    FEDonnanEquilibriumMaterialPoint* pd = pt.GetPointData(k)->ExtractData<FEDonnanEquilibriumMaterialPoint>();
                    if (pd) ew += pd->m_cFr;
                }
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicDomain*       pd = dynamic_cast<FEBiphasicDomain*      >(&dom);
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain*      ptd = dynamic_cast<FETriphasicDomain*     >(&dom);
	FEMultiphasicDomain*    pmd = dynamic_cast<FEMultiphasicDomain*   >(&dom);

	// special handling of mixed biphasic formulation
	if (pd)
	{
		// get the pressure dof index
		int dof_p = GetFEModel()->GetDOFIndex("p");
		if (dof_p == -1) return false;

		DOFS& dofs = GetFEModel()->GetDOFS();
		int varU = dofs.GetVariableIndex("displacement");
		int varP = dofs.GetVariableIndex("fluid pressure");

		int kd = dofs.GetVariableInterpolationOrder(varU);
		int kp = dofs.GetVariableInterpolationOrder(varP);
		if ((kd != 1) && (kp == 1))
		{
			int N = dom.Nodes();
			vector<double> p(N, 0.0);
			for (int i = 0; i < dom.Nodes(); ++i)
			{
				FENode& node = dom.Node(i);
				p[i] = node.get(dof_p);
			}

			FEEdgeList EL;
			if (EL.Create(&dom) == false) return false;

			for (int i = 0; i < EL.Edges(); ++i)
			{
				const FEEdgeList::EDGE& edge = EL.Edge(i);
				assert(edge.ntype == 3);
				p[edge.node[2]] = 0.5*(p[edge.node[0]] + p[edge.node[1]]);
			}

			a << p;

			return true;
		}
	}

	if (pd || psd || ptd || pmd)
	{
		// get the pressure dof index
		int dof_p = GetFEModel()->GetDOFIndex("p");
		if (dof_p == -1) return false;

		// write the nodal values
		writeNodalValues<double>(dom, a, [=, &dom](int i) {
			FENode& node = dom.Node(i);
			return node.get(dof_p);
		});

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveShellFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicShellDomain* pbsd = dynamic_cast<FEBiphasicShellDomain*>(&dom);
	FEBiphasicSoluteShellDomain* pbssd = dynamic_cast<FEBiphasicSoluteShellDomain*>(&dom);
	FEMultiphasicShellDomain* pmpsd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (pbsd || pbssd || pmpsd)
	{
		// get the pressure dof index
		int dof_q = GetFEModel()->GetDOFIndex("q");
		assert(dof_q != -1);

		// write the nodal values
		writeNodalValues<double>(dom, a, [=, &dom](int i) {
			FENode& node = dom.Node(i);
			return node.get(dof_q);
		});

		return true;
	}
	return false;
}

//=================================================================================================
// FEPlotEffectiveSoluteConcentration
//=================================================================================================

FEPlotEffectiveSoluteConcentration::FEPlotEffectiveSoluteConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_NODE)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i=0; i<ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
			m_sol.push_back(ps->GetID());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);
    SetUnits(UNIT_CONCENTRATION);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	// get the dof
	DOFS& dofs = GetFEModel()->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	if (nsol == -1) return false;

	// get the start index
	const int dof_C = GetFEModel()->GetDOFIndex("concentration", 0);

	FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int negs = 0;
	for (int i = 0; i<(int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
		if (lid[i] < 0) negs++;
	}
	if (negs == nsol) return false;

	// save the concentrations
	int N = dom.Nodes();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = dom.Node(i);
		for (int j=0; j<nsol; ++j)
		{
			double c = (lid[j] >= 0 ? node.get(dof_C + j) : 0.0);
			a << c;
		}
	}
	return true;
}

//=================================================================================================

//-----------------------------------------------------------------------------
FEPlotEffectiveShellSoluteConcentration::FEPlotEffectiveShellSoluteConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_NODE)
{
	m_nsol = 0;
}

//-----------------------------------------------------------------------------
// Resolve solute by name
bool FEPlotEffectiveShellSoluteConcentration::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*GetFEModel(), sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotEffectiveShellSoluteConcentration::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*GetFEModel(), nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveShellSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if (pm == 0) return false;

	// make sure we have a valid index
	int nsid = pm->FindLocalSoluteID(m_nsol);
	if (nsid == -1) return false;

	// get the dof
	const int dof_D = GetFEModel()->GetDOFIndex("shell concentration", nsid);
	if (dof_D == -1) return false;

	int N = dom.Nodes();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = dom.Node(i);
		a << node.get(dof_D);
	}
	return true;
}

//=================================================================================================
bool FEPlotReceptorLigandConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicSoluteSolidDomain* pbd = dynamic_cast<FEBiphasicSoluteSolidDomain*>(&dom);
    FEBiphasicSoluteShellDomain* psd = dynamic_cast<FEBiphasicSoluteShellDomain*>(&dom);
	if (pbd || psd)
	{
		writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
			const FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
			return (pt ? pt->m_sbmr[0] : 0.0);
		});
		return true;
	}
	return false;
}

//=================================================================================================
FEPlotSBMRefAppDensity::FEPlotSBMRefAppDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	// count SBMs
	int sbms = 0;
	int ndata = pfem->GlobalDataItems();
	vector<string> names;
	for (int i = 0; i<ndata; ++i)
	{
		FESBMData* sbm = dynamic_cast<FESBMData*>(pfem->GetGlobalData(i));
		if (sbm)
		{
			names.push_back(sbm->GetName());
			m_sbm.push_back(sbm->GetID());
			sbms++;
		}
	}

	SetArraySize(sbms);
	SetArrayNames(names);
    SetUnits(UNIT_DENSITY);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMRefAppDensity::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local SBM IDs. This depend on the material
	int nsbm = (int)m_sbm.size();
	vector<int> lid(nsbm, -1);
	for (int i = 0; i<(int)m_sbm.size(); ++i)
	{
		lid[i] = GetLocalSBMID(pm, m_sbm[i]);
	}

	// figure out the sbm ID to export. This depends on the material type.
	// loop over all elements
	for (int i = 0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k=0; k<nsbm; ++k)
		{
			int nsid = lid[k];
			if (nsid == -1) a << 0.f;
			else
			{
				// calculate average concentration
				double ew = 0;
				for (int j = 0; j<el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());

					if (st) ew += st->m_sbmr[nsid];
				}
				ew /= el.GaussPoints();

				a << ew;
			}
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveElasticity::Save(FEDomain &dom, FEDataStream& a)
{
    tens4ds c;
    
	if ((dom.Class() != FE_DOMAIN_SOLID)
        && (dom.Class() != FE_DOMAIN_SHELL)) return false;
	FESolidDomain* pbd = static_cast<FESolidDomain*>(&dom);
    FEShellDomain* psd = static_cast<FEShellDomain*>(&dom);
    
    FEBiphasic*       pb  = dynamic_cast<FEBiphasic      *> (dom.GetMaterial());
    FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
    FETriphasic*      ptp = dynamic_cast<FETriphasic     *> (dom.GetMaterial());
    FEMultiphasic*    pmp = dynamic_cast<FEMultiphasic   *> (dom.GetMaterial());
    if ((pb == 0) && (pbs == 0) && (ptp == 0) && (pmp == 0)) return false;

    if (pbd) {
        for (int i=0; i<pbd->Elements(); ++i)
        {
            FESolidElement& el = pbd->Element(i);
            
            int nint = el.GaussPoints();
            double f = 1.0 / (double) nint;
            
            // since the PLOT file requires floats we need to convert
            // the doubles to single precision
            // we output the average stress values of the gauss points
            tens4ds s(0.0);
            for (int j=0; j<nint; ++j)
            {
                FEMaterialPoint& pt = (*el.GetMaterialPoint(j)->ExtractData<FEMaterialPoint>());
                if (pb) c = pb->Tangent(pt);
                else if (pbs) c = pbs->Tangent(pt);
                else if (ptp) c = ptp->Tangent(pt);
                else if (pmp) c = pmp->Tangent(pt);
                
                s += c;
            }
            s *= f;
            
            // store average elasticity
            a << s;
        }
    }
    else if (psd) {
        for (int i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            int nint = el.GaussPoints();
            double f = 1.0 / (double) nint;
            
            // since the PLOT file requires floats we need to convert
            // the doubles to single precision
            // we output the average stress values of the gauss points
            tens4ds s(0.0);
            for (int j=0; j<nint; ++j)
            {
                FEMaterialPoint& pt = (*el.GetMaterialPoint(j)->ExtractData<FEMaterialPoint>());
                if (pb) c = pb->Tangent(pt);
                else if (pbs) c = pbs->Tangent(pt);
                else if (ptp) c = ptp->Tangent(pt);
                else if (pmp) c = pmp->Tangent(pt);
                
                s += c;
            }
            s *= f;
            
            // store average elasticity
            a << s;
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotOsmoticCoefficient::Save(FEDomain &dom, FEDataStream& a)
{
    if ((dom.Class() != FE_DOMAIN_SOLID) && (dom.Class() != FE_DOMAIN_SHELL)) return false;

	FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
	if (pm == nullptr) return false;

	FEOsmoticCoefficient* osm = pm->GetOsmoticCoefficient();
	if (osm == nullptr) return false;
    
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		double c = osm->OsmoticCoefficient(const_cast<FEMaterialPoint&>(mp));
		return c;
		});

    return true;
}
