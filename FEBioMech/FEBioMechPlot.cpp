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
#include "FEBioMechPlot.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FEReactiveFatigue.h"
#include "FEReactivePlasticity.h"
#include "FEReactivePlasticDamage.h"
#include "FERemodelingElasticMaterial.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEElasticShellDomainOld.h"
#include "FEElasticEASShellDomain.h"
#include "FEElasticANSShellDomain.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FEUT4Domain.h"
#include "FEContactSurface.h"
#include "FERigidBody.h"
#include <FECore/FESPRProjection.h>
#include "FEUncoupledElasticMixture.h"
#include "FERigidMaterial.h"
#include "FEVolumeConstraint.h"
#include "FEFacet2FacetSliding.h"
#include "FEMortarSlidingContact.h"
#include "FEMechModel.h"
#include "FEPreStrainElastic.h"
#include <FECore/writeplot.h>
#include <FECore/FEDomainParameter.h>
#include <FECore/FEModel.h>
#include "FEDiscreteElasticMaterial.h"
#include "FEDiscreteElasticDomain.h"
#include "FEContinuousElasticDamage.h"
#include <FECore/FEMeshAdaptor.h> // for projectToNodes
#include "FESlidingInterface.h"
#include "FESlidingElasticInterface.h"
#include "FETiedContactSurface.h"
#include "FEReactiveVEMaterialPoint.h"
#include <FECore/FESurface.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FETrussDomain.h>
#include <FECore/FEElement.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotNodeDisplacement::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
	const int dof_X = fem->GetDOFIndex("x");
	const int dof_Y = fem->GetDOFIndex("y");
	const int dof_Z = fem->GetDOFIndex("z");

	writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
		return node.get_vec3d(dof_X, dof_Y, dof_Z);
		});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeVelocity::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
	const int dof_VX = fem->GetDOFIndex("vx");
	const int dof_VY = fem->GetDOFIndex("vy");
	const int dof_VZ = fem->GetDOFIndex("vz");

	writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
		return node.get_vec3d(dof_VX, dof_VY, dof_VZ);
	});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeAcceleration::Save(FEMesh& m, FEDataStream& a)
{
	writeNodalValues<vec3d>(m, a, [](const FENode& node) {
		return node.m_at;
	});
	return true;
}

//-----------------------------------------------------------------------------
//! Store nodal reaction forces
bool FEPlotNodeReactionForces::Save(FEMesh& m, FEDataStream& a)
{
	// NOTE: Currently, for nodes attached to rigid bodies the reaction forces 
	//       will actually be the rigid reaction forces. 
	FEModel& fem = *GetFEModel();
	int dofX = fem.GetDOFIndex("x");
	int dofY = fem.GetDOFIndex("y");
	int dofZ = fem.GetDOFIndex("z");
	if ((dofX >= 0) && (dofY >= 0) && (dofZ >= 0))
	{
		writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
			return node.get_load3(dofX, dofY, dofZ);
			});
		return true;
	}
	else
		return false;
}

//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
// Plot contact gap
bool FEPlotContactGap::Save(FESurface& surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		// NOTE: the sliding surface does not use material points, so we need this little hack. 
		FESlidingSurface* ss = dynamic_cast<FESlidingSurface*>(pcs);
		if (ss)
		{
			for (int i = 0; i < ss->Elements(); ++i)
			{
				FEElement& el = ss->Element(i);
				double g = 0.0;
				for (int j = 0; j < el.Nodes(); ++j)
				{
					double gj = ss->m_data[el.m_lnode[j]].m_gap;
					g += gj;
				}
				g /= el.Nodes();
				a << g;
			}
			return true;
		}

		FETiedContactSurface* ts = dynamic_cast<FETiedContactSurface*>(pcs);
		if (ts)
		{
			for (int i = 0; i < ts->Elements(); ++i)
			{
				FEElement& el = ts->Element(i);
				double g = 0.0;
				for (int j = 0; j < el.Nodes(); ++j)
				{
					double gj = ts->m_data[el.m_lnode[j]].m_gap;
					g += gj;
				}
				g /= el.Nodes();
				a << g;
			}
			return true;
		}

		writeAverageElementValue<double>(surf, a, [=](const FEMaterialPoint& mp) {
			const FEContactMaterialPoint* pt = dynamic_cast<const FEContactMaterialPoint*>(&mp);
			return (pt ? pt->m_gap : 0.0);
			});

		return true;
	}
	else return false;
}

//-----------------------------------------------------------------------------
// Plot vector gap
bool FEPlotVectorGap::Save(FESurface& surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		writeElementValue<vec3d>(surf, a, [=](int nface) {
			vec3d gn;
			pcs->GetVectorGap(nface, gn);
			return gn;
			});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot contact pressure
bool FEPlotContactPressure::Save(FESurface &surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		// NOTE: the sliding surface does not use material points, so we need this little hack. 
		FESlidingSurface* ss = dynamic_cast<FESlidingSurface*>(pcs);
		if (ss)
		{
			for (int i = 0; i < ss->Elements(); ++i)
			{
				FEElement& el = ss->Element(i);
				double Lm = 0.0;
				for (int j = 0; j < el.Nodes(); ++j)
				{
					double Lmj = ss->m_data[el.m_lnode[j]].m_Ln;
					Lm += Lmj;
				}
				Lm /= el.Nodes();
				a << Lm;
			}
			return true;
		}

		writeAverageElementValue<double>(surf, a, [](const FEMaterialPoint& mp) {
			const FEContactMaterialPoint* pt = dynamic_cast<const FEContactMaterialPoint*>(&mp);
			return (pt ? pt->m_Ln : 0.0);
		});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot contact traction
bool FEPlotContactTraction::Save(FESurface &surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		writeElementValue<vec3d>(surf, a, [=](int nface) {
			vec3d tn;
			pcs->GetContactTraction(nface, tn);
			return tn;
		});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot nodal contact gap
bool FEPlotNodalContactGap::Save(FESurface& surf, FEDataStream& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		// NOTE: the sliding surface does not use material points, so we need this little hack. 
		FESlidingSurface* ss = dynamic_cast<FESlidingSurface*>(pcs);
		if (ss)
		{
			for (int i = 0; i < ss->Elements(); ++i)
			{
				FEElement& el = ss->Element(i);
				for (int j = 0; j < el.Nodes(); ++j)
				{
					double gap = ss->m_data[el.m_lnode[j]].m_gap;
					a << gap;
				}
			}
			return true;
		}

		writeNodalProjectedElementValues<double>(surf, a, [](const FEMaterialPoint& mp) {
			const FEContactMaterialPoint* pt = dynamic_cast<const FEContactMaterialPoint*>(&mp);
			return (pt ? pt->m_gap : 0.0);
			});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot nodal vector gap
bool FEPlotNodalVectorGap::Save(FESurface &surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		int NF = pcs->Elements();
		vec3d gn[FEElement::MAX_NODES];
		for (int j = 0; j < NF; ++j)
		{
			FESurfaceElement& el = pcs->Element(j);
			pcs->GetNodalVectorGap(j, gn);

			// store in archive
			int ne = el.Nodes();
			for (int k = 0; k < ne; ++k) a << gn[k];
		}

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot nodal contact pressure
bool FEPlotNodalContactPressure::Save(FESurface &surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		writeNodalProjectedElementValues<double>(surf, a, [](const FEMaterialPoint& mp) {
			const FEContactMaterialPoint* pt = dynamic_cast<const FEContactMaterialPoint*>(&mp);
			return (pt ? pt->m_Ln : 0.0);
			});

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot nodal contact traction
bool FEPlotNodalContactTraction::Save(FESurface &surf, FEDataStream& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		int NF = pcs->Elements();
		vec3d tn[FEElement::MAX_NODES];
		for (int j = 0; j < NF; ++j)
		{
			FESurfaceElement& el = pcs->Element(j);
			pcs->GetNodalContactTraction(j, tn);

			// store in archive
			int ne = el.Nodes();
			for (int k = 0; k < ne; ++k) a << tn[k];
		}

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot surface traction
bool FEPlotSurfaceTraction::Save(FESurface &surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{

		writeElementValue<vec3d>(surf, a, [=](int nface) {
			vec3d tn;
			pcs->GetSurfaceTraction(nface, tn);
			return tn;
			});

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot nodal contact traction
bool FEPlotNodalSurfaceTraction::Save(FESurface &surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		int NF = pcs->Elements();
		vec3d tn[FEElement::MAX_NODES];
		for (int j = 0; j < NF; ++j)
		{
			FESurfaceElement& el = pcs->Element(j);
			pcs->GetNodalSurfaceTraction(j, tn);

			// store in archive
			int ne = el.Nodes();
			for (int k = 0; k < ne; ++k) a << tn[k];
		}

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot stick status
bool FEPlotStickStatus::Save(FESurface& surf, FEDataStream& a)
{
    FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
    if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		writeElementValue<double>(surf, a, [=](int nface) {
			double gn;
			pcs->GetStickStatus(nface, gn);
			return gn;
			});

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotContactForce::Save(FESurface &surf, FEDataStream &a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		vec3d fn = pcs->GetContactForce();
		a << fn;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot contact area
bool FEPlotContactArea::Save(FESurface &surf, FEDataStream& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		double area = pcs->GetContactArea();
		a << area;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Plot contact penalty parameter
bool FEPlotContactPenalty::Save(FESurface& surf, FEDataStream& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		FEFacetSlidingSurface* ps = dynamic_cast<FEFacetSlidingSurface*>(&surf);
		if (ps)
		{
			writeAverageElementValue<double>(surf, a, [](const FEMaterialPoint& mp) {
				const FEFacetSlidingSurface::Data* pt = dynamic_cast<const FEFacetSlidingSurface::Data*>(&mp);
				return (pt ? pt->m_eps : 0);
				});
			return true;
		}

		FESlidingElasticSurface* pse = dynamic_cast<FESlidingElasticSurface*>(&surf);
		if (pse)
		{
			writeAverageElementValue<double>(surf, a, [](const FEMaterialPoint& mp) {
				const FESlidingElasticSurface::Data* pt = dynamic_cast<const FESlidingElasticSurface::Data*>(&mp);
				return (pt ? pt->m_epsn : 0);
				});
			return true;
		}
	}

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotContactStatus::Save(FESurface& surf, FEDataStream& a)
{
	FEFacetSlidingSurface* ps = dynamic_cast<FEFacetSlidingSurface*>(&surf);
	if (ps == nullptr) return false;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = ps->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		int NF = ps->Elements();
		for (int i = 0; i < NF; ++i)
		{
			FESurfaceElement& el = ps->Element(i);
			double nc = 0.0;
			int nint = el.GaussPoints();
			for (int j = 0; j < nint; ++j)
			{
				FEContactMaterialPoint* mp = dynamic_cast<FEContactMaterialPoint*>(el.GetMaterialPoint(j));
				if (mp && mp->m_pme) nc++;
			}

			a << nc;
		}

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotMortarContactGap::Save(FESurface& S, FEDataStream& a)
{
	FEMortarSlidingSurface* ps = dynamic_cast<FEMortarSlidingSurface*>(&S);
	if (ps)
	{
		writeNodalValues<double>(S, a, [=](int i) {
			vec3d vA = ps->m_nu[i];
			vec3d gA = ps->m_gap[i];
			return gA*vA;
		});
		return true;
	}
	else return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEnclosedVolume::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
	writeIntegratedElementValue<double>(surf, a, [=](const FEMaterialPoint& mp) {
		FESurfaceElement& el = static_cast<FESurfaceElement&>(*mp.m_elem);
		int n = mp.m_index;
		vec3d xi = pcs->Local2Global(el, n);
		vec3d g[2];
		pcs->CoBaseVectors(el, n, g);
		return xi*(g[0] ^ g[1]) / 3;
	});
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSurfaceArea::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    writeIntegratedElementValue<double>(surf, a, [=](const FEMaterialPoint& mp) {
        FESurfaceElement& el = static_cast<FESurfaceElement&>(*mp.m_elem);
        int n = mp.m_index;
        vec3d g[2];
        pcs->CoBaseVectors(el, n, g);
        return (g[0] ^ g[1]).norm();
    });
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFacetArea::Save(FESurface& surf, FEDataStream& a)
{
	FESurface* pcs = &surf;
	if (pcs == 0) return false;

	writeElementValue<double>(surf, a, [=](int nface) {
		double A = pcs->CurrentFaceArea(pcs->Element(nface));
		return A;
		});
	return true;
}

//-----------------------------------------------------------------------------
// Plot scalar surface load
bool FEPlotScalarSurfaceLoad::Save(FESurface &surf, FEDataStream& a)
{
    FEModel* fem = GetFEModel();
    int nsl = fem->ModelLoads();
    FESurfaceLoad* psl = nullptr;
    for (int i = 0; i<nsl; ++i)
    {
        psl = dynamic_cast<FESurfaceLoad*>(fem->ModelLoad(i));
        if (psl && (&psl->GetSurface() == &surf)) break;
    }

    if (psl == nullptr) return false;
    
    if (psl->IsActive()) {
        writeAverageElementValue<double>(surf, a, [=](const FEMaterialPoint& mp) {
			const FESurfaceMaterialPoint* pt = dynamic_cast<const FESurfaceMaterialPoint*>(&mp);
			return (pt ? psl->ScalarLoad(*const_cast<FESurfaceMaterialPoint*>(pt)) : 0.0);
        });
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotNetSurfaceReactionForce::Save(FESurface& surf, FEDataStream& a)
{
	// NOTE: Currently, for nodes attached to rigid bodies the reaction forces 
	//       will actually be the rigid reaction forces. 
	FEModel& fem = *GetFEModel();
	int dofX = fem.GetDOFIndex("x");
	int dofY = fem.GetDOFIndex("y");
	int dofZ = fem.GetDOFIndex("z");
	if ((dofX < 0) || (dofY < 0) || (dofZ < 0)) return false;

	vec3d F(0, 0, 0);
	for (int i = 0; i < surf.Nodes(); ++i)
	{
		FENode& n = surf.Node(i);
		vec3d Fi = n.get_load3(dofX, dofY, dofZ);
		F += Fi;
	}

	a << F;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNetSurfaceReactionMoment::Save(FESurface& surf, FEDataStream& a)
{
	// NOTE: Currently, for nodes attached to rigid bodies the reaction forces 
	//       will actually be the rigid reaction forces. 
	FEModel& fem = *GetFEModel();
	int dofX = fem.GetDOFIndex("x");
	int dofY = fem.GetDOFIndex("y");
	int dofZ = fem.GetDOFIndex("z");
	if ((dofX < 0) || (dofY < 0) || (dofZ < 0)) return false;

	vec3d M(0, 0, 0);
	for (int i = 0; i < surf.Nodes(); ++i)
	{
		FENode& n = surf.Node(i);
		vec3d Fi = n.get_load3(dofX, dofY, dofZ);
		vec3d r = n.m_rt;
		vec3d Mi = r ^ Fi;
		M += Mi;
	}

	a << M;

	return true;
}


//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotElementVelocity::Save(FEDomain &dom, FEDataStream& a)
{
	FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
	if ((pme == 0) || pme->IsRigid()) return false;

	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		return pt.m_v;
	});

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotElementAcceleration::Save(FEDomain &dom, FEDataStream& a)
{
	FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
    if ((pme == 0) || pme->IsRigid()) return false;

	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		return pt.m_a;
	});

	return true;
}

//=============================================================================
class FEStress
{
public:
	mat3ds operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
		return (pt ? pt->m_s : mat3ds(0));
	}
};

//-----------------------------------------------------------------------------
//! Store the average stresses for each element. 
bool FEPlotElementStress::Save(FEDomain& dom, FEDataStream& a)
{
	FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
	if ((pme == 0) || pme->IsRigid()) return false;

	FEDomainParameter* var = pme->FindDomainParameter("stress");
	if (var == nullptr) return false;

	writeAverageElementValue<mat3ds>(dom, a, var);

	return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element. 
bool FEPlotElementPK2Stress::Save(FEDomain& dom, FEDataStream& a)
{
	FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
	if ((pme == 0) || pme->IsRigid()) return false;

	writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
		const FEElasticMaterialPoint& ep = *mp.ExtractData< FEElasticMaterialPoint>();
		mat3ds s = ep.m_s;
		mat3ds S = ep.pull_back(s);
		return S;
		});

	return true;
}

//-----------------------------------------------------------------------------
//! Store the average PK1 stress for each element. 
bool FEPlotElementPK1Stress::Save(FEDomain& dom, FEDataStream& a)
{
	FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
	if ((pme == 0) || pme->IsRigid()) return false;

	writeAverageElementValue<mat3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEElasticMaterialPoint& ep = *mp.ExtractData< FEElasticMaterialPoint>();
		mat3d  F = ep.m_F;
		double J = F.det();
		mat3ds s = ep.m_s;	// Cauchy stress

		mat3d P = (s * F.transinv()) * J;
		return P;
		});

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRStresses::Save(FEDomain& dom, FEDataStream& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	writeSPRElementValueMat3ds(sd, a, FEStress());

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRLinearStresses::Save(FEDomain& dom, FEDataStream& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	writeSPRElementValueMat3ds(sd, a, FEStress(), 1);

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodalStresses::Save(FEDomain& dom, FEDataStream& a)
{
	writeNodalProjectedElementValues<mat3ds>(dom, a, FEStress());
	return true;
}

//=============================================================================
FEPlotElementMixtureStress::FEPlotElementMixtureStress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) 
{
	m_comp = -1;
	SetUnits(UNIT_PRESSURE);
}

bool FEPlotElementMixtureStress::SetFilter(const char* szfilter)
{
	sscanf(szfilter, "solid[%d]", &m_comp);
	return true;
}

bool FEPlotElementMixtureStress::Save(FEDomain& dom, FEDataStream& a)
{
	// make sure we start from the elastic component
	FEElasticMaterial* pmat = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pmat == nullptr) return false;

	// make sure this is a mixture
	FEElasticMixture* pmm = dynamic_cast<FEElasticMixture*>(pmat);
	FEUncoupledElasticMixture* pum = dynamic_cast<FEUncoupledElasticMixture*>(pmat);
	if ((pmm == nullptr) && (pum == nullptr)) return false;

	// get the mixture component
	if (m_comp < 0) return false;

	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		mat3ds savg; savg.zero();
		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
			if (mmp && (m_comp < mmp->Components()))
			{
				FEElasticMaterialPoint& ep = *mmp->GetPointData(m_comp)->ExtractData<FEElasticMaterialPoint>();
				savg += ep.m_s;
			}
		}
		savg /= (double)el.GaussPoints();

		a << savg;
	}

	return true;
}

//=============================================================================
//! Store the uncoupled pressure for each element.
bool FEPlotElementUncoupledPressure::Save(FEDomain& dom, FEDataStream& a)
{
	FEUncoupledMaterial* pmu = dom.GetMaterial()->ExtractProperty<FEUncoupledMaterial>();
    if (pmu == 0) return false;
   
    // write element data
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
		if (pt == 0) return 0.0;
		return -pt->m_p;   // use negative sign to get positive pressure in compression
	});
    
    return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average Cauchy stress for each element. 
bool FEPlotElementsnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
	if ((pme == 0) || pme->IsRigid()) return false;

	writeAverageElementValue<mat3ds, double>(dom, a, FEStress(), [](const mat3ds& s) { return sqrt(s.dotdot(s)); });
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average elasticity for each element.

class FEElementElasticity
{
public:
	FEElementElasticity(FESolidMaterial* pm) : m_mat(pm) {}
	tens4ds operator()(const FEMaterialPoint& mp)
	{
		return m_mat->Tangent(const_cast<FEMaterialPoint&>(mp));
	}
private:
	FESolidMaterial*	m_mat;
};

bool FEPlotElementElasticity::Save(FEDomain& dom, FEDataStream& a)
{
    FESolidMaterial* pme = dom.GetMaterial()->ExtractProperty<FESolidMaterial>();
    if ((pme == 0) || pme->IsRigid()) return false;

	writeAverageElementValue<tens4ds>(dom, a, FEElementElasticity(pme));
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average deviatoric elasticity for each element.

class FEElementDevElasticity
{
public:
    FEElementDevElasticity(FEUncoupledMaterial* pm) : m_mat(pm) {}
    tens4ds operator()(const FEMaterialPoint& mp)
    {
        return m_mat->DevTangent(const_cast<FEMaterialPoint&>(mp));
    }
private:
    FEUncoupledMaterial*    m_mat;
};

bool FEPlotElementDevElasticity::Save(FEDomain& dom, FEDataStream& a)
{
    FEUncoupledMaterial* pme = dom.GetMaterial()->ExtractProperty<FEUncoupledMaterial>();
    if ((pme == 0) || pme->IsRigid()) return false;
    
    writeAverageElementValue<tens4ds>(dom, a, FEElementDevElasticity(pme));
    return true;
}

//-----------------------------------------------------------------------------
class FEStrainEnergy
{
public:
	FEStrainEnergy(FEElasticMaterial* pm) : m_mat(pm) {}
	double operator()(const FEMaterialPoint& mp)
	{
		return m_mat->StrainEnergyDensity(const_cast<FEMaterialPoint&>(mp));
	}
private:
	FEElasticMaterial*	m_mat;
};

bool FEPlotStrainEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == 0) return false;
    
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FEStrainEnergy W(pme);
		writeAverageElementValue<double>(dom, a, W);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
class FEDevStrainEnergy
{
public:
	FEDevStrainEnergy(FEUncoupledMaterial* pm) : m_mat(pm) {}
	double operator()(const FEMaterialPoint& mp)
	{
		return m_mat->DevStrainEnergyDensity(const_cast<FEMaterialPoint&>(mp));
	}
private:
	FEUncoupledMaterial*	m_mat;
};

bool FEPlotDevStrainEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEUncoupledMaterial* pmu = dom.GetMaterial()->ExtractProperty<FEUncoupledMaterial>();
    if (pmu == 0) return false;
    
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FEDevStrainEnergy devW(pmu);
		writeAverageElementValue<double>(dom, a, devW);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
class FESpecificStrainEnergy
{
public:
	double operator()(const FEMaterialPoint& mp)
	{
		const FERemodelingMaterialPoint* rpt = mp.ExtractData<FERemodelingMaterialPoint>();
		return (rpt ? rpt->m_sed / rpt->m_rhor : 0.0);
	}
};

bool FEPlotSpecificStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESpecificStrainEnergy E;
	writeAverageElementValue<double>(dom, a, E);

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotKineticEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    const int dof_VX = GetFEModel()->GetDOFIndex("vx");
    const int dof_VY = GetFEModel()->GetDOFIndex("vy");
    const int dof_VZ = GetFEModel()->GetDOFIndex("vz");
    const int dof_VU = GetFEModel()->GetDOFIndex("vu");
    const int dof_VV = GetFEModel()->GetDOFIndex("vv");
    const int dof_VW = GetFEModel()->GetDOFIndex("vw");
    
    FEMesh& mesh = *dom.GetMesh();
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            double *H;
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d vt[FEElement::MAX_NODES];
            vec3d vn[FEElement::MAX_NODES];
            for (int j=0; j<el.Nodes(); ++j) {
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                H = el.H(j);
                vn[j] = vec3d(0, 0, 0);
                for (int k=0; k<el.Nodes(); ++k)
                    vn[j] += vt[k]*H[k];
            }
            
            // integrate kinetic energy
            double ew = 0;
            double V = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double detJ = bd.detJ0(el, j)*gw[j];
                V += detJ;
                ew += vn[j]*vn[j]*(pme->Density(mp)/2*detJ);
            }
            
            a << ew/V;
        }
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
        for (int i=0; i<bd->Elements(); ++i)
        {
            FEShellElement& el = bd->Element(i);
            double* gw = el.GaussWeights();

			double ew = 0.0;
			if ((dof_VU >= 0) && (dof_VV >= 0) && (dof_VW >= 0))
			{
				// get nodal velocities
				vec3d vt[FEElement::MAX_NODES];
				vec3d wt[FEElement::MAX_NODES];
				vec3d vn[FEElement::MAX_NODES];
				for (int j = 0; j < el.Nodes(); ++j) {
					vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
					wt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VU, dof_VV, dof_VW);
				}

				// evaluate velocities at integration points
				for (int j = 0; j < el.GaussPoints(); ++j)
					vn[j] = bd->evaluate(el, vt, wt, j);

				// integrate kinetic energy
				double ew = 0;
				double V = 0;
				for (int j = 0; j < el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);

					double detJ = bd->detJ0(el, j) * gw[j];
					V += detJ;
					ew += vn[j] * vn[j] * (pme->Density(mp) / 2 * detJ);
				}

				// normalize by volume
				ew /= V;
			}
            a << ew;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
// TODO: Should I call the density for remodeling materials something else? 
//       Or maybe the FEElasticMaterialPoint should define a density parameter
//       that will be updated by the materials to define the current density?

class FEDensity
{
public:
	FEDensity(FEElasticMaterial* pm) : m_mat(pm) {}
	double operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
		double J = ep.m_F.det();
		return m_mat->Density(const_cast<FEMaterialPoint&>(mp)) / J;
	}
private:
	FEElasticMaterial*	m_mat;
};

class FERemodelingDensity
{
public:
	double operator()(const FEMaterialPoint& mp)
	{
		const FERemodelingMaterialPoint* pt = (mp.ExtractData<FERemodelingMaterialPoint>());
		return (pt ? pt->m_rhor : 0.0);
	}
};

bool FEPlotDensity::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		FEElasticMaterial* em = dynamic_cast<FEElasticMaterial*>(bd.GetMaterial());
		if (em == 0) return false;

		FERemodelingElasticMaterial* rm = dynamic_cast<FERemodelingElasticMaterial*>(em);
		if (rm)
		{
			FERemodelingDensity dens;
			writeAverageElementValue<double>(dom, a, dens);
			return true;
		}
		else
		{
			FEDensity dens(em);
			writeAverageElementValue<double>(dom, a, dens);
			return true;
		}
	}
	else if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEElasticMaterial* em = dynamic_cast<FEElasticMaterial*>(dom.GetMaterial());
		if (em == 0) return false;
		FEDensity dens(em);
		writeAverageElementValue<double>(dom, a, dens);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElementStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<double>(bd, a, FEStrainEnergy(pme));
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
		writeIntegratedElementValue(*bd, a, FEStrainEnergy(pme));
		return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
// integrated element kinetic energy
class FEKineticEnergyDensity
{
public:
	FEKineticEnergyDensity(FEElasticMaterial* pm) : m_mat(pm) {}
	double operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
		return 0.5*(ep.m_v*ep.m_v)*m_mat->Density(const_cast<FEMaterialPoint&>(mp));
	}
private:
	FEElasticMaterial*	m_mat;
};

bool FEPlotElementKineticEnergy::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<double>(bd, a, FEKineticEnergyDensity(pme));
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
		writeIntegratedElementValue(*bd, a, FEKineticEnergyDensity(pme));
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElementCenterOfMass::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		for (int i = 0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			double* gw = el.GaussWeights();

			// integrate zeroth and first mass moments
			vec3d ew = vec3d(0, 0, 0);
			double m = 0;
			for (int j = 0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				double detJ = bd.detJ0(el, j)*gw[j];
				ew += mp.m_rt*(pme->Density(mp)*detJ);
				m += pme->Density(mp)*detJ;
			}

			a << ew / m;
		}
		return true;
	}
	else if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
		if (bd == 0) return false;
		for (int i = 0; i<bd->Elements(); ++i)
		{
			FEShellElement& el = bd->Element(i);
			double* gw = el.GaussWeights();

			// integrate zeroth and first mass moments
			vec3d ew = vec3d(0, 0, 0);
			double m = 0;
			for (int j = 0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& pt = *el.GetMaterialPoint(j);
				double detJ = bd->detJ0(el, j)*gw[j];
				ew += pt.m_rt*(pme->Density(pt)*detJ);
				m += pme->Density(pt)*detJ;
			}

			a << ew / m;
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
class FEElementLinearMomentum
{
public:
	FEElementLinearMomentum(FEElasticMaterial* pm) : m_mat(pm) {}
	vec3d operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		return pt.m_v*m_mat->Density(const_cast<FEMaterialPoint&>(mp));
	}

private:
	FEElasticMaterial*	m_mat;
};

bool FEPlotElementLinearMomentum::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<vec3d>(bd, a, FEElementLinearMomentum(pme));
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
		writeIntegratedElementValue(*bd, a, FEElementLinearMomentum(pme));
		return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
// Integrated element angular momentum
class FEElementAngularMomentum
{
public:
	FEElementAngularMomentum(FEElasticMaterial* pm) : m_mat(pm) {}
	vec3d operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		return (mp.m_rt ^ pt.m_v)*m_mat->Density(const_cast<FEMaterialPoint&>(mp));
	}

private:
	FEElasticMaterial*	m_mat;
};

bool FEPlotElementAngularMomentum::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<vec3d>(bd, a, FEElementAngularMomentum(pme));
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
		writeIntegratedElementValue(*bd, a, FEElementAngularMomentum(pme));
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
class FEElementStressPower
{
public:
	double operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
		return  ep.m_s.dotdot(ep.m_L.sym())*ep.m_J;
	}
};

bool FEPlotElementStressPower::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<double>(bd, a, FEElementStressPower());
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
		writeIntegratedElementValue(*bd, a, FEElementStressPower());
		return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
class FECurrentElementStrainEnergy
{
public:
	double operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
		return ep.m_Wt;
	}
};

bool FEPlotCurrentElementStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<double>(bd, a, FECurrentElementStrainEnergy());
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
		writeIntegratedElementValue(*bd, a, FECurrentElementStrainEnergy());
		return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentElementKineticEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    const int dof_VX = GetFEModel()->GetDOFIndex("vx");
    const int dof_VY = GetFEModel()->GetDOFIndex("vy");
    const int dof_VZ = GetFEModel()->GetDOFIndex("vz");
    const int dof_VU = GetFEModel()->GetDOFIndex("vu");
    const int dof_VV = GetFEModel()->GetDOFIndex("vv");
    const int dof_VW = GetFEModel()->GetDOFIndex("vw");
    
    FEMesh& mesh = *dom.GetMesh();
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    const int NELN = FEElement::MAX_NODES;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d vt[NELN], vn[NELN];
            for (int j=0; j<el.Nodes(); ++j)
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
                vn[j] = el.Evaluate(vt, j);
            
            // integrate kinetic energy
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);

                double detJ = bd.detJ0(el, j)*gw[j]* pme->Density(mp)/2;
                ew += vn[j]*vn[j]*detJ;
            }
            
            a << ew;
        }
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
        for (int i=0; i<bd->Elements(); ++i)
        {
            FEShellElement& el = bd->Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d vt[NELN], wt[NELN], vn[NELN];
            for (int j=0; j<el.Nodes(); ++j) {
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
                wt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VU, dof_VV, dof_VW);
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
                vn[j] = bd->evaluate(el, vt, wt, j);
            
            // integrate kinetic energy
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);

                double detJ = bd->detJ0(el, j)*gw[j]* pme->Density(mp)/2;
                ew += vn[j]*vn[j]*detJ;
            }
            
            a << ew;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentElementCenterOfMass::Save(FEDomain &dom, FEDataStream& a)
{
    FEMesh& mesh = *dom.GetMesh();
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    const int NELN = FEElement::MAX_NODES;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal positions and velocities
            vec3d rt[NELN], rn[NELN];
            for (int j=0; j<el.Nodes(); ++j)
                rt[j] = mesh.Node(el.m_node[j]).m_rt;
            
            // evaluate positions at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
                rn[j] = el.Evaluate(rt, j);
            
            // integrate zeroth and first mass moment
            double ez = 0;
            vec3d ef = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);

                double detJ = bd.detJ0(el, j)*gw[j]* pme->Density(mp);
                ez += detJ;
                ef += rn[j]*detJ;
            }
            
            a << ef/ez;
        }
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
        for (int i=0; i<bd->Elements(); ++i)
        {
            FEShellElement& el = bd->Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d rt[NELN], st[NELN], rn[NELN];
            for (int j=0; j<el.Nodes(); ++j) {
                rt[j] = mesh.Node(el.m_node[j]).m_rt;
                st[j] = mesh.Node(el.m_node[j]).st();
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
                rn[j] = bd->evaluate(el, rt, st, j);
            
            // integrate zeroth and first mass moment
            double ez = 0;
            vec3d ef = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);

                double detJ = bd->detJ0(el, j)*gw[j]* pme->Density(mp);
                ez += detJ;
                ef += rn[j]*detJ;
            }
            
            a << ef/ez;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentElementLinearMomentum::Save(FEDomain &dom, FEDataStream& a)
{
    const int dof_VX = GetFEModel()->GetDOFIndex("vx");
    const int dof_VY = GetFEModel()->GetDOFIndex("vy");
    const int dof_VZ = GetFEModel()->GetDOFIndex("vz");
    const int dof_VU = GetFEModel()->GetDOFIndex("vu");
    const int dof_VV = GetFEModel()->GetDOFIndex("vv");
    const int dof_VW = GetFEModel()->GetDOFIndex("vw");
    
    FEMesh& mesh = *dom.GetMesh();
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    const int NELN = FEElement::MAX_NODES;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d vt[NELN], vn[NELN];
            for (int j=0; j<el.Nodes(); ++j) {
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
                vn[j] = el.Evaluate(vt, j);
            
            // integrate linear momentum
            vec3d ew = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double detJ = bd.detJ0(el, j)*gw[j];
                ew += vn[j]*(pme->Density(mp)*detJ);
            }
            
            a << ew;
        }
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
        for (int i=0; i<bd->Elements(); ++i)
        {
            FEShellElement& el = bd->Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d vt[NELN], wt[NELN], vn[NELN];
            for (int j=0; j<el.Nodes(); ++j) {
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
                wt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VU, dof_VV, dof_VW);
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j)
                vn[j] = bd->evaluate(el, vt, wt, j);
            
            // integrate linear momentum
            vec3d ew = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double detJ = bd->detJ0(el, j)*gw[j];
                ew += vn[j]*(pme->Density(mp)*detJ);
            }
            
            a << ew;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentElementAngularMomentum::Save(FEDomain &dom, FEDataStream& a)
{
    const int dof_VX = GetFEModel()->GetDOFIndex("vx");
    const int dof_VY = GetFEModel()->GetDOFIndex("vy");
    const int dof_VZ = GetFEModel()->GetDOFIndex("vz");
    const int dof_SVX = GetFEModel()->GetDOFIndex("svx");
    const int dof_SVY = GetFEModel()->GetDOFIndex("svy");
    const int dof_SVZ = GetFEModel()->GetDOFIndex("svz");
    
    FEMesh& mesh = *dom.GetMesh();
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

    const int NELN = FEElement::MAX_NODES;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal positions and velocities
            vec3d rt[NELN], rn[NELN];
            vec3d vt[NELN], vn[NELN];
            for (int j=0; j<el.Nodes(); ++j) {
                rt[j] = mesh.Node(el.m_node[j]).m_rt;
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j) {
                rn[j] = el.Evaluate(rt, j);
                vn[j] = el.Evaluate(vt, j);
            }
            
            // integrate angular momentum
            vec3d ew = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double detJ = bd.detJ0(el, j)*gw[j];
                ew += (rn[j] ^ vn[j])*(pme->Density(mp)*detJ);
            }
            
            a << ew;
        }
        return true;
    }
    else if (dom.Class() == FE_DOMAIN_SHELL)
    {
        FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
        if (bd == 0) return false;
        for (int i=0; i<bd->Elements(); ++i)
        {
            FEShellElement& el = bd->Element(i);
            double* gw = el.GaussWeights();
            
            // get nodal velocities
            vec3d rt[NELN], st[NELN], rn[NELN];
            vec3d vt[NELN], wt[NELN], vn[NELN];
            for (int j=0; j<el.Nodes(); ++j) {
                rt[j] = mesh.Node(el.m_node[j]).m_rt;
                st[j] = mesh.Node(el.m_node[j]).st();
                vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_VX, dof_VY, dof_VZ);
                wt[j] = mesh.Node(el.m_node[j]).get_vec3d(dof_SVX, dof_SVY, dof_SVZ);
            }
            
            // evaluate velocities at integration points
            for (int j=0; j<el.GaussPoints(); ++j) {
                rn[j] = bd->evaluate(el, rt, st, j);
                vn[j] = bd->evaluate(el, vt, wt, j);
            }
            
            // integrate angular momentum
            vec3d ew = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double detJ = bd->detJ0(el, j)*gw[j];
                ew += (rn[j] ^ vn[j])*(pme->Density(mp)*detJ);
            }
            
            a << ew;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotRelativeVolume::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
			const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
			return (pt ? pt->m_J : 0.0);
			});
	}
	else if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEShellDomain* sd = dynamic_cast<FEShellDomain*>(&dom); assert(sd);

		// a filter to get J from a strain tensor
		auto getJfromE = [](const mat3ds& E) {
			mat3ds C = mat3dd(1) + E * 2;
			return sqrt(C.det());
		};

		FEShellDomainNew* newsd = dynamic_cast<FEShellDomainNew*>(sd);
		FEElasticEASShellDomain* easd = dynamic_cast<FEElasticEASShellDomain*>(newsd);
		FEElasticANSShellDomain* ansd = dynamic_cast<FEElasticANSShellDomain*>(newsd);
		if (easd || ansd) {
			writeAverageElementValue<mat3ds, double>(dom, a, [](FEElement& el, int ip) {
				FEShellElementNew& se = static_cast<FEShellElementNew&>(el);
				return se.m_E[ip];
				}, getJfromE);
		}
		else {
			writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
				const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
				return (pt ? pt->m_J : 0.0);
				});
		}
		return true;
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotShellRelativeVolume::Save(FEDomain& dom, FEDataStream& a)
{
	FEShellDomain* sd = dynamic_cast<FEShellDomain*>(&dom); assert(sd);
	if (sd == nullptr) return false;

	// a filter to get J from a strain tensor
	auto getJfromE = [](const mat3ds& E) {
		mat3ds C = mat3dd(1) + E * 2;
		return sqrt(C.det());
	};

	FEShellDomainNew* newsd = dynamic_cast<FEShellDomainNew*>(sd);
	FEElasticEASShellDomain* easd = dynamic_cast<FEElasticEASShellDomain*>(newsd);
	FEElasticANSShellDomain* ansd = dynamic_cast<FEElasticANSShellDomain*>(newsd);
	if (easd || ansd) {
		writeAverageElementValue<mat3ds, double>(dom, a, [](FEElement& el, int ip) {
			FEShellElementNew& se = static_cast<FEShellElementNew&>(el);
			return se.m_E[ip];
			}, getJfromE);
	}
	else {
		writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
			const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
			return (pt ? pt->m_J : 0.0);
			});
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFiberStretch::SetFilter(const char* szfilter)
{
	m_matComp = szfilter;
	return true;
}

bool FEPlotFiberStretch::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (m_matComp.empty() == false)
	{
		pme = dynamic_cast<FEElasticMaterial*>(pme->GetProperty(m_matComp.c_str()));
		if (pme == nullptr) return false;
	}

	// get the fiber property
	FEVec3dValuator* vec = dynamic_cast<FEVec3dValuator*>(pme->GetProperty("fiber"));
	if (vec == 0) return false;

	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	writeAverageElementValue<double>(dom, a, [&](const FEMaterialPoint& mp) -> double { 
		const FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
		mat3d Q = pme->GetLocalCS(mp);
		vec3d a0 = vec->unitVector(mp);
		vec3d ar = Q * a0;
		mat3d F = ep.m_F;
		vec3d a = F*ar;
		return a.norm(); 
	});

	return true;
}

//-----------------------------------------------------------------------------
class FEFiberVector
{
public:
	FEFiberVector(FEMaterial* pm, FEVec3dValuator& vec) : m_pm(pm), m_vec(vec) {}
	vec3d operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint* pt = mp.ExtractData<const FEElasticMaterialPoint>();
		if (pt)
		{
			mat3d Q = m_pm->GetLocalCS(mp);
			mat3d F = pt->m_F;
			vec3d a0 = m_vec.unitVector(mp);
			vec3d ar = Q * a0;
			vec3d a = F * ar; a.unit();
			return a;
		}
		else
			return m_vec(mp);
	}
private:
	FEMaterial*		m_pm;
	FEVec3dValuator&	m_vec;
};


FEPlotFiberVector::FEPlotFiberVector(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM) 
{
}

bool FEPlotFiberVector::SetFilter(const char* szfilter)
{
	m_matComp = szfilter;
	return true;
}

bool FEPlotFiberVector::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (m_matComp.empty() == false)
	{
		pme = dynamic_cast<FEElasticMaterial*>(pme->GetProperty(m_matComp.c_str()));
		if (pme == nullptr) return false;
	}

	// get the fiber property
	FEVec3dValuator* vec = dynamic_cast<FEVec3dValuator*>(pme->GetProperty("fiber"));
	if (vec == 0) return false;

	writeAverageElementValue<vec3d, vec3d>(dom, a, FEFiberVector(pme, *vec), [](const vec3d& r) -> vec3d { vec3d n(r); n.unit(); return n; });

	return true;
}

//-----------------------------------------------------------------------------

bool FEPlotMaterialAxes::SetFilter(const char* szfilter)
{
	m_matComp = szfilter;
	return true;
}

bool FEPlotMaterialAxes::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (m_matComp.empty() == false)
	{
		pme = dynamic_cast<FEElasticMaterial*>(pme->GetProperty(m_matComp.c_str()));
		if (pme == nullptr) return false;
	}

	int BE = dom.Elements();
	for (int i = 0; i<BE; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		// I cannot average the material axes since the average may not be orthogonal
		// Until I find a better option, I'll just export the first integration point.
		FEMaterialPoint& mp = *el.GetMaterialPoint(0);
		mat3d Q = pme->GetLocalCS(mp);
		a << Q;
	}
	return true;
}

//-----------------------------------------------------------------------------
// TODO: The factor Jm13 is not used. This doesn't look correct
class FEDevFiberStretch
{
public:
	FEDevFiberStretch(FEElasticMaterial* mat) : m_mat(mat) {}
public:
	double operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

		mat3d Q = m_mat->GetLocalCS(mp);

		// get the material fiber axis
		vec3d r0 = Q.col(0);

		// apply deformation
		vec3d r = pt.m_F*r0;

		// calculate the deviatoric fiber stretch
		double lam = r.norm();
		return lam;
	}

private:
	FEElasticMaterial*	m_mat;
};

bool FEPlotDevFiberStretch::SetFilter(const char* szfilter)
{
	m_matComp = szfilter;
	return true;
}

bool FEPlotDevFiberStretch::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (m_matComp.empty() == false)
	{
		pme = dynamic_cast<FEElasticMaterial*>(pme->GetProperty(m_matComp.c_str()));
		if (pme == nullptr) return false;
	}

	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FEDevFiberStretch lam(pme);
	writeAverageElementValue<double>(dom, a, lam);
	return true;
}


//=============================================================================
// Principal components of stress

class FEPrincStresses
{
public:
	mat3dd operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
		const mat3ds& s = ep.m_s;
		double l[3];
		s.exact_eigen(l);
		return mat3dd(l[0], l[1], l[2]);
	}
};

bool FEPlotSPRPrincStresses::Save(FEDomain& dom, FEDataStream& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	writeSPRElementValueMat3dd(sd, a, FEPrincStresses());

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotDeformationGradient::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	writeAverageElementValue<mat3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		return pt.m_F;
	});

	return true;
}

//=============================================================================
//! Store the average Lagrangian strain
class FELagrangeStrain
{
public:
	mat3ds operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
		if (pt == 0) return mat3ds(0, 0, 0, 0, 0, 0);

		mat3d C = pt->RightCauchyGreen();
		mat3ds E = ((C - mat3dd(1.0))*0.5).sym();
		return E;
	}
};

//-----------------------------------------------------------------------------
bool FEPlotLagrangeStrain::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
			const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			return pt.Strain();
			});
	}
	else if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEShellDomain* sd = dynamic_cast<FEShellDomain*>(&dom); assert(sd);

		FEElasticEASShellDomain* easd = dynamic_cast<FEElasticEASShellDomain*>(&dom);
		FEElasticANSShellDomain* ansd = dynamic_cast<FEElasticANSShellDomain*>(&dom);
		if (easd || ansd)
		{
			writeAverageElementValue<mat3ds>(dom, a, [](FEElement& el, int ip) {
				FEShellElementNew& se = static_cast<FEShellElementNew&>(el);
				return se.m_E[ip];
				});
		}
		else
		{
			writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
				const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				return pt.Strain();
				});
		}
	}
	else return false;

	return true;
}

bool FEPlotShellStrain::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	FEShellDomain* sd = dynamic_cast<FEShellDomain*>(&dom);
	if (sd == nullptr) return false;

	FEElasticEASShellDomain* easd = dynamic_cast<FEElasticEASShellDomain*>(&dom);
	FEElasticANSShellDomain* ansd = dynamic_cast<FEElasticANSShellDomain*>(&dom);
	if (easd || ansd)
	{
		writeAverageElementValue<mat3ds>(dom, a, [](FEElement& el, int ip) {
			FEShellElementNew& se = static_cast<FEShellElementNew&>(el);
			return se.m_E[ip];
			});
	}
	else
	{
		writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
			const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			return pt.Strain();
			});
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotInfStrain::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
			const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

			// displacement tensor
			mat3d U = pt.m_F - mat3dd(1.0);

			// evaluate small strain tensor eij = 0.5*(Uij + Uji)
			mat3ds e = U.sym();

			return e;
			});
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRLagrangeStrain::Save(FEDomain& dom, FEDataStream& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	writeSPRElementValueMat3ds(sd, a, FELagrangeStrain());
	return true;
}

//=============================================================================
//! Store the average right stretch
class FERightStretch
{
public:
    mat3ds operator()(const FEMaterialPoint& mp)
    {
        const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
        if (pt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
            
        return pt->RightStretch();
    }
};

//-----------------------------------------------------------------------------
bool FEPlotRightStretch::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    writeAverageElementValue<mat3ds>(dom, a, FERightStretch());
    return true;
}

//=============================================================================
//! Store the average right stretch
class FELeftStretch
{
public:
    mat3ds operator()(const FEMaterialPoint& mp)
    {
        const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
        if (pt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
            
            return pt->LeftStretch();
            }
};

//-----------------------------------------------------------------------------
bool FEPlotLeftStretch::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    writeAverageElementValue<mat3ds>(dom, a, FELeftStretch());
    return true;
}

//=============================================================================
//! Store the average right Hencky
class FERightHencky
{
public:
    mat3ds operator()(const FEMaterialPoint& mp)
    {
        const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
        if (pt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
            
        return pt->RightHencky();
    }
};

//-----------------------------------------------------------------------------
bool FEPlotRightHencky::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    writeAverageElementValue<mat3ds>(dom, a, FERightHencky());
    return true;
}

//=============================================================================
//! Store the average left Hencky
class FELeftHencky
{
public:
    mat3ds operator()(const FEMaterialPoint& mp)
    {
        const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
        if (pt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
            
        return pt->LeftHencky();
    }
};

//-----------------------------------------------------------------------------
bool FEPlotLeftHencky::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    writeAverageElementValue<mat3ds>(dom, a, FELeftHencky());
    return true;
}

//=============================================================================
//! Store the average rate of deformation
class FERateOfDeformation
{
public:
    mat3ds operator()(const FEMaterialPoint& mp)
    {
        const FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
        if (pt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
            
            return pt->RateOfDeformation();
            }
};

//-----------------------------------------------------------------------------
bool FEPlotRateOfDeformation::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == nullptr) return false;
    writeAverageElementValue<mat3ds>(dom, a, FERateOfDeformation());
    return true;
}

//-----------------------------------------------------------------------------
//! Store shell thicknesses
bool FEPlotShellThickness::Save(FEDomain &dom, FEDataStream &a)
{
	if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
		int NS = sd.Elements();
		for (int i=0; i<NS; ++i)
		{	
			FEShellElement& e = sd.Element(i);
			int n = e.Nodes();
			for (int j=0; j<n; ++j) a << e.m_ht[j];
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store shell directors
bool FEPlotShellDirector::Save(FEDomain &dom, FEDataStream &a)
{
    const int dof_X = GetFEModel()->GetDOFIndex("x");
    const int dof_Y = GetFEModel()->GetDOFIndex("y");
    const int dof_Z = GetFEModel()->GetDOFIndex("z");
	const int dof_U = GetFEModel()->GetDOFIndex("u");
	const int dof_V = GetFEModel()->GetDOFIndex("v");
	const int dof_W = GetFEModel()->GetDOFIndex("w");
    const int dof_SX = GetFEModel()->GetDOFIndex("sx");
    const int dof_SY = GetFEModel()->GetDOFIndex("sy");
    const int dof_SZ = GetFEModel()->GetDOFIndex("sz");
	if (dom.Class() == FE_DOMAIN_SHELL)
	{
		if (dynamic_cast<FEElasticShellDomainOld*>(&dom))
		{
			FEShellDomainOld& sd = static_cast<FEShellDomainOld&>(dom);
			int NS = sd.Elements();
			FEMesh& mesh = *sd.GetMesh();
			for (int i = 0; i<NS; ++i)
			{
				FEShellElementOld& e = sd.ShellElement(i);
				int n = e.Nodes();
				for (int j = 0; j<n; ++j)
				{
					FENode& nj = mesh.Node(e.m_node[j]);
					vec3d D = e.m_D0[j] + nj.get_vec3d(dof_U, dof_V, dof_W);
					a << D;
				}
			}
			return true;
		}
        else if (dynamic_cast<FESSIShellDomain*>(&dom))
        {
            FESSIShellDomain* bd = dynamic_cast<FESSIShellDomain*>(&dom);
            int NS = bd->Elements();
            FEMesh& mesh = *bd->GetMesh();
            for (int i=0; i<NS; ++i)
            {
                FEShellElement& e = bd->Element(i);
                int n = e.Nodes();
                for (int j=0; j<n; ++j)
                {
                    FENode& nj = mesh.Node(e.m_node[j]);
                    vec3d D;
                    if (bd->m_bnodalnormals) {
                        D = nj.m_d0;
                    }
                    else {
                        D = e.m_d0[j];
                    }
                    D += nj.get_vec3d(dof_X, dof_Y, dof_Z) - nj.get_vec3d(dof_SX, dof_SY, dof_SZ);
                    a << D;
                }
            }
            return true;
        }
	}
	return false;
}

//=============================================================================
FEPlotDamage::FEPlotDamage(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
    m_comp = -1;
}

//-----------------------------------------------------------------------------
bool FEPlotDamage::SetFilter(const char* szfilter)
{
    sscanf(szfilter, "solid[%d]", &m_comp);
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotDamage::Save(FEDomain &dom, FEDataStream& a)
{
    if (m_comp == -1) {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
            const FEDamageMaterialPoint* ppd = mp.ExtractData<FEDamageMaterialPoint>();
            double D = 0.0;
            if (ppd) D = (float)ppd->m_D;
            return D;
        });
        return true;
    }
    else {
        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);

            float D = 0.0;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
                if (mmp && (m_comp < mmp->Components()))
                {
                    FEDamageMaterialPoint* dp = mmp->GetPointData(m_comp)->ExtractData<FEDamageMaterialPoint>();
                    if (dp) D += dp->m_D;
                }
            }
            D /= (float)el.GaussPoints();

            a << D;
        }
        return true;
    }
}

//=============================================================================
FEPlotIntactBondFraction::FEPlotIntactBondFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
    m_comp = -1;
}

//-----------------------------------------------------------------------------
bool FEPlotIntactBondFraction::SetFilter(const char* szfilter)
{
    sscanf(szfilter, "solid[%d]", &m_comp);
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotIntactBondFraction::Save(FEDomain &dom, FEDataStream& a)
{
    if (m_comp == -1) {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& pt) {
			FEMaterialPoint& mp = const_cast<FEMaterialPoint&>(pt);
            double D = 0.0;
            FEReactiveFatigueMaterialPoint* ppr = mp.ExtractData<FEReactiveFatigueMaterialPoint>();
            FEReactivePlasticityMaterialPoint* prp = mp.ExtractData<FEReactivePlasticityMaterialPoint>();
            FEReactivePlasticDamageMaterialPoint* prd = mp.ExtractData<FEReactivePlasticDamageMaterialPoint>();
            FEDamageMaterialPoint* ppd = mp.ExtractData<FEDamageMaterialPoint>();
            FEReactiveViscoelasticMaterialPoint* pve = mp.ExtractData<FEReactiveViscoelasticMaterialPoint>();
            if (prp) D = (float) (1-prp->YieldedBonds());
            else if (prd) D = (float) prd->IntactBonds();
            else if (ppd) D = (float) (1 - ppd->m_D);
            else if (ppr) D = (float) ppr->m_wit;
            else if (pve) {
                const FEReactiveFatigueMaterialPoint* pr = pve->GetPointData(0)->ExtractData< FEReactiveFatigueMaterialPoint>();
                if (pr) D = (float) pr->m_wit;
            }
            return D;
        });
        return true;
    }
    else {
        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);

            float D = 0.0;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
                if (mmp && (m_comp < mmp->Components()))
                {
                    FEReactiveFatigueMaterialPoint* ppr = mmp->GetPointData(m_comp)->ExtractData<FEReactiveFatigueMaterialPoint>();
                    FEReactivePlasticityMaterialPoint* prp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticityMaterialPoint>();
                    FEReactivePlasticDamageMaterialPoint* prd = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticDamageMaterialPoint>();
                    FEDamageMaterialPoint* ppd = mmp->GetPointData(m_comp)->ExtractData<FEDamageMaterialPoint>();
                    FEReactiveViscoelasticMaterialPoint* pve = mmp->GetPointData(m_comp)->ExtractData<FEReactiveViscoelasticMaterialPoint>();
                    if (prp) D += (float) (1-prp->YieldedBonds());
                    else if (prd) D += (float) prd->IntactBonds();
                    else if (ppd) D += (float) (1 - ppd->m_D);
                    else if (ppr) D += (float) ppr->m_wit;
                    else if (pve) {
                        FEReactiveFatigueMaterialPoint* pr = pve->GetPointData(0)->ExtractData<FEReactiveFatigueMaterialPoint>();
                        if (pr) D += (float) pr->m_wit;
                    }
                }
            }
            D /= (float)el.GaussPoints();

            a << D;
        }
        return true;
    }
}

//=============================================================================
FEPlotFatigueBondFraction::FEPlotFatigueBondFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
    m_comp = -1;
}

//-----------------------------------------------------------------------------
bool FEPlotFatigueBondFraction::SetFilter(const char* szfilter)
{
    sscanf(szfilter, "solid[%d]", &m_comp);
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFatigueBondFraction::Save(FEDomain &dom, FEDataStream& a)
{
    if (m_comp == -1) {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& pt) {
			FEMaterialPoint& mp = const_cast<FEMaterialPoint&>(pt);
            float wf = 0.0;
            FEReactiveFatigueMaterialPoint* ppr = mp.ExtractData<FEReactiveFatigueMaterialPoint>();
            FEReactiveViscoelasticMaterialPoint* pve = mp.ExtractData<FEReactiveViscoelasticMaterialPoint>();
            if (ppr) wf = (float) ppr->m_wft;
            else if (pve) {
                FEReactiveFatigueMaterialPoint* pr = pve->GetPointData(0)->ExtractData<FEReactiveFatigueMaterialPoint>();
                if (pr) wf = (float) pr->m_wft;
            }
            return wf;
        });
        return true;
    }
    else {
        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);

            float wf = 0.0;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
                if (mmp && (m_comp < mmp->Components()))
                {
                    FEReactiveFatigueMaterialPoint* ppr = mmp->GetPointData(m_comp)->ExtractData<FEReactiveFatigueMaterialPoint>();
                    FEReactiveViscoelasticMaterialPoint* pve = mmp->GetPointData(m_comp)->ExtractData<FEReactiveViscoelasticMaterialPoint>();
                    if (ppr) wf += (float) ppr->m_wft;
                    else if (pve) {
                        FEReactiveFatigueMaterialPoint* pr = pve->GetPointData(0)->ExtractData<FEReactiveFatigueMaterialPoint>();
                        if (pr) wf += (float) pr->m_wft;
                    }
                }
            }
            wf /= (float)el.GaussPoints();

            a << wf;
        }
        return true;
    }
}

//=============================================================================
FEPlotYieldedBondFraction::FEPlotYieldedBondFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
    m_comp = -1;
}

//-----------------------------------------------------------------------------
bool FEPlotYieldedBondFraction::SetFilter(const char* szfilter)
{
    sscanf(szfilter, "solid[%d]", &m_comp);
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotYieldedBondFraction::Save(FEDomain &dom, FEDataStream& a)
{
    if (m_comp == -1) {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
            float wy = 0.0;
            const FEReactivePlasticityMaterialPoint* prp = mp.ExtractData<FEReactivePlasticityMaterialPoint>();
            const FEReactivePlasticDamageMaterialPoint* ppp = mp.ExtractData<FEReactivePlasticDamageMaterialPoint>();
            if (prp) wy = (float) prp->YieldedBonds();
            else if (ppp) wy = (float) ppp->YieldedBonds();
            return wy;
        });
        return true;
    }
    else {
        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);

            float wy = 0.0;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
                if (mmp && (m_comp < mmp->Components()))
                {
                    const FEReactivePlasticityMaterialPoint* prp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticityMaterialPoint>();
                    const FEReactivePlasticDamageMaterialPoint* ppp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticDamageMaterialPoint>();
                    if (prp) wy += (float) prp->YieldedBonds();
                    else if (ppp) wy += (float) ppp->YieldedBonds();
                }
            }
            wy /= (float)el.GaussPoints();

            a << wy;
        }
        return true;
    }
}

//=============================================================================
FEPlotReactivePlasticityHeatSupply::FEPlotReactivePlasticityHeatSupply(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
    m_comp = -1;
}

//-----------------------------------------------------------------------------
bool FEPlotReactivePlasticityHeatSupply::SetFilter(const char* szfilter)
{
    sscanf(szfilter, "solid[%d]", &m_comp);
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotReactivePlasticityHeatSupply::Save(FEDomain &dom, FEDataStream& a)
{
    if (m_comp == -1) {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
            float Rhat = 0.0;
            const FEReactivePlasticityMaterialPoint* prp = mp.ExtractData<FEReactivePlasticityMaterialPoint>();
            const FEReactivePlasticDamageMaterialPoint* ppp = mp.ExtractData<FEReactivePlasticDamageMaterialPoint>();
            if (prp) Rhat = (float) prp->m_Rhat;
            else if (ppp) Rhat = (float) ppp->m_Rhat;
            return Rhat;
        });
        return true;
    }
    else {
        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);

            float Rhat = 0.0;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
                if (mmp && (m_comp < mmp->Components()))
                {
                    const FEReactivePlasticityMaterialPoint* prp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticityMaterialPoint>();
                    const FEReactivePlasticDamageMaterialPoint* ppp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticDamageMaterialPoint>();
                    if (prp) Rhat += (float) prp->m_Rhat;
                    else if (ppp) Rhat += (float) ppp->m_Rhat;
                }
            }
            Rhat /= (float)el.GaussPoints();

            a << Rhat;
        }
        return true;
    }
}


//=============================================================================
FEPlotOctahedralPlasticStrain::FEPlotOctahedralPlasticStrain(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
    m_comp = -1;
}

//-----------------------------------------------------------------------------
bool FEPlotOctahedralPlasticStrain::SetFilter(const char* szfilter)
{
    sscanf(szfilter, "solid[%d]", &m_comp);
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotOctahedralPlasticStrain::Save(FEDomain &dom, FEDataStream& a)
{
    if (m_comp == -1) {
        writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
            float gp = 0.0;
            const FEReactivePlasticityMaterialPoint* prp = mp.ExtractData<FEReactivePlasticityMaterialPoint>();
            const FEReactivePlasticDamageMaterialPoint* ppp = mp.ExtractData<FEReactivePlasticDamageMaterialPoint>();
            if (prp && prp->m_gp.size()) gp = (float) prp->m_gp[0];
            else if (ppp && ppp->m_gp.size()) gp = (float) ppp->m_gp[0];
            return gp;
        });
        return true;
    }
    else {
        for (int i = 0; i < dom.Elements(); ++i)
        {
            FEElement& el = dom.ElementRef(i);

            float gp = 0.0;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
                if (mmp && (m_comp < mmp->Components()))
                {
                    const FEReactivePlasticityMaterialPoint* prp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticityMaterialPoint>();
                    const FEReactivePlasticDamageMaterialPoint* ppp = mmp->GetPointData(m_comp)->ExtractData<FEReactivePlasticDamageMaterialPoint>();
                    if (prp && prp->m_gp.size()) gp += (float) prp->m_gp[0];
                    else if (ppp && ppp->m_gp.size()) gp += (float) ppp->m_gp[0];
                }
            }
            gp /= (float)el.GaussPoints();

            a << gp;
        }
        return true;
    }
}

//-----------------------------------------------------------------------------
bool FEPlotMixtureVolumeFraction::Save(FEDomain &dom, FEDataStream &a)
{
	// extract the mixture material
	FEMaterial* pmat = dom.GetMaterial();
	FEElasticMixture* pm = dynamic_cast<FEElasticMixture*>(pmat);
	if (pm == 0) return false;

	writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
		const FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();
		return pt.m_w[0];
	});

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotUT4NodalStresses::Save(FEDomain& dom, FEDataStream& a)
{
	// make sure this is a UT4 domain
	FEUT4Domain* pd = dynamic_cast<FEUT4Domain*>(&dom);
	if (pd == 0) return false;

	// write the nodal values
	writeNodalValues<mat3ds>(dom, a, [=](int i) {
		FEUT4Domain::UT4NODE& n = pd->UT4Node(i);
		return n.si;
	});

	return true;
}

//==============================================================================
//                  R I G I D   B O D Y   D A T A
//==============================================================================

//-----------------------------------------------------------------------------
bool FEPlotRigidDisplacement::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;
    
	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());

	// store the rigid body position
	// TODO: why do we not store displacement?
	a << rb.m_rt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidVelocity::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    
	// store the rigid velocity
	a << rb.m_vt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAcceleration::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    
	// store rigid body acceleration
	a << rb.m_at;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidRotation::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
	vec3d q = rb.GetRotation().GetRotationVector();
    
	// store rotation vector
	a << q;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularVelocity::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    
	// store rigid angular velocity
	a << rb.m_wt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularAcceleration::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    
	// store angular acceleration
	a << rb.m_alt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidKineticEnergy::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    vec3d v = rb.m_vt;
    double m = rb.m_mass;
    vec3d w = rb.m_wt;
	mat3d Rt = rb.GetRotation().RotationMatrix();
    mat3ds Jt = (Rt*rb.m_moi*Rt.transpose()).sym();
    double ke = ((v*v)*m + w*(Jt*w))/2;
    
	// store kinetic energy
	a << ke;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidLinearMomentum::Save(FEDomain& dom, FEDataStream& a)
{
    // get the rigid material
    FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

    // get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    
    // store linear momentum (mass x velocity)
    a << rb.m_vt*rb.m_mass;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularMomentum::Save(FEDomain& dom, FEDataStream& a)
{
    // get the rigid material
    FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

    // get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());
    
    // store angular momentum (mass moment of inertia x angular velocity)
	mat3d Rt = rb.GetRotation().RotationMatrix();
    mat3ds Jt = (Rt*rb.m_moi*Rt.transpose()).sym();
    
    a << Jt*rb.m_wt;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidEuler::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());

	// get the Euler angles
	double E[3];
	quat2euler(rb.GetRotation(), E);
    
	// store Euler
	a << E[0] << E[1] << E[2];
    
	return true;
}

//-----------------------------------------------------------------------------
// TODO: I think this already gets stored somewhere
bool FEPlotRigidRotationVector::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pm);
	if (prm == 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidBody& rb = *fem.GetRigidBody(prm->GetRigidBodyID());

	// get the rotation vector and angle
	vec3d r = rb.GetRotation().GetRotationVector();
    
	// store rotation vector
	a << r;
    
	return true;
}

//=============================================================================
bool FEPlotRigidReactionForce::Save(FEDomain& dom, FEDataStream& a)
{
	// get the material
	FEMaterial* pmat = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pmat);
	if (prm == 0) return false;

	// get the rigid body ID
	int nrid = prm->GetRigidBodyID();
	if (nrid < 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& rb = *fem.GetRigidBody(nrid);

	a << rb.m_Fr;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidReactionTorque::Save(FEDomain& dom, FEDataStream& a)
{
	// get the material
	FEMaterial* pmat = dom.GetMaterial();
	FERigidMaterial* prm = dynamic_cast<FERigidMaterial*>(pmat);
	if (prm == 0) return false;

	// get the rigid body ID
	int nrid = prm->GetRigidBodyID();
	if (nrid < 0) return false;

	// get the rigid body
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& rb = *fem.GetRigidBody(nrid);

	a << rb.m_Mr;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotStressError::Save(FEDomain& dom, FEDataStream& a)
{
	writeRelativeError(dom, a, [](FEMaterialPoint& mp) {
		FEElasticMaterialPoint* ep = mp.ExtractData<FEElasticMaterialPoint>();
		mat3ds& s = ep->m_s;
		double v = s.effective_norm();
		return v;
	});
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFiberTargetStretch::Save(FEDomain& dom, FEDataStream& a)
{
	FEMaterial* mat = dom.GetMaterial();
	FEPrestrainMaterial * pmat = dynamic_cast<FEPrestrainMaterial*>(mat);
	if (pmat == nullptr) return false;

	// get the elastic component
	FEProperty* prop = mat->FindProperty("elastic");
	if (prop == nullptr) return false;

	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(prop->get(0));
	if (pme == 0) return false;

	// get the fiber property
	FEVec3dValuator* vec = dynamic_cast<FEVec3dValuator*>(pme->GetProperty("fiber"));
	if (vec == 0) return false;

	// we're good so store the in-situ stretch
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int nint = e.GaussPoints();
		double lam = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(j)->GetPointData(0);
			FEPrestrainMaterialPoint& pp = *mp.ExtractData<FEPrestrainMaterialPoint>();

			mat3d Fp = pp.initialPrestrain();
			mat3d Q = mat->GetLocalCS(mp);
			vec3d a0 = vec->unitVector(mp);
			vec3d ar = Q * a0;
			vec3d a = Fp*ar;
			double lamp = a.norm();

			lam += lamp;
		}
		lam /= (double)nint;

		a << lam;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPreStrainStretch::Save(FEDomain& dom, FEDataStream& a)
{
	FEMaterial* mat = dom.GetMaterial();
	FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(mat);
	if (pmat == 0) return false;

	// get the elastic component
	FEProperty* prop = mat->FindProperty("elastic");
	if (prop == nullptr) return false;

	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(prop->get(0));
	if (pme== 0) return false;

	// get the fiber property
	FEVec3dValuator* vec = dynamic_cast<FEVec3dValuator*>(pme->GetProperty("fiber"));
	if (vec == 0) return false;

	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int nint = e.GaussPoints();
		double lam = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(j)->GetPointData(0);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEPrestrainMaterialPoint& pp = *mp.ExtractData<FEPrestrainMaterialPoint>();

			mat3d& F = pt.m_F;
			mat3d Fp = pp.prestrain();
			mat3d Ft = F*Fp;

			mat3d Q = mat->GetLocalCS(mp);
			vec3d a0 = vec->unitVector(mp);
			vec3d ar = Q * a0;
			vec3d a = Ft*ar;

			double lambda = a.norm();
			lam += lambda;
		}
		lam /= (double)nint;

		a << lam;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPreStrainStretchError::Save(FEDomain& dom, FEDataStream& a)
{
	FEMaterial* mat = dom.GetMaterial();
	FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(mat);
	if (pmat == 0) return false;

	// get the elastic component
	FEProperty* prop = mat->FindProperty("elastic");
	if (prop == nullptr) return false;

	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(prop->get(0));
	if (pme == 0) return false;

	// get the fiber property
	FEVec3dValuator* vec = dynamic_cast<FEVec3dValuator*>(pme->GetProperty("fiber"));
	if (vec == 0) return false;

	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int nint = e.GaussPoints();
		double err = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(j)->GetPointData(0);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEPrestrainMaterialPoint& pp = *mp.ExtractData<FEPrestrainMaterialPoint>();

			// initial fiber vector
			mat3d Q = mat->GetLocalCS(mp);
			vec3d a0 = vec->unitVector(mp);
			vec3d ar = Q * a0;

			// target stretch
			mat3d Fp = pp.initialPrestrain();
			vec3d a = Fp*ar;
			double lam_trg = a.norm();

			// current stretch
			mat3d& F = pt.m_F;
			Fp = pp.prestrain();
			mat3d Ft = F*Fp;
			a = Ft*ar;

			double lam_cur = a.norm();

			err += fabs(lam_cur / lam_trg - 1.0);
		}
		err /= (double)nint;

		a << err;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPreStrainCorrection::Save(FEDomain& dom, FEDataStream& a)
{
	FEMaterial* mat = dom.GetMaterial();
	FEElasticMaterial* solidMat = mat->ExtractProperty<FEElasticMaterial>();
	FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(solidMat);
	if (pmat == 0) return false;

	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int nint = e.GaussPoints();
		mat3d Fc; Fc.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(j)->GetPointData(0);
			FEPrestrainMaterialPoint& pt = *mp.ExtractData<FEPrestrainMaterialPoint>();

			Fc += pt.PrestrainCorrection();
		}
		Fc *= 1.0 / (double)nint;

		a << Fc;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRPreStrainCorrection::Save(FEDomain& dom, FEDataStream& a)
{
	const int LUT[9][2] = { { 0,0 },{ 0,1 },{ 0,2 },{ 1,0 },{ 1,1 },{ 1,2 },{ 2,0 },{ 2,1 },{ 2,2 } };

	FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(dom.GetMaterial());
	if (pmat == 0) return false;

	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	int NN = sd.Nodes();
	int NE = sd.Elements();

	// build the element data array
	vector< vector<double> > ED;
	ED.resize(NE);
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		ED[i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	vector<double> val[9];

	// loop over stress components
	for (int n = 0; n<9; ++n)
	{
		// fill the ED array
		for (int i = 0; i<NE; ++i)
		{
			FESolidElement& el = sd.Element(i);
			int nint = el.GaussPoints();
			for (int j = 0; j<nint; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j)->GetPointData(0);
				FEPrestrainMaterialPoint& pt = *mp.ExtractData<FEPrestrainMaterialPoint>();
				const mat3d& F = pt.PrestrainCorrection();
				ED[i][j] = F(LUT[n][0], LUT[n][1]);
			}
		}

		// project to nodes
		map.Project(sd, ED, val[n]);
	}

	// copy results to archive
	for (int i = 0; i<NN; ++i)
	{
		a << val[0][i];
		a << val[1][i];
		a << val[2][i];
		a << val[3][i];
		a << val[4][i];
		a << val[5][i];
		a << val[6][i];
		a << val[7][i];
		a << val[8][i];
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPreStrainCompatibility::Save(FEDomain& dom, FEDataStream& a)
{
	const int LUT[9][2] = { { 0,0 },{ 0,1 },{ 0,2 },{ 1,0 },{ 1,1 },{ 1,2 },{ 2,0 },{ 2,1 },{ 2,2 } };

	// make sure this is a pre-strain material
	FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(dom.GetMaterial());
	if (pmat == 0) return false;

	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	int NE = sd.Elements();

	// STEP 1 - first we do an SPR recovery of the pre-strain gradient

	// build the element data array
	vector< vector<double> > ED;
	ED.resize(NE);
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		ED[i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	vector<double> val[9];

	// create a global-to-local node list
	FEMesh& mesh = *dom.GetMesh();
	vector<int> g2l; g2l.assign(mesh.Nodes(), -1);
	int nn = 0;
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd.Element(i);
		int neln = el.Nodes();
		for (int j = 0; j<neln; ++j)
		{
			if (g2l[el.m_node[j]] == -1) g2l[el.m_node[j]] = nn++;
		}
	}

	// loop over tensor components
	for (int n = 0; n<9; ++n)
	{
		// fill the ED array
		for (int i = 0; i<NE; ++i)
		{
			FESolidElement& el = sd.Element(i);
			int nint = el.GaussPoints();
			for (int j = 0; j<nint; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j)->GetPointData(0);
				FEPrestrainMaterialPoint& pt = *mp.ExtractData<FEPrestrainMaterialPoint>();
				mat3d Fp = pt.prestrain();
				ED[i][j] = Fp(LUT[n][0], LUT[n][1]);
			}
		}

		// project to nodes
		map.Project(sd, ED, val[n]);
	}

	// STEP 2 - now we calculate the gradient of the nodal values at the integration points
	vector<double> vn(FEElement::MAX_NODES);
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd.Element(i);
		int neln = el.Nodes();
		int nint = el.GaussPoints();

		vector<mat3d> gF[3];
		gF[0].assign(nint, mat3dd(0));
		gF[1].assign(nint, mat3dd(0));
		gF[2].assign(nint, mat3dd(0));
		for (int n = 0; n<9; ++n)
		{
			// get the nodal values
			for (int m = 0; m<neln; ++m) vn[m] = val[n][g2l[el.m_node[m]]];

			// calculate the gradient at the integration points
			for (int j = 0; j<nint; ++j)
			{
				vec3d g = sd.gradient(el, vn, j);

				gF[0][j](LUT[n][0], LUT[n][1]) = g.x;
				gF[1][j](LUT[n][0], LUT[n][1]) = g.y;
				gF[2][j](LUT[n][0], LUT[n][1]) = g.z;
			}
		}

		double c = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			mat3d C;
			C(0, 0) = gF[1][j](0, 2) - gF[2][j](0, 1);
			C(0, 1) = gF[1][j](1, 2) - gF[2][j](1, 1);
			C(0, 2) = gF[1][j](2, 2) - gF[2][j](2, 1);

			C(1, 0) = gF[2][j](0, 0) - gF[0][j](0, 2);
			C(1, 1) = gF[2][j](1, 0) - gF[0][j](1, 2);
			C(1, 2) = gF[2][j](2, 0) - gF[0][j](2, 2);

			C(2, 0) = gF[0][j](0, 1) - gF[1][j](0, 0);
			C(2, 1) = gF[0][j](1, 1) - gF[1][j](1, 0);
			C(2, 2) = gF[0][j](2, 1) - gF[1][j](2, 0);

			c += sqrt(C.dotdot(C));
		}
		c /= (double)nint;

		// store the compatibility
		a << c;
	}

	return true;
}

bool FEPlotDiscreteElementStretch::Save(FEDomain& dom, FEDataStream& a)
{
	FEDiscreteDomain* pdiscreteDomain = dynamic_cast<FEDiscreteDomain*>(&dom);
	if (pdiscreteDomain == nullptr) return false;
	FEDiscreteDomain& discreteDomain = *pdiscreteDomain;

	FEMesh& mesh = *dom.GetMesh();
	int NE = discreteDomain.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEDiscreteElement& el = discreteDomain.Element(i);
		
		vec3d ra0 = mesh.Node(el.m_node[0]).m_r0;
		vec3d ra1 = mesh.Node(el.m_node[0]).m_rt;
		vec3d rb0 = mesh.Node(el.m_node[1]).m_r0;
		vec3d rb1 = mesh.Node(el.m_node[1]).m_rt;

		double L0 = (rb0 - ra0).norm();
		double Lt = (rb1 - ra1).norm();

		double l = Lt / L0;
		a << l;
	}

	return true;
}

bool FEPlotDiscreteElementElongation::Save(FEDomain& dom, FEDataStream& a)
{
	FEDiscreteDomain* pdiscreteDomain = dynamic_cast<FEDiscreteDomain*>(&dom);
	if (pdiscreteDomain == nullptr) return false;
	FEDiscreteDomain& discreteDomain = *pdiscreteDomain;

	FEMesh& mesh = *dom.GetMesh();
	int NE = discreteDomain.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEDiscreteElement& el = discreteDomain.Element(i);
		
		vec3d ra0 = mesh.Node(el.m_node[0]).m_r0;
		vec3d ra1 = mesh.Node(el.m_node[0]).m_rt;
		vec3d rb0 = mesh.Node(el.m_node[1]).m_r0;
		vec3d rb1 = mesh.Node(el.m_node[1]).m_rt;

		double L0 = (rb0 - ra0).norm();
		double Lt = (rb1 - ra1).norm();

		double l = Lt - L0;
		a << l;
	}

	return true;
}

bool FEPlotDiscreteElementPercentElongation::Save(FEDomain& dom, FEDataStream& a)
{
	FEDiscreteDomain* pdiscreteDomain = dynamic_cast<FEDiscreteDomain*>(&dom);
	if (pdiscreteDomain == nullptr) return false;
	FEDiscreteDomain& discreteDomain = *pdiscreteDomain;

	FEMesh& mesh = *dom.GetMesh();
	int NE = discreteDomain.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEDiscreteElement& el = discreteDomain.Element(i);
		
		vec3d ra0 = mesh.Node(el.m_node[0]).m_r0;
		vec3d ra1 = mesh.Node(el.m_node[0]).m_rt;
		vec3d rb0 = mesh.Node(el.m_node[1]).m_r0;
		vec3d rb1 = mesh.Node(el.m_node[1]).m_rt;

		double L0 = (rb0 - ra0).norm();
		double Lt = (rb1 - ra1).norm();

		double l = (Lt - L0)/L0;
		a << l;
	}

	return true;
}

bool FEPlotDiscreteElementForce::Save(FEDomain& dom, FEDataStream& a)
{
	FEDiscreteElasticDomain* pdiscreteDomain = dynamic_cast<FEDiscreteElasticDomain*>(&dom);
	if (pdiscreteDomain == nullptr) return false;
	FEDiscreteElasticDomain& discreteDomain = *pdiscreteDomain;

	int NE = discreteDomain.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEDiscreteElement& el = discreteDomain.Element(i);

		// get the (one) material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(0);
		FEDiscreteElasticMaterialPoint& dmp = *mp.ExtractData<FEDiscreteElasticMaterialPoint>();

		// write the force
		a << dmp.m_Ft;
	}

	return true;
}

bool FEPlotDiscreteElementSignedForce::Save(FEDomain& dom, FEDataStream& a)
{
	FEDiscreteElasticDomain* pdiscreteDomain = dynamic_cast<FEDiscreteElasticDomain*>(&dom);
	if (pdiscreteDomain == nullptr) return false;
	FEDiscreteElasticDomain& discreteDomain = *pdiscreteDomain;

	int NE = discreteDomain.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEDiscreteElement& el = discreteDomain.Element(i);

		// get the (one) material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(0);
		FEDiscreteElasticMaterialPoint& dmp = *mp.ExtractData<FEDiscreteElasticMaterialPoint>();

		vec3d ra1 = dom.Node(el.m_lnode[0]).m_rt;
		vec3d rb1 = dom.Node(el.m_lnode[1]).m_rt;
		vec3d e = rb1 - ra1; e.unit();

		vec3d F = dmp.m_Ft;

		double Fm = F * e;

		// write the force
		a << Fm;
	}

	return true;
}

bool FEPlotDiscreteElementStrainEnergy::Save(FEDomain& dom, FEDataStream& a)
{
	FEDiscreteElasticDomain* pdiscreteDomain = dynamic_cast<FEDiscreteElasticDomain*>(&dom);
	if (pdiscreteDomain == nullptr) return false;
	FEDiscreteElasticDomain& discreteDomain = *pdiscreteDomain;

	FEDiscreteElasticMaterial* discreteMaterial = dynamic_cast<FEDiscreteElasticMaterial*>(discreteDomain.GetMaterial());

	int NE = discreteDomain.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEDiscreteElement& el = discreteDomain.Element(i);

		// get the (one) material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(0);
		FEDiscreteElasticMaterialPoint& dmp = *mp.ExtractData<FEDiscreteElasticMaterialPoint>();

		// write the strain energy
		a << discreteMaterial->StrainEnergy(dmp);
	}

	return true;

}

//=================================================================================================
FEPlotContinuousDamage_::FEPlotContinuousDamage_(FEModel* fem, int n) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM)
{
	m_propIndex = 0;
	m_comp = n;
}

bool FEPlotContinuousDamage_::SetFilter(const char* sz)
{
	m_prop = sz;
	return true;
}

bool FEPlotContinuousDamage_::Save(FEDomain& dom, FEDataStream& a)
{
	// get the material
	FEMaterial* domMat = dom.GetMaterial();
	if (domMat == nullptr) return false;

	// get the fiber damage component
	FEDamageElasticFiber* mat = nullptr;
	if (m_prop.empty()) mat = dynamic_cast<FEDamageElasticFiber*>(domMat);
	else
	{
		ParamString ps(m_prop.c_str());
		m_propIndex = ps.Index();
		mat = dynamic_cast<FEDamageElasticFiber*>(domMat->GetProperty(ps));
	}
	if (mat == nullptr) return false;

	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		double D = 0.0;
		int nint = el.GaussPoints();
		for (int j = 0; j < nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEMaterialPoint* pt = mp.GetPointData(m_propIndex);
			double Dj = mat->Damage(*pt, m_comp);

			D += Dj;
		}
		D /= (double)nint;

		a << D;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRVEgenerations::Save(FEDomain& dom, FEDataStream& a)
{
    int N = dom.Elements();
    FEElasticMaterial* pmat = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pmat == nullptr) return false;
    FEReactiveViscoelasticMaterial* rvmat = dynamic_cast<FEReactiveViscoelasticMaterial*>(pmat);
    FEUncoupledReactiveViscoelasticMaterial* rumat = dynamic_cast<FEUncoupledReactiveViscoelasticMaterial*>(pmat);
    if (rvmat) {
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int nint = el.GaussPoints();
            double bmf = 0;
            for (int j=0; j<nint; ++j)
                bmf += rvmat->RVEGenerations(*el.GetMaterialPoint(j));
            a << bmf/nint;
        }
    }
    else if (rumat) {
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int nint = el.GaussPoints();
            double bmf = 0;
            for (int j=0; j<nint; ++j)
                bmf += rumat->RVEGenerations(*el.GetMaterialPoint(j));
            a << bmf/nint;
        }
    }
    else {
        int NC = pmat->Properties();
        // check all elements
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int n = 0;
            double bmf = 0;
            // check all properties
            for (int ic=0; ic < NC; ++ic) {
                FEReactiveViscoelasticMaterial* rvmat = pmat->GetProperty(ic)->ExtractProperty<FEReactiveViscoelasticMaterial>();
                FEUncoupledReactiveViscoelasticMaterial* rumat = dynamic_cast<FEUncoupledReactiveViscoelasticMaterial*>(pmat);
                if (rvmat) {
                    int nint = el.GaussPoints();
                    for (int j=0; j<nint; ++j)
                    {
						FEMaterialPoint& mp = *el.GetMaterialPoint(j)->GetPointData(ic);
                        bmf += rvmat->RVEGenerations(mp);
                        ++n;
                    }
                }
                else if (rumat) {
                    int nint = el.GaussPoints();
                    for (int j=0; j<nint; ++j)
                    {
						FEMaterialPoint& mp = *el.GetMaterialPoint(j)->GetPointData(ic);
						bmf += rumat->RVEGenerations(mp);
                        ++n;
                    }
                }
            }
            if (n > 0) bmf /= n;
            a << bmf;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRVEbonds::Save(FEDomain& dom, FEDataStream& a)
{
    int N = dom.Elements();
    FEElasticMaterial* pmat = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pmat == nullptr) return false;
    FEReactiveViscoelasticMaterial* rvmat = dynamic_cast<FEReactiveViscoelasticMaterial*>(pmat);
    FEUncoupledReactiveViscoelasticMaterial* rumat = dynamic_cast<FEUncoupledReactiveViscoelasticMaterial*>(pmat);
    if (rvmat) {
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int nint = el.GaussPoints();
            double bmf = 0;
            for (int j=0; j<nint; ++j)
                bmf += rvmat->ReformingBondMassFraction(*rvmat->GetBondMaterialPoint(*el.GetMaterialPoint(j)));
            a << bmf/nint;
        }
    }
    else if (rumat) {
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int nint = el.GaussPoints();
            double bmf = 0;
            for (int j=0; j<nint; ++j)
                bmf += rumat->ReformingBondMassFraction(*rumat->GetBondMaterialPoint(*el.GetMaterialPoint(j)));
            a << bmf/nint;
        }
    }
    else {
        int NC = pmat->Properties();
        // check all elements
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int n = 0;
            double bmf = 0;
            // check all properties
            for (int ic=0; ic < NC; ++ic) {
                FEReactiveViscoelasticMaterial* rvmat = pmat->GetProperty(ic)->ExtractProperty<FEReactiveViscoelasticMaterial>();
                FEUncoupledReactiveViscoelasticMaterial* rumat = dynamic_cast<FEUncoupledReactiveViscoelasticMaterial*>(pmat);
                if (rvmat) {
                    int nint = el.GaussPoints();
                    for (int j=0; j<nint; ++j)
                    {
                        bmf += rvmat->ReformingBondMassFraction(*rvmat->GetBondMaterialPoint(*el.GetMaterialPoint(j)->GetPointData(ic)));
                        ++n;
                    }
                }
                else if (rumat) {
                    int nint = el.GaussPoints();
                    for (int j=0; j<nint; ++j)
                    {
                        bmf += rumat->ReformingBondMassFraction(*rumat->GetBondMaterialPoint(*el.GetMaterialPoint(j)->GetPointData(ic)));
                        ++n;
                    }
                }
            }
            if (n > 0) bmf /= n;
            a << bmf;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRVEstrain::Save(FEDomain& dom, FEDataStream& a)
{
    int N = dom.Elements();
    FEElasticMaterial* pmat = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pmat == nullptr) return false;
    FEReactiveViscoelasticMaterial* rvmat = dynamic_cast<FEReactiveViscoelasticMaterial*>(pmat);
    FEUncoupledReactiveViscoelasticMaterial* rumat = dynamic_cast<FEUncoupledReactiveViscoelasticMaterial*>(pmat);
    if (rvmat) {
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int nint = el.GaussPoints();
            double bmf = 0;
            for (int j=0; j<nint; ++j)
                bmf += rvmat->ScalarStrain(*rvmat->GetBondMaterialPoint(*el.GetMaterialPoint(j)));
            a << bmf/nint;
        }
    }
    else if (rumat) {
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int nint = el.GaussPoints();
            double bmf = 0;
            for (int j=0; j<nint; ++j)
                bmf += rumat->ScalarStrain(*rumat->GetBondMaterialPoint(*el.GetMaterialPoint(j)));
            a << bmf/nint;
        }
    }
    else {
        int NC = pmat->Properties();
        // check all elements
        for (int iel=0; iel<N; ++iel)
        {
            FEElement& el = dom.ElementRef(iel);
            
            int n = 0;
            double bmf = 0;
            // check all properties
            for (int ic=0; ic < NC; ++ic) {
                FEReactiveViscoelasticMaterial* rvmat = pmat->GetProperty(ic)->ExtractProperty<FEReactiveViscoelasticMaterial>();
                FEUncoupledReactiveViscoelasticMaterial* rumat = dynamic_cast<FEUncoupledReactiveViscoelasticMaterial*>(pmat);
                if (rvmat) {
                    int nint = el.GaussPoints();
                    for (int j=0; j<nint; ++j)
                    {
                        bmf += rvmat->ScalarStrain(*rvmat->GetBondMaterialPoint(*el.GetMaterialPoint(j)->GetPointData(ic)));
                        ++n;
                    }
                }
                else if (rumat) {
                    int nint = el.GaussPoints();
                    for (int j=0; j<nint; ++j)
                    {
                        bmf += rumat->ScalarStrain(*rumat->GetBondMaterialPoint(*el.GetMaterialPoint(j)->GetPointData(ic)));
                        ++n;
                    }
                }
            }
            if (n > 0) bmf /= n;
            a << bmf;
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
class FEStrongBondSED
{
public:
    FEStrongBondSED(FEElasticMaterial* pm) : m_mat(pm) {}
    double operator()(const FEMaterialPoint& mp)
    {
        return m_mat->StrongBondSED(const_cast<FEMaterialPoint&>(mp));
    }
private:
    FEElasticMaterial*    m_mat;
};

bool FEPlotStrongBondSED::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FEStrongBondSED W(pme);
        writeAverageElementValue<double>(dom, a, W);
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
class FEWeakBondSED
{
public:
    FEWeakBondSED(FEElasticMaterial* pm) : m_mat(pm) {}
    double operator()(const FEMaterialPoint& mp)
    {
        return m_mat->WeakBondSED(const_cast<FEMaterialPoint&>(mp));
    }
private:
    FEElasticMaterial*    m_mat;
};

bool FEPlotWeakBondSED::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FEWeakBondSED W(pme);
        writeAverageElementValue<double>(dom, a, W);
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
class FEStrongBondDevSED
{
public:
    FEStrongBondDevSED(FEUncoupledMaterial* pm) : m_mat(pm) {}
    double operator()(const FEMaterialPoint& mp)
    {
        return m_mat->StrongBondDevSED(const_cast<FEMaterialPoint&>(mp));
    }
private:
    FEUncoupledMaterial*    m_mat;
};

bool FEPlotStrongBondDevSED::Save(FEDomain& dom, FEDataStream& a)
{
    FEUncoupledMaterial* pme = dom.GetMaterial()->ExtractProperty<FEUncoupledMaterial>();
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FEStrongBondDevSED W(pme);
        writeAverageElementValue<double>(dom, a, W);
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
class FEWeakBondDevSED
{
public:
    FEWeakBondDevSED(FEUncoupledMaterial* pm) : m_mat(pm) {}
    double operator()(const FEMaterialPoint& mp)
    {
        return m_mat->WeakBondDevSED(const_cast<FEMaterialPoint&>(mp));
    }
private:
    FEUncoupledMaterial*    m_mat;
};

bool FEPlotWeakBondDevSED::Save(FEDomain& dom, FEDataStream& a)
{
    FEUncoupledMaterial* pme = dom.GetMaterial()->ExtractProperty<FEUncoupledMaterial>();
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FEWeakBondDevSED W(pme);
        writeAverageElementValue<double>(dom, a, W);
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotTrussStretch::Save(FEDomain& dom, FEDataStream& a)
{
	FETrussDomain* td = dynamic_cast<FETrussDomain*>(&dom);
	if (td == nullptr) return false;

	for (int i = 0; i < td->Elements(); ++i)
	{
		FETrussElement& el = td->Element(i);
		a << el.m_lam;
	}
	return true;
}
