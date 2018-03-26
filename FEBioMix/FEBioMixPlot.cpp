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
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FEBiphasicContactSurface.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include <FECore/FEModel.h>

//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotMixtureFluidFlowRate::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    // Evaluate this field only for a specific domain, by checking domain name
    if (pcs->GetName() != m_szdom) return false;
    
    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    FEMesh* m_pMesh = pcs->GetMesh();
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        m_elem.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
            m_elem[j] = m_pMesh->FindElementFromID(pcs->FindElement(el));
        }
        m_binit = false;
    }
    
    // calculate net flow rate normal to this surface
    for (int j=0; j<NF; ++j)
    {
        // get the element this surface element belongs to
        FEElement* pe = m_elem[j];
        if (pe)
        {
            // evaluate the average fluid flux in this element
            int nint = pe->GaussPoints();
            vec3d w(0,0,0);
            for (int n=0; n<nint; ++n)
            {
                FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                FEBiphasicMaterialPoint* ptb = mp.ExtractData<FEBiphasicMaterialPoint>();
                if (ptb) w += ptb->m_w;
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

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

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
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average pressure
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				
				if (pt) ew += pt->m_pa;
			}
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if ((dynamic_cast<FEBiphasicShellDomain*>(&bsd)) ||
             (dynamic_cast<FEBiphasicSoluteShellDomain*>(&bsd)) ||
             (dynamic_cast<FEMultiphasicShellDomain*>(&bsd))
             )
    {
        for (int i=0; i<bsd.Elements(); ++i)
        {
            FEShellElement& el = bsd.Element(i);
            
            // calculate average pressure
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
                
                if (pt) ew += pt->m_pa;
            }
            ew /= el.GaussPoints();
            
            a << ew;
        }
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
		for (int i=0; i<bd->Elements(); ++i)
		{
			FESolidElement& el = bd->Element(i);

			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());

				if (pt) ew += pt->m_w;
			}
			ew /= el.GaussPoints();

			a << ew;
		}
		return true;
	}
    else if (bsd && (
                     (dynamic_cast<FEBiphasicShellDomain*>(bsd)) ||
                     (dynamic_cast<FEBiphasicSoluteShellDomain*>(bsd)) ||
                     (dynamic_cast<FEMultiphasicShellDomain*>(bsd))
                     ))
    {
        for (int i=0; i<bsd->Elements(); ++i)
        {
            FEShellElement& el = bsd->Element(i);
            
            // calculate average flux
            vec3d ew = vec3d(0,0,0);
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
                
                if (pt) ew += pt->m_w;
            }
            ew /= el.GaussPoints();
            
            a << ew;
        }
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
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);

			int nint = el.GaussPoints();
			int neln = el.Nodes();

			// fluid flux at gauss points
			double vi[3][FEElement::MAX_NODES];
			for (int j=0; j<nint; ++j)
			{
				FEBiphasicMaterialPoint* pt = el.GetMaterialPoint(j)->ExtractData<FEBiphasicMaterialPoint>(); assert(pt);
				vi[0][j] = pt->m_w.x;
				vi[1][j] = pt->m_w.y;
				vi[2][j] = pt->m_w.z;
			}

			// project to nodes
			double vn[3][FEElement::MAX_NODES];
			el.project_to_nodes(vi[0], vn[0]);
			el.project_to_nodes(vi[1], vn[1]);
			el.project_to_nodes(vi[2], vn[2]);

			// output data
			for (int j=0; j<neln; ++j)
			{
				a.push_back((float)vn[0][j]);
				a.push_back((float)vn[1][j]);
				a.push_back((float)vn[2][j]);
			}
		}
		return true;
	}
	return false;
}

/*
//-----------------------------------------------------------------------------
bool FEPlotActualSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicSoluteDomain* pbd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	if (pbd)
	{
		for (int i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_ca[0];
			}
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}
*/

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
			if (psd->GetID() == nsol) return psd->GetID();
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
			if (psd->GetID() == nsol) return psd->GetID();
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
// find the local solute ID, given a global ID. If the material is not a 
// biphasic-solute, triphasic, or multiphasic material, this returns -1.
int GetLocalSoluteID(FEMaterial* pm, int nsol)
{
	// figure out the solute ID to export. This depends on the material type.
	int nsid = -1;
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (pm);
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == nsol);
		if (!present) return false;
		nsid = 0;
	}

	FETriphasic* ptm = dynamic_cast<FETriphasic*> (pm);
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		if      (ptm->m_pSolute[0]->GetSoluteID() == nsol) nsid = 0;
		else if (ptm->m_pSolute[1]->GetSoluteID() == nsol) nsid = 1;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (pm);
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		nsid = pmm->FindLocalSoluteID(nsol);
	}
	return nsid;
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
FEPlotActualSoluteConcentration::FEPlotActualSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_ARRAY, FMT_ITEM)
{
	m_pfem = pfem;
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
bool FEPlotActualSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEMaterial* pm = dom.GetMaterial();
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int negs = 0;
	for (int i = 0; i<(int)m_sol.size(); ++i)
	{
		lid[i] = GetLocalSoluteID(pm, m_sol[i]);
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
					FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());

					if (pt) ew += pt->m_ca[nsid];
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
FEPlotActualSoluteConcentration_old::FEPlotActualSoluteConcentration_old(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsol = 0;
}

//-----------------------------------------------------------------------------
// Resolve solute by name
bool FEPlotActualSoluteConcentration_old::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotActualSoluteConcentration_old::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotActualSoluteConcentration_old::Save(FEDomain &dom, FEDataStream& a)
{
	// figure out the solute ID to export. This depends on the material type.
	int nsid = GetLocalSoluteID(dom.GetMaterial(), m_nsol);

	// make sure we have a valid index
	if (nsid == -1) return false;

	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
			
		// calculate average concentration
		double ew = 0;
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
			if (pt) ew += pt->m_ca[nsid];
		}
			
		ew /= el.GaussPoints();
			
		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotActualSolConcentration_::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == m_nsol);
		if (!present) return false;

		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_ca[0];
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}

	FETriphasic* ptm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		int sid = -1;
		if (ptm->m_pSolute[0]->GetSoluteID() == m_nsol) sid = 0;
		else if (ptm->m_pSolute[1]->GetSoluteID() == m_nsol) sid = 1;
		if (sid == -1) return false;

		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_ca[sid];
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		int sid = -1;
		for (int i=0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == m_nsol) {sid = i; break;}
		if (sid == -1) return false;
		
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_ca[sid];
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
	return false;
}

//=================================================================================================
//-----------------------------------------------------------------------------
FEPlotSoluteFlux::FEPlotSoluteFlux(FEModel* pfem) : FEDomainData(PLT_ARRAY_VEC3F, FMT_ITEM)
{
	m_pfem = pfem;
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
bool FEPlotSoluteFlux::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	for (int i = 0; i<(int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
	}

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

					if (pt) ew += pt->m_j[nsid];
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
FEPlotSoluteFlux_old::FEPlotSoluteFlux_old(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM)
{
	m_nsol = 0;
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSoluteFlux_old::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve sbm by solute ID
bool FEPlotSoluteFlux_old::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotSoluteFlux_old::Save(FEDomain &dom, FEDataStream& a)
{
	// figure out the solute ID to export. This depends on the material type.
	int nsid = GetLocalSoluteID(dom.GetMaterial(), m_nsol);

	// make sure we have a valid index
	if (nsid == -1) return false;

	for (int i=0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
			
		// calculate average flux
		vec3d ew = vec3d(0,0,0);
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
			
			if (pt) ew += pt->m_j[nsid];
		}
			
		ew /= el.GaussPoints();
			
		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSolFlux_::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == m_nsol);
		if (!present) return false;
		
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_j[0];
			}
			
			ew /= el.GaussPoints();

			a << ew;
		}
		return true;
	}

	FETriphasic* ptm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		int sid = -1;
		if (ptm->m_pSolute[0]->GetSoluteID() == m_nsol) sid = 0;
		else if (ptm->m_pSolute[1]->GetSoluteID() == m_nsol) sid = 1;
		if (sid == -1) return false;
		
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_j[sid];
			}
			
			ew /= el.GaussPoints();

			a << ew;
		}
		return true;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		int sid = -1;
		for (int i=0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == m_nsol) {sid = i; break;}
		if (sid == -1) return false;
		
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			
			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_j[sid];
			}
			
			ew /= el.GaussPoints();

			a << ew;
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotOsmolarity::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
    
    FESolidDomain* sdom = dynamic_cast<FESolidDomain*>(&dom);
    FEShellDomain* ldom = dynamic_cast<FEShellDomain*>(&dom);
    if (sdom && (
        dynamic_cast<FEBiphasicSoluteSolidDomain*>(&dom) ||
        dynamic_cast<FETriphasicDomain*>(&dom) ||
        dynamic_cast<FEMultiphasicSolidDomain*>(&dom))) {
        for (i=0; i<sdom->Elements(); ++i)
        {
            FESolidElement& el = sdom->Element(i);
            
            // calculate average concentration
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (pt)
                    for (int isol=0; isol<(int)pt->m_ca.size(); ++isol)
                        ew += pt->m_ca[isol];
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
    else if (ldom && (
             dynamic_cast<FEBiphasicSoluteShellDomain*>(&dom) ||
             dynamic_cast<FEMultiphasicShellDomain*>(&dom))) {
        for (i=0; i<ldom->Elements(); ++i)
        {
            FEShellElement& el = ldom->Element(i);
            
            // calculate average concentration
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (pt)
                    for (int isol=0; isol<(int)pt->m_ca.size(); ++isol)
                        ew += pt->m_ca[isol];
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//=================================================================================================
// FEPlotSBMConcentration
//=================================================================================================

//-----------------------------------------------------------------------------
FEPlotSBMConcentration::FEPlotSBMConcentration(FEModel* pfem) : FEDomainData(PLT_ARRAY, FMT_ITEM)
{
	m_pfem = pfem;

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
// FEPlotSBMConcentration_old
//=================================================================================================

//-----------------------------------------------------------------------------
FEPlotSBMConcentration_old::FEPlotSBMConcentration_old(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsbm = 0;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBMConcentration_old::SetFilter(const char* sz)
{
	m_nsbm = GetSBMID(*m_pfem, sz);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
// Resolve sbm by solute ID
bool FEPlotSBMConcentration_old::SetFilter(int nsol)
{
	m_nsbm = GetSBMID(*m_pfem, nsol);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMConcentration_old::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the sbm ID to export. This depends on the material type.
	int nsid = GetLocalSBMID(pm, m_nsbm);

	// make sure we have a valid index
	if (nsid == -1) return false;

	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
			
		// calculate average concentration
		double ew = 0;
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
			if (pt) ew += pm->SBMConcentration(mp, nsid);
		}
			
		ew /= el.GaussPoints();
			
		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSBMConcentration_::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (pmd)
	{
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
		// Check if this solid-bound molecule is present in this specific multiphasic mixture
		int sid = -1;
		for (i=0; i<pm->SBMs(); ++i)
			if (pm->GetSBM(i)->GetSBMID() == m_nsbm) {sid = i; break;}
		if (sid == -1) return false;
		
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += pm->SBMConcentration(mp,sid);
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
        // Check if this solid-bound molecule is present in this specific multiphasic mixture
        int sid = -1;
        for (i=0; i<pm->SBMs(); ++i)
            if (pm->GetSBM(i)->GetSBMID() == m_nsbm) {sid = i; break;}
        if (sid == -1) return false;
        
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average concentration
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (st) ew += pm->SBMConcentration(mp,sid);
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElectricPotential::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average electric potential
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_psi;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
	else if (pmd)
	{
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average electric potential
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_psi;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average electric potential
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (pt) ew += pt->m_psi;
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentDensity::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	vec3d ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_Ie;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
	else if (pmd)
	{
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_Ie;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average flux
            ew = vec3d(0,0,0);
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (pt) ew += pt->m_Ie;
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialSolidVolumeFraction::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (pmd)
	{
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				
				if (pt) ew += pt->m_phi0;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average concentration
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
                
                if (pt) ew += pt->m_phi0;
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFixedChargeDensity::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average electric potential
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_cF;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
	else if (pmd)
	{
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average electric potential
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_cF;
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average electric potential
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (pt) ew += pt->m_cF;
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
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average electric potential
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMaterialPoint*  ept = (mp.ExtractData<FEElasticMaterialPoint >());
                FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				FESolutesMaterialPoint*  spt = (mp.ExtractData<FESolutesMaterialPoint >());
				
				if (spt) ew += (ept->m_J - bpt->m_phi0)*spt->m_cF/(1 - bpt->m_phi0);
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
	else if (pmd)
	{
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average electric potential
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMaterialPoint*  ept = (mp.ExtractData<FEElasticMaterialPoint >());
                FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				FESolutesMaterialPoint*  spt = (mp.ExtractData<FESolutesMaterialPoint >());
				
				if (spt) ew += (ept->m_J - bpt->m_phi0)*spt->m_cF/(1 - bpt->m_phi0);
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average electric potential
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEElasticMaterialPoint*  ept = (mp.ExtractData<FEElasticMaterialPoint >());
                FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
                FESolutesMaterialPoint*  spt = (mp.ExtractData<FESolutesMaterialPoint >());
                
                if (spt) ew += (ept->m_J - bpt->m_phi0)*spt->m_cF/(1 - bpt->m_phi0);
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
	if (pd || psd || ptd || pmd)
	{
		// get the pressure dof index
		int dof_p = GetFEModel()->GetDOFIndex("p");
		if (dof_p == -1) return false;

		int N = dom.Nodes();
		for (int i = 0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.get(dof_p);
		}
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

		int N = dom.Nodes();
		for (int i = 0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.get(dof_q);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
FEPlotEffectiveSoluteConcentration_old::FEPlotEffectiveSoluteConcentration_old(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_NODE)
{
	m_pfem = pfem;
	m_nsol = 0;
}

//-----------------------------------------------------------------------------
// Resolve solute by name
bool FEPlotEffectiveSoluteConcentration_old::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotEffectiveSoluteConcentration_old::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSoluteConcentration_old::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	int nsid = GetLocalSoluteID(dom.GetMaterial(), m_nsol);
	if (nsid == -1) return false;

	// get the dof
	const int dof_C = GetFEModel()->GetDOFIndex("concentration", nsid);

	int N = dom.Nodes();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = dom.Node(i);
		a << node.get(dof_C);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSolConcentration_::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicSolute* pbm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
	if (pbm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (pbm->GetSolute()->GetSoluteID() == m_nsol);
		if (!present) return false;

		// get the dof
		const int dof_C = GetFEModel()->GetDOFIndex("concentration", m_nsol);

		int N = dom.Nodes();
		for (int i = 0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.get(dof_C);
		}
		return true;
	}

	FETriphasic* ptm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		bool present = (ptm->m_pSolute[0]->GetSoluteID() == m_nsol) || (ptm->m_pSolute[1]->GetSoluteID() == m_nsol);
		if (!present) return false;

		// get the dof
		const int dof_C = GetFEModel()->GetDOFIndex("concentration", m_nsol);

		int N = dom.Nodes();
		for (int i = 0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.get(dof_C);
		}
		return true;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		bool present = false;
		for (int i = 0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == m_nsol) { present = true; break; }
		if (!present) return false;

		// get the dof
		const int dof_C = GetFEModel()->GetDOFIndex("concentration", m_nsol);

		int N = dom.Nodes();
		for (int i = 0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.get(dof_C);
		}
		return true;
	}
	return false;
}

//=================================================================================================
// FEPlotEffectiveSoluteConcentration
//=================================================================================================

FEPlotEffectiveSoluteConcentration::FEPlotEffectiveSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_ARRAY, FMT_NODE)
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

	FEMaterial* pm = dom.GetMaterial();
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int negs = 0;
	for (int i = 0; i<(int)m_sol.size(); ++i)
	{
		lid[i] = GetLocalSoluteID(pm, m_sol[i]);
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
FEPlotEffectiveShellSoluteConcentration::FEPlotEffectiveShellSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_NODE)
{
	m_pfem = pfem;
	m_nsol = 0;
}

//-----------------------------------------------------------------------------
// Resolve solute by name
bool FEPlotEffectiveShellSoluteConcentration::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotEffectiveShellSoluteConcentration::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveShellSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicDomain* pbsd = dynamic_cast<FEBiphasicDomain*>(&dom);
	FEBiphasicSoluteDomain* pbssd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FEMultiphasicDomain* pmpsd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (pbsd || pbssd || pmpsd)
	{
		// make sure we have a valid index
		int nsid = GetLocalSoluteID(dom.GetMaterial(), m_nsol);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveShellSolConcentration_::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicDomain* pbsd = dynamic_cast<FEBiphasicDomain*>(&dom);
	FEBiphasicSoluteDomain* pbssd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FEMultiphasicDomain* pmpsd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (pbsd || pbssd || pmpsd)
	{
		FEBiphasicSolute* pbm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
		if (pbm)
		{
			// Check if this solute is present in this specific biphasic-solute mixture
			bool present = (pbm->GetSolute()->GetSoluteID() == m_nsol);
			if (!present) return false;

			// get the dof
			const int dof_D = GetFEModel()->GetDOFIndex("shell concentration", m_nsol);
			if (dof_D == -1) return false;

			int N = dom.Nodes();
			for (int i = 0; i<N; ++i)
			{
				FENode& node = dom.Node(i);
				a << node.get(dof_D);
			}
			return true;
		}

		FETriphasic* ptm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
		if (ptm)
		{
			// Check if this solute is present in this specific triphasic mixture
			bool present = (ptm->m_pSolute[0]->GetSoluteID() == m_nsol) || (ptm->m_pSolute[1]->GetSoluteID() == m_nsol);
			if (!present) return false;

			// get the dof
			const int dof_D = GetFEModel()->GetDOFIndex("shell concentration", m_nsol);
			if (dof_D == -1) return false;

			int N = dom.Nodes();
			for (int i = 0; i<N; ++i)
			{
				FENode& node = dom.Node(i);
				a << node.get(dof_D);
			}
			return true;
		}

		FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
		if (pmm)
		{
			// Check if this solute is present in this specific multiphasic mixture
			bool present = false;
			for (int i = 0; i<pmm->Solutes(); ++i)
				if (pmm->GetSolute(i)->GetSoluteID() == m_nsol) { present = true; break; }
			if (!present) return false;

			// get the dof
			const int dof_D = GetFEModel()->GetDOFIndex("shell concentration", m_nsol);
			if (dof_D == -1) return false;

			int N = dom.Nodes();
			for (int i = 0; i<N; ++i)
			{
				FENode& node = dom.Node(i);
				a << node.get(dof_D);
			}
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReceptorLigandConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FEBiphasicSoluteSolidDomain* pbd = dynamic_cast<FEBiphasicSoluteSolidDomain*>(&dom);
    FEBiphasicSoluteShellDomain* psd = dynamic_cast<FEBiphasicSoluteShellDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_sbmr[0];
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average concentration
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (pt) ew += pt->m_sbmr[0];
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
}

//=================================================================================================

FEPlotSBMRefAppDensity::FEPlotSBMRefAppDensity(FEModel* pfem) : FEDomainData(PLT_ARRAY, FMT_ITEM)
{
	m_pfem = pfem;

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

//=================================================================================================

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBMRefAppDensity_old::SetFilter(const char* sz)
{
	m_nsbm = GetSBMID(*m_pfem, sz);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
// Resolve sbm by solute ID
bool FEPlotSBMRefAppDensity_old::SetFilter(int nsol)
{
	m_nsbm = GetSBMID(*m_pfem, nsol);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMRefAppDensity_old::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the sbm ID to export. This depends on the material type.
	int nsid = GetLocalSBMID(pm, m_nsbm);
	if (nsid == -1) return false;

	for (int i=0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
			
		// calculate average concentration
		double ew = 0;
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
			
			if (st) ew += st->m_sbmr[nsid];
		}
		ew /= el.GaussPoints();
			
		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSBMRefAppDensity_::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&dom);
    FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&dom);
	if (pmd)
	{
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
		// Check if this solid-bound molecule is present in this specific multiphasic mixture
		int sid = -1;
		for (i=0; i<pm->SBMs(); ++i)
			if (pm->GetSBM(i)->GetSBMID() == m_nsbm) {sid = i; break;}
		if (sid == -1) return false;
		
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_sbmr[sid];
			}
			
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
    else if (psd)
    {
        FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
        // Check if this solid-bound molecule is present in this specific multiphasic mixture
        int sid = -1;
        for (i=0; i<pm->SBMs(); ++i)
            if (pm->GetSBM(i)->GetSBMID() == m_nsbm) {sid = i; break;}
        if (sid == -1) return false;
        
        for (i=0; i<psd->Elements(); ++i)
        {
            FEShellElement& el = psd->Element(i);
            
            // calculate average concentration
            ew = 0;
            for (j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
                
                if (st) ew += st->m_sbmr[sid];
            }
            
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
	return false;
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


//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
// Plot contact gap
bool FEPlotPressureGap::Save(FESurface& surf, FEDataStream& a)
{
	FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	double gn[MFN];
	a.assign(MFN*NF, 0.f);
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& f = pcs->Element(i);
		pcs->GetNodalPressureGap(i, gn);
		int ne = (int)f.m_lnode.size();
		for (int j = 0; j< ne; ++j) a[MFN*i + j] = (float) gn[j];
	}
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
bool FEPlotFluidForce2::Save(FESurface &surf, FEDataStream &a)
{
	// get number of facets
	int NF = surf.Elements();
	if (NF == 0) return false;

	// this assumes that the surface sits on top of a single domain
	// so that we can figure out the domain from a single element
	FESurfaceElement& ref = surf.Element(0);
	if (ref.m_elem[0] <= 0) return false;

	// get the element
	FEMesh& mesh = *surf.GetMesh();
	FEElement* el = mesh.FindElementFromID(ref.m_elem[0]);
	if (el == 0) return false;

	// get the domain this element belongs to
	FEDomain* dom = el->GetDomain();
	if (dom == 0) return false;

	// see if this is a biphasic domain
	FEBiphasicSolidDomain* biphasicDomain = dynamic_cast<FEBiphasicSolidDomain*>(dom);
	if (biphasicDomain == 0) return false;

	// The biphasic solid domain contains actual nodal pressures.
	// we want to evaluate this over the surface. 
	vector<double> nodalPressures;
	biphasicDomain->GetNodalPressures(nodalPressures);

	// loop over all surfaces facets
	double pn[FEElement::MAX_NODES];
	vec3d F(0, 0, 0);
	for (int i = 0; i<NF; ++i)
	{
		FESurfaceElement& face = surf.Element(i);

		int nint = face.GaussPoints();
		int neln = face.Nodes();

		// nodal pressures
		for (int j = 0; j<neln; ++j) pn[j] = nodalPressures[face.m_node[j]];

		// evaluate the fluid force for that element
		for (int j = 0; j<nint; ++j)
		{
			// get the base vectors
			vec3d g[2];
			surf.CoBaseVectors(face, j, g);

			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];

			// gauss weight
			double w = face.GaussWeights()[j];

			// fluid pressure
			double p = face.eval(pn, j);
			
			// contact force
			F += n*(w*p);
		}
	}

	// store the force
	a << F;

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
