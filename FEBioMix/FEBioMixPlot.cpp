#include "stdafx.h"
#include "FEBioMixPlot.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicSoluteDomain.h"
#include "FETriphasicDomain.h"
#include "FEMultiphasicDomain.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FEBiphasicContactSurface.h"
#include "FEBioPlot/FEBioPlotFile.h"

//-----------------------------------------------------------------------------
bool FEPlotActualFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&bd)) || 
		(dynamic_cast<FEBiphasicSoluteDomain*>(&bd)) ||
		(dynamic_cast<FETriphasicDomain*     >(&bd)) ||
		(dynamic_cast<FEMultiphasicDomain*   >(&bd)))
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidFlux::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&bd)) || 
		(dynamic_cast<FEBiphasicSoluteDomain*>(&bd)) ||
		(dynamic_cast<FETriphasicDomain*     >(&bd)) ||
		(dynamic_cast<FEMultiphasicDomain*   >(&bd)))
	{
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);

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
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&bd)) ||
		(dynamic_cast<FEBiphasicSoluteDomain*>(&bd)) ||
		(dynamic_cast<FETriphasicDomain*     >(&bd)) ||
		(dynamic_cast<FEMultiphasicDomain*   >(&bd)))
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
			if (psd->m_nID == nsol) return psd->m_nID;
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
			if (psd->m_nID == nsol) return psd->m_nID;
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
// Finds the solute ID given the name of the solute
int GetSoluteID(FEModel& fem, const char* sz)
{
	// find the solute with that name
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd)
		{
			if (strcmp(psd->m_szname, sz) == 0) return psd->m_nID;
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
// Finds the sbm ID given the name of the sbm
int GetSBMID(FEModel& fem, const char* sz)
{
	// find the solute with that name
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESBMData* psd = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
		if (psd)
		{
			if (strcmp(psd->m_szname, sz) == 0) return psd->m_nID;
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
		for (int i=0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == nsol) {nsid = i; break;}
	}
	return nsid;
}

//-----------------------------------------------------------------------------
// find the local SBM ID, given a global ID. If the material is not a 
// multiphasic material, this returns -1.
int GetLocalSBMID(FEMaterial* pm, int nsbm)
{
	// figure out the SBM ID to export. This depends on the material type.
	int nsid = -1;
	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (pm);
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		for (int i=0; i<pmm->SBMs(); ++i)
			if (pmm->GetSBM(i)->GetSBMID() == nsbm) {nsid = i; break;}
	}
	return nsid;
}

//-----------------------------------------------------------------------------
FEPlotActualSoluteConcentration::FEPlotActualSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsol = 0;
}

//-----------------------------------------------------------------------------
// Resolve solute by name
bool FEPlotActualSoluteConcentration::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotActualSoluteConcentration::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotActualSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
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

//-----------------------------------------------------------------------------
FEPlotSoluteFlux::FEPlotSoluteFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM)
{
	m_nsol = 0;
	m_pfem = 0;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSoluteFlux::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve sbm by solute ID
bool FEPlotSoluteFlux::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotSoluteFlux::Save(FEDomain &dom, FEDataStream& a)
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
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (psd)
	{
		for (i=0; i<psd->Elements(); ++i)
		{
			FESolidElement& el = psd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
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
	else if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_ca[0] + pt->m_ca[1];
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

//-----------------------------------------------------------------------------
FEPlotSBMConcentration::FEPlotSBMConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsbm = 0;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBMConcentration::SetFilter(const char* sz)
{
	m_nsbm = GetSBMID(*m_pfem, sz);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
// Resolve sbm by solute ID
bool FEPlotSBMConcentration::SetFilter(int nsol)
{
	m_nsbm = GetSBMID(*m_pfem, nsol);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the sbm ID to export. This depends on the material type.
	int nsid = GetLocalSBMID(dom.GetMaterial(), m_nsbm);

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
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElectricPotential::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentDensity::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	vec3d ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialSolidVolumeFraction::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFixedChargeDensity::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialFixedChargeDensity::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	FEBiphasicSolidDomain*  pd  = dynamic_cast<FEBiphasicSolidDomain* >(&dom);
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain*      ptd = dynamic_cast<FETriphasicDomain*     >(&dom);
	FEMultiphasicDomain*    pmd = dynamic_cast<FEMultiphasicDomain*   >(&dom);
	if (pd || psd || ptd || pmd)
	{
		int N = dom.Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.m_pt;
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
FEPlotEffectiveSoluteConcentration::FEPlotEffectiveSoluteConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_NODE)
{
	m_pfem = pfem;
	m_nsol = 0;
}

//-----------------------------------------------------------------------------
// Resolve solute by name
bool FEPlotEffectiveSoluteConcentration::SetFilter(const char* sz)
{
	m_nsol = GetSoluteID(*m_pfem, sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotEffectiveSoluteConcentration::SetFilter(int nsol)
{
	m_nsol = GetSoluteID(*m_pfem, nsol);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSoluteConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	int nsid = GetLocalSoluteID(dom.GetMaterial(), m_nsol);

	// make sure we have a valid index
	if (nsid == -1) return false;
		
	int N = dom.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = dom.Node(i);
		a << node.m_ct[nsid];
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
		
		int N = dom.Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.m_ct[m_nsol];
		}
		return true;
	}
	
	FETriphasic* ptm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		bool present = (ptm->m_pSolute[0]->GetSoluteID() == m_nsol) || (ptm->m_pSolute[1]->GetSoluteID() == m_nsol);
		if (!present) return false;
		
		int N = dom.Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.m_ct[m_nsol];
		}
		return true;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		bool present = false;
		for (int i=0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == m_nsol) {present = true; break;}
		if (!present) return false;
		
		int N = dom.Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = dom.Node(i);
			a << node.m_ct[m_nsol];
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReceptorLigandConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	int i, j;
	double ew;
	FEBiphasicSoluteDomain* pbd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBMRefAppDensity::SetFilter(const char* sz)
{
	m_nsbm = GetSBMID(*m_pfem, sz);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
// Resolve sbm by solute ID
bool FEPlotSBMRefAppDensity::SetFilter(int nsol)
{
	m_nsbm = GetSBMID(*m_pfem, nsol);
	return (m_nsbm != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotSBMRefAppDensity::Save(FEDomain &dom, FEDataStream& a)
{
	FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the sbm ID to export. This depends on the material type.
	int nsid = GetLocalSBMID(dom.GetMaterial(), m_nsbm);
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
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveElasticity::Save(FEDomain &dom, FEDataStream& a)
{
    tens4ds c;
    
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESolidDomain* pbd = static_cast<FESolidDomain*>(&dom);
    
    FEBiphasic*       pb  = dynamic_cast<FEBiphasic      *> (dom.GetMaterial());
    FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
    FETriphasic*      ptp = dynamic_cast<FETriphasic     *> (dom.GetMaterial());
    FEMultiphasic*    pmp = dynamic_cast<FEMultiphasic   *> (dom.GetMaterial());
    if ((pb == 0) && (pbs == 0) && (ptp == 0) && (pmp == 0)) return false;

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
