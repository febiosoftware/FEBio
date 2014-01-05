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
bool FEPlotActualFluidPressure::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&dom)) || 
		(dynamic_cast<FEBiphasicSoluteDomain*>(&dom)) ||
		(dynamic_cast<FETriphasicDomain*     >(&dom)) ||
		(dynamic_cast<FEMultiphasicDomain*   >(&dom)))
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average pressure
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				
				if (pt) ew += pt->m_pa;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidFlux::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
	vec3d ew;
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&dom)) || 
		(dynamic_cast<FEBiphasicSoluteDomain*>(&dom)) ||
		(dynamic_cast<FETriphasicDomain*     >(&dom)) ||
		(dynamic_cast<FEMultiphasicDomain*   >(&dom)))
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);

			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());

				if (pt) ew += pt->m_w;
			}

			ew /= el.GaussPoints();

			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;

			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
		}
		return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotNodalFluidFlux::Save(FEDomain &dom, vector<float>& a)
{
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if ((dynamic_cast<FEBiphasicSolidDomain* >(&dom)) ||
		(dynamic_cast<FEBiphasicSoluteDomain*>(&dom)) ||
		(dynamic_cast<FETriphasicDomain*     >(&dom)) ||
		(dynamic_cast<FEMultiphasicDomain*   >(&dom)))
	{
		for (int i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);

			int nint = el.GaussPoints();
			int neln = el.Nodes();
			assert(nint == neln); // TODO: just for now

			// fluid flux at gauss points
			int j;
			double vi[3][FEElement::MAX_NODES];
			for (j=0; j<nint; ++j)
			{
				FEBiphasicMaterialPoint* pt = el.m_State[j]->ExtractData<FEBiphasicMaterialPoint>(); assert(pt);
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
			for (j=0; j<neln; ++j)
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

//-----------------------------------------------------------------------------
bool FEPlotActualSoluteConcentration::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
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

//-----------------------------------------------------------------------------
bool FEPlotActualSolConcentration_::Save(FEDomain &dom, vector<float>& a)
{
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == m_nsol);
		if (!present) return false;

		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
		for (int i=0; i<sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_ca[0];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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

		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
		for (int i=0; i<sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_ca[sid];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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
		
		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
		for (int i=0; i<sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			
			// calculate average concentration
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_ca[sid];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSoluteFlux::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
	vec3d ew;
	FEBiphasicSoluteDomain* pbd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_j[0];
			}
			
			ew /= el.GaussPoints();
			
			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;
			
			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSolFlux_::Save(FEDomain &dom, vector<float>& a)
{
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == m_nsol);
		if (!present) return false;
		
		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
		for (int i=0; i<sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			
			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_j[0];
			}
			
			ew /= el.GaussPoints();
			
			float af[3];
			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;
			
			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
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
		
		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
		for (int i=0; i<sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			
			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_j[sid];
			}
			
			ew /= el.GaussPoints();
			
			float af[3];
			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;
			
			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
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
		
		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
		for (int i=0; i<sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			
			// calculate average flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_j[sid];
			}
			
			ew /= el.GaussPoints();
			
			float af[3];
			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;
			
			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotOsmolarity::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_ca[0];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_ca[0] + pt->m_ca[1];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt)
                    for (int isol=0; isol<(int)pt->m_ca.size(); ++isol)
                        ew += pt->m_ca[isol];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSBMConcentration_::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += pm->SBMConcentration(mp,sid);
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElectricPotential::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_psi;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_psi;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentDensity::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_Ie;
			}
			
			ew /= el.GaussPoints();
			
			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;
			
			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_Ie;
			}
			
			ew /= el.GaussPoints();
			
			af[0] = (float) ew.x;
			af[1] = (float) ew.y;
			af[2] = (float) ew.z;
			
			a.push_back(af[0]);
			a.push_back(af[1]);
			a.push_back(af[2]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialSolidVolumeFraction::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				
				if (pt) ew += pt->m_phi0;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFixedChargeDensity::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_cF;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_cF;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReferentialFixedChargeDensity::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
                FEElasticMaterialPoint* ept = (mp.ExtractData<FEElasticMaterialPoint>());
                FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				FESolutesMaterialPoint* spt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (spt) ew += (ept->m_J - bpt->m_phi0)*spt->m_cF/(1 - bpt->m_phi0);
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
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
				FEMaterialPoint& mp = *el.m_State[j];
                FEElasticMaterialPoint* ept = (mp.ExtractData<FEElasticMaterialPoint>());
                FEBiphasicMaterialPoint* bpt = (mp.ExtractData<FEBiphasicMaterialPoint>());
				FESolutesMaterialPoint* spt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (spt) ew += (ept->m_J - bpt->m_phi0)*spt->m_cF/(1 - bpt->m_phi0);
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveFluidPressure::Save(FEDomain &dom, vector<float>& a)
{
	FEBiphasicSolidDomain* pd = dynamic_cast<FEBiphasicSolidDomain*>(&dom);
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (pd)
	{
		int N = pd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pd->Node(i);
			a.push_back((float) node.m_pt);
		}
		return true;
	}
	else if (psd)
	{
		int N = psd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = psd->Node(i);
			a.push_back((float) node.m_pt);
		}
		return true;
	}
	else if (ptd)
	{
		int N = ptd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = ptd->Node(i);
			a.push_back((float) node.m_pt);
		}
		return true;
	}
	else if (pmd)
	{
		int N = pmd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pmd->Node(i);
			a.push_back((float) node.m_pt);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSoluteConcentration::Save(FEDomain &dom, vector<float>& a)
{
	FEBiphasicSoluteDomain* pd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	if (pd)
	{
		int N = pd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pd->Node(i);
			a.push_back((float) node.m_ct[0]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSolConcentration_::Save(FEDomain &dom, vector<float>& a)
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
			a.push_back((float) node.m_ct[m_nsol]);
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
			a.push_back((float) node.m_ct[m_nsol]);
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
			a.push_back((float) node.m_ct[m_nsol]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotReceptorLigandConcentration::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (pt) ew += pt->m_sbmr[0];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSBMRefAppDensity_::Save(FEDomain &dom, vector<float>& a)
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
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_sbmr[sid];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}


//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
// Plot contact gap
bool FEPlotPressureGap::Save(FESurface& surf, vector<float>& a)
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
bool FEPlotFluidForce::Save(FESurface &surf, std::vector<float> &a)
{
	FEBiphasicContactSurface* pcs = dynamic_cast<FEBiphasicContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(3*MFN*NF, 0.f);
	vec3d fn = pcs->GetFluidForce();
	for (int j=0; j<NF; ++j)
	{
		FESurfaceElement& el = pcs->Element(j);
        
		// store in archive
		int ne = el.Nodes();
		for (int k=0; k<ne; ++k)
		{
			a[3*MFN*j +3*k   ] = (float) fn.x;
			a[3*MFN*j +3*k +1] = (float) fn.y;
			a[3*MFN*j +3*k +2] = (float) fn.z;
		}
	}
    
	return true;
}
