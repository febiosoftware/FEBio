#include "stdafx.h"
#include "FEPlotDomainData.h"
#include "FEBioLib/FEDamageNeoHookean.h"
#include "FEBioLib/FEDamageTransIsoMooneyRivlin.h"
#include "FEBioLib/FEBiphasicSoluteDomain.h"
#include "FEBioLib/FEBiphasicSolidDomain.h"
#include "FEBioLib/FETriphasicDomain.h"
#include "FEBioLib/FERigidSolidDomain.h"
#include "FEBioLib/FERigidShellDomain.h"
#include "FEBioLib/FEElasticMixture.h"
#include "FEBioLib/FEBiphasicSolute.h"
#include "FEBioLib/FETriphasic.h"
#include "FEBioLib/FEMultiphasicDomain.h"
#include "FEBioLib/FEMultiphasic.h"
#include "FEBioLib/FEUT4Domain.h"
#include <FEBioLib/FEHeatSolidDomain.h>
#include <FEBioLib/FEHeatTransferMaterial.h>

//-----------------------------------------------------------------------------
//! Store the average stresses for each element. 
bool FEPlotElementStress::Save(FEDomain& dom, vector<float>& a)
{
	// write solid stresses
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (pbd) return WriteSolidStress(*pbd, a);

	FELinearSolidDomain* pbl = dynamic_cast<FELinearSolidDomain*>(&dom);
	if (pbl) return WriteLinearSolidStress(*pbl, a);

	// write shell stresses
	FEElasticShellDomain* pbs = dynamic_cast<FEElasticShellDomain*>(&dom);
	if (pbs) return WriteShellStress(*pbs, a);

	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElementStress::WriteSolidStress(FEElasticSolidDomain& d, vector<float>& a)
{
	// make sure this is not a rigid body
	if (dynamic_cast<FERigidSolidDomain*>(&d)) return false;

	// write solid element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FESolidElement& el = d.Element(i);

		float s[6] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

			if (ppt)
			{
				FEElasticMaterialPoint& pt = *ppt;
				s[0] += (float) (f*pt.s.xx());
				s[1] += (float) (f*pt.s.yy());
				s[2] += (float) (f*pt.s.zz());
				s[3] += (float) (f*pt.s.xy());
				s[4] += (float) (f*pt.s.yz());
				s[5] += (float) (f*pt.s.xz());
			}
		}

		a.push_back(s[0]);
		a.push_back(s[1]);
		a.push_back(s[2]);
		a.push_back(s[3]);
		a.push_back(s[4]);
		a.push_back(s[5]);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotElementStress::WriteShellStress(FEElasticShellDomain& d, vector<float>& a)
{
	// make sure this is not a rigid body
	if (dynamic_cast<FERigidShellDomain*>(&d)) return false;

	// write shell element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FEShellElement& el = d.Element(i);

		float s[6] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

			if (ppt)
			{
				FEElasticMaterialPoint& pt = *ppt;
				s[0] += (float) (f*pt.s.xx());
				s[1] += (float) (f*pt.s.yy());
				s[2] += (float) (f*pt.s.zz());
				s[3] += (float) (f*pt.s.xy());
				s[4] += (float) (f*pt.s.yz());
				s[5] += (float) (f*pt.s.xz());
			}
		}

		a.push_back(s[0]);
		a.push_back(s[1]);
		a.push_back(s[2]);
		a.push_back(s[3]);
		a.push_back(s[4]);
		a.push_back(s[5]);
	}

	return true;
}


//-----------------------------------------------------------------------------
bool FEPlotElementStress::WriteLinearSolidStress(FELinearSolidDomain& d, vector<float>& a)
{
	// write solid element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FESolidElement& el = d.Element(i);

		float s[6] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

			if (ppt)
			{
				FEElasticMaterialPoint& pt = *ppt;
				mat3ds& es = pt.s;
				s[0] += (float) (f*es.xx());
				s[1] += (float) (f*es.yy());
				s[2] += (float) (f*es.zz());
				s[3] += (float) (f*es.xy());
				s[4] += (float) (f*es.yz());
				s[5] += (float) (f*es.xz());
			}
		}

		a.push_back(s[0]);
		a.push_back(s[1]);
		a.push_back(s[2]);
		a.push_back(s[3]);
		a.push_back(s[4]);
		a.push_back(s[5]);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRelativeVolume::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average flux
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEElasticMaterialPoint* pt = (mp.ExtractData<FEElasticMaterialPoint>());
				
				if (pt) ew += pt->J;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

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
				FESoluteMaterialPoint* pt = (mp.ExtractData<FESoluteMaterialPoint>());
				
				if (pt) ew += pt->m_ca;
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
				FESoluteMaterialPoint* pt = (mp.ExtractData<FESoluteMaterialPoint>());
				
				if (pt) ew += pt->m_j;
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
bool FEPlotActualSolConcentration_::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (psd)
	{
		FEBiphasicSolute* pm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (pm->m_pSolute->GetSoluteID() == m_nsol);
		if (!present) return false;

		for (i=0; i<psd->Elements(); ++i)
		{
			FESolidElement& el = psd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESoluteMaterialPoint* pt = (mp.ExtractData<FESoluteMaterialPoint>());
				
				if (pt) ew += pt->m_ca;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	else if (ptd)
	{
		FETriphasic* pm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
		// Check if this solute is present in this specific triphasic mixture
		int sid = -1;
		if (pm->m_pSolute[0]->GetSoluteID() == m_nsol) sid = 0;
		else if (pm->m_pSolute[1]->GetSoluteID() == m_nsol) sid = 1;
		if (sid == -1) return false;

		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESaltMaterialPoint* st = (mp.ExtractData<FESaltMaterialPoint>());
				
				if (st) ew += st->m_ca[sid];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	else if (pmd)
	{
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
		// Check if this solute is present in this specific multiphasic mixture
		int sid = -1;
		for (i=0; i<(int)pm->m_pSolute.size(); ++i)
			if (pm->m_pSolute[i]->GetSoluteID() == m_nsol) {sid = i; break;}
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
bool FEPlotSolFlux_::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
	vec3d ew;
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (psd)
	{
		FEBiphasicSolute* pm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (pm->m_pSolute->GetSoluteID() == m_nsol);
		if (!present) return false;
		
		for (i=0; i<psd->Elements(); ++i)
		{
			FESolidElement& el = psd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESoluteMaterialPoint* pt = (mp.ExtractData<FESoluteMaterialPoint>());
				
				if (pt) ew += pt->m_j;
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
	else if (ptd)
	{
		FETriphasic* pm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
		// Check if this solute is present in this specific triphasic mixture
		int sid = -1;
		if (pm->m_pSolute[0]->GetSoluteID() == m_nsol) sid = 0;
		else if (pm->m_pSolute[1]->GetSoluteID() == m_nsol) sid = 1;
		if (sid == -1) return false;
		
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESaltMaterialPoint* st = (mp.ExtractData<FESaltMaterialPoint>());
				
				if (st) ew += st->m_j[sid];
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
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
		// Check if this solute is present in this specific multiphasic mixture
		int sid = -1;
		for (i=0; i<(int)pm->m_pSolute.size(); ++i)
			if (pm->m_pSolute[i]->GetSoluteID() == m_nsol) {sid = i; break;}
		if (sid == -1) return false;
		
		for (i=0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());
				
				if (st) ew += st->m_j[sid];
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
				FESaltMaterialPoint* pt = (mp.ExtractData<FESaltMaterialPoint>());
				
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
				FESaltMaterialPoint* pt = (mp.ExtractData<FESaltMaterialPoint>());
				
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
				FESaltMaterialPoint* pt = (mp.ExtractData<FESaltMaterialPoint>());
				
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
bool FEPlotFiberVector::Save(FEDomain &dom, vector<float>& a)
{
	int i, j, n;
	float f[3];
	vec3d r;
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (pbd)
	{
		int BE = pbd->Elements();
		for (i=0; i<BE; ++i)
		{
			FESolidElement& el = pbd->Element(i);
			n = el.GaussPoints();
			r = vec3d(0,0,0);
			for (j=0; j<n; ++j)
			{
				FEElasticMaterialPoint& pt = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
				r.x += pt.Q[0][0];
				r.y += pt.Q[1][0];
				r.z += pt.Q[2][0];
			}
			r /= n;
			f[0] = (float) r.x;
			f[1] = (float) r.y;
			f[2] = (float) r.z;

			a.push_back(f[0]);
			a.push_back(f[1]);
			a.push_back(f[2]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store shell thicknesses
bool FEPlotShellThickness::Save(FEDomain &dom, vector<float> &a)
{
	FEShellDomain* pbs = dynamic_cast<FEShellDomain*>(&dom);
	if (pbs)
	{
		int NS = pbs->Elements();
		for (int i=0; i<NS; ++i)
		{
			FEShellElement& e = pbs->Element(i);
			int n = e.Nodes();
			for (int j=0; j<n; ++j)
			{
				vec3d D = pbs->GetMesh()->Node(e.m_node[j]).m_Dt;
				double h = e.m_h0[j] * D.norm();
				a.push_back((float) h);
			}
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
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	FEMultiphasicDomain* pmd = dynamic_cast<FEMultiphasicDomain*>(&dom);
	if (psd)
	{
		FEBiphasicSolute* pm = dynamic_cast<FEBiphasicSolute*> (dom.GetMaterial());
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (pm->m_pSolute->GetSoluteID() == m_nsol);
		if (!present) return false;
		
		int N = psd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = psd->Node(i);
			a.push_back((float) node.m_ct[m_nsol]);
		}
		return true;
	}
	else if (ptd)
	{
		FETriphasic* pm = dynamic_cast<FETriphasic*> (dom.GetMaterial());
		// Check if this solute is present in this specific triphasic mixture
		bool present = (pm->m_pSolute[0]->GetSoluteID() == m_nsol) || (pm->m_pSolute[1]->GetSoluteID() == m_nsol);
		if (!present) return false;
		
		int N = ptd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = ptd->Node(i);
			a.push_back((float) node.m_ct[m_nsol]);
		}
		return true;
	}
	else if (pmd)
	{
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (dom.GetMaterial());
		// Check if this solute is present in this specific multiphasic mixture
		bool present = false;
		for (int i=0; i<(int)pm->m_pSolute.size(); ++i)
			if (pm->m_pSolute[i]->GetSoluteID() == m_nsol) {present = true; break;}
		if (!present) return false;
		
		int N = pmd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pmd->Node(i);
			a.push_back((float) node.m_ct[m_nsol]);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotDamage::Save(FEDomain &m, vector<float>& a)
{
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&m);
	if (pbd)
	{
		FESolidDomain& d = *pbd;
		for (int i=0; i<d.Elements(); ++i)
		{
			FESolidElement& el = d.Element(i);

			float D = 0.f;
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEDamageMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEDamageMaterialPoint>());

				if (ppt)
				{
					FEDamageMaterialPoint& pt = *ppt;
					D += (float) pt.m_D;
				}

				FETIMRDamageMaterialPoint* pt2 = (el.m_State[j]->ExtractData<FETIMRDamageMaterialPoint>());
				if (pt2)
				{
					D += (float) pt2->m_Df;
				}
			}
			D /= (float) nint;
			a.push_back(1.f - D);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotMixtureVolumeFraction::Save(FEDomain &m, std::vector<float> &a)
{
	// extract the mixture material
	FEMaterial* pmat = m.GetMaterial();
	FEElasticMixture* pm = dynamic_cast<FEElasticMixture*>(pmat);
	if (pm == 0) return false;

	// store the volume fraction of the first material
	int N = m.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& e = m.ElementRef(i);

		float s = 0.f;
		int nint = e.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			FEElasticMixtureMaterialPoint& pt = *e.m_State[n]->ExtractData<FEElasticMixtureMaterialPoint>();
			s += (float) pt.m_w[0];
		}

		a.push_back(s / (float) nint);
	}

	return true;
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
				FESoluteMaterialPoint* pt = (mp.ExtractData<FESoluteMaterialPoint>());
				
				if (pt) ew += pt->m_crc;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}


//-----------------------------------------------------------------------------
bool FEPlotUT4NodalStresses::Save(FEDomain& dom, vector<float>& a)
{
	FEUT4Domain* pd = dynamic_cast<FEUT4Domain*>(&dom);
	if (pd == 0) return false;

	int N = pd->Nodes();
	for (int i=0; i<N; ++i)
	{
		FEUT4Domain::UT4NODE& n = pd->UT4Node(i);
		mat3ds& s = n.si;
		a.push_back((float) s.xx());
		a.push_back((float) s.yy());
		a.push_back((float) s.zz());
		a.push_back((float) s.xy());
		a.push_back((float) s.yz());
		a.push_back((float) s.xz());
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotHeatFlux::Save(FEDomain &dom, vector<float>& a)
{
	FEHeatSolidDomain* pbd = dynamic_cast<FEHeatSolidDomain*>(&dom);
	if (pbd)
	{
		// loop over all elements
		for (int i=0; i<pbd->Elements(); ++i)
		{
			// get the next element
			FESolidElement& el = pbd->Element(i);

			// calculate average heat flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEHeatMaterialPoint* pt = (mp.ExtractData<FEHeatMaterialPoint>());
				if (pt) ew += pt->m_q;
			}
			ew /= el.GaussPoints();

			// store to buffer
			a.push_back((float) ew.x);
			a.push_back((float) ew.y);
			a.push_back((float) ew.z);
		}
		return true;
	}

	return false;
}
