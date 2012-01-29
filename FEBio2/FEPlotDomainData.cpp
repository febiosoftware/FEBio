#include "stdafx.h"
#include "FEPlotDomainData.h"
#include "FEBioLib/FEDamageNeoHookean.h"
#include "FEBioLib/FEDamageTransIsoMooneyRivlin.h"
#include "FEBiphasicSoluteDomain.h"
#include "FEBiphasicSolidDomain.h"
#include "FETriphasicDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEBioLib/FEElasticMixture.h"
#include "FEBioLib/FEBiphasicSolute.h"
#include "FEBioLib/FETriphasic.h"

//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotEffectiveFluidPressure		, FEPlotData, "effective fluid pressure"        );
REGISTER_FEBIO_CLASS(FEPlotActualFluidPressure          , FEPlotData, "fluid pressure"                  );
REGISTER_FEBIO_CLASS(FEPlotElementStress                , FEPlotData, "stress"                          );
REGISTER_FEBIO_CLASS(FEPlotRelativeVolume               , FEPlotData, "relative volume"                 );
REGISTER_FEBIO_CLASS(FEPlotFluidFlux                    , FEPlotData, "fluid flux"                      );
REGISTER_FEBIO_CLASS(FEPlotFiberVector                  , FEPlotData, "fiber vector"                    );
REGISTER_FEBIO_CLASS(FEPlotEffectiveSoluteConcentration , FEPlotData, "effective solute concentration"  );
REGISTER_FEBIO_CLASS(FEPlotShellThickness               , FEPlotData, "shell thickness"                 );
REGISTER_FEBIO_CLASS(FEPlotActualSoluteConcentration    , FEPlotData, "solute concentration"            );
REGISTER_FEBIO_CLASS(FEPlotSoluteFlux                   , FEPlotData, "solute flux"                     );
REGISTER_FEBIO_CLASS(FEPlotDamage                       , FEPlotData, "damage"                          );
REGISTER_FEBIO_CLASS(FEPlotMixtureVolumeFraction        , FEPlotData, "volume fraction"                 );
REGISTER_FEBIO_CLASS(FEPlotReceptorLigandConcentration  , FEPlotData, "receptor-ligand concentration"   );
REGISTER_FEBIO_CLASS(FEPlotEffectiveSol1Concentration   , FEPlotData, "effective solute 1 concentration");
REGISTER_FEBIO_CLASS(FEPlotActualSol1Concentration      , FEPlotData, "solute 1 concentration"          );
REGISTER_FEBIO_CLASS(FEPlotSol1Flux                     , FEPlotData, "solute 1 flux"                   );
REGISTER_FEBIO_CLASS(FEPlotEffectiveSol2Concentration   , FEPlotData, "effective solute 2 concentration");
REGISTER_FEBIO_CLASS(FEPlotActualSol2Concentration      , FEPlotData, "solute 2 concentration"          );
REGISTER_FEBIO_CLASS(FEPlotSol2Flux                     , FEPlotData, "solute 2 flux"                   );
REGISTER_FEBIO_CLASS(FEPlotElectricPotential            , FEPlotData, "electric potential"              );
REGISTER_FEBIO_CLASS(FEPlotCurrentDensity               , FEPlotData, "current density"                 );
REGISTER_FEBIO_CLASS(FEPlotFixedChargeDensity           , FEPlotData, "fixed charge density"            );
REGISTER_FEBIO_CLASS(FEPlotNodalFluidFlux               , FEPlotData, "nodal fluid flux"                );
//-----------------------------------------------------------------------------

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
	if ((dynamic_cast<FEBiphasicSolidDomain*>(&dom))	|| 
		(dynamic_cast<FEBiphasicSoluteDomain*>(&dom)) ||
		(dynamic_cast<FETriphasicDomain*>(&dom)))
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
	if ((dynamic_cast<FEBiphasicSolidDomain*>(&dom)) || 
		(dynamic_cast<FEBiphasicSoluteDomain*>(&dom)) ||
		(dynamic_cast<FETriphasicDomain*>(&dom)))
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
	if ((dynamic_cast<FEBiphasicSolidDomain*>(&dom))
		|| (dynamic_cast<FEBiphasicSoluteDomain*>(&dom))
		|| (dynamic_cast<FETriphasicDomain*>(&dom)))
	{
		for (int i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);

			int nint = el.GaussPoints();
			int neln = el.Nodes();
			assert(nint == neln); // TODO: just for now

			// fluid flux at gauss points
			int j;
			vec3d vi[8];
			for (j=0; j<nint; ++j)
			{
				FEBiphasicMaterialPoint* pt = el.m_State[j]->ExtractData<FEBiphasicMaterialPoint>(); assert(pt);
				vi[j] = pt->m_w;
			}

			// project to nodes
			vec3d vn[8];
			matrix& Hi = el.m_pT->Hi;
			for (j=0; j<neln; ++j)
			{
				vn[j] = 0;
				for (int k=0; k<nint; ++k) 
				{
					vn[j].x += Hi[j][k]*vi[k].x;
					vn[j].y += Hi[j][k]*vi[k].y;
					vn[j].z += Hi[j][k]*vi[k].z;
				}
			}

			// output data
			for (j=0; j<neln; ++j)
			{
				a.push_back((float)vn[j].x);
				a.push_back((float)vn[j].y);
				a.push_back((float)vn[j].z);
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
bool FEPlotActualSol1Concentration::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESaltMaterialPoint* pt = (mp.ExtractData<FESaltMaterialPoint>());
				
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
bool FEPlotSol1Flux::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
	vec3d ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
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
bool FEPlotActualSol2Concentration::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
	if (ptd)
	{
		for (i=0; i<ptd->Elements(); ++i)
		{
			FESolidElement& el = ptd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESaltMaterialPoint* pt = (mp.ExtractData<FESaltMaterialPoint>());
				
				if (pt) ew += pt->m_ca[1];
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSol2Flux::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
	vec3d ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
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
				
				if (pt) ew += pt->m_j[1];
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotCurrentDensity::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	float af[3];
	vec3d ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
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
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFixedChargeDensity::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FETriphasicDomain* ptd = dynamic_cast<FETriphasicDomain*>(&dom);
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
bool FEPlotEffectiveSol1Concentration::Save(FEDomain &dom, vector<float>& a)
{
	FETriphasicDomain* pd = dynamic_cast<FETriphasicDomain*>(&dom);
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
bool FEPlotEffectiveSol2Concentration::Save(FEDomain &dom, vector<float>& a)
{
	FETriphasicDomain* pd = dynamic_cast<FETriphasicDomain*>(&dom);
	if (pd)
	{
		int N = pd->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pd->Node(i);
			a.push_back((float) node.m_ct[1]);
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
	if (pm == 0)
	{
		FENestedMaterial* pnm = dynamic_cast<FENestedMaterial*>(pmat);
		if (pnm)
		{
			pmat = pnm->m_pBase;
			FEElasticMixture* pm = dynamic_cast<FEElasticMixture*>(pmat);
			if (pm == 0) return false;
		}
		else return false;
	}

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
