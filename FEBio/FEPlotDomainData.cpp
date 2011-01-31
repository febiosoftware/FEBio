#include "stdafx.h"
#include "FEPlotDomainData.h"
#include "FEPlotDataFactory.h"

//-----------------------------------------------------------------------------
REGISTER_PLOTDATA(FEPlotEffectiveFluidPressure       , "effective fluid pressure"      );
REGISTER_PLOTDATA(FEPlotActualFluidPressure          , "fluid pressure"                );
REGISTER_PLOTDATA(FEPlotElementStress                , "stress"                        );
REGISTER_PLOTDATA(FEPlotRelativeVolume               , "relative volume"               );
REGISTER_PLOTDATA(FEPlotFluidFlux                    , "fluid flux"                    );
REGISTER_PLOTDATA(FEPlotFiberVector                  , "fiber vector"                  );
REGISTER_PLOTDATA(FEPlotEffectiveSoluteConcentration , "effective solute concentration");
REGISTER_PLOTDATA(FEPlotShellThickness               , "shell thickness"               );
REGISTER_PLOTDATA(FEPlotActualSoluteConcentration    , "solute concentration"          );
REGISTER_PLOTDATA(FEPlotSoluteFlux                   , "solute flux"                   );
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//! Store the average stresses for each element. 
bool FEPlotElementStress::Save(FEDomain& dom, vector<float>& a)
{
	// write solid stresses
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (pbd) return WriteSolidStress(*pbd, a);

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
	if ((dynamic_cast<FEPoroSolidDomain*>(&dom)) || (dynamic_cast<FEBiphasicDomain*>(&dom)))
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEPoroElasticMaterialPoint* pt = (mp.ExtractData<FEPoroElasticMaterialPoint>());
				
				if (pt) ew += pt->m_pa;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	else if (dynamic_cast<FEBiphasicSoluteDomain*>(&dom))
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average concentration
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutePoroElasticMaterialPoint* pt = (mp.ExtractData<FESolutePoroElasticMaterialPoint>());
				
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
	if ((dynamic_cast<FEPoroSolidDomain*>(&dom)) || (dynamic_cast<FEBiphasicDomain*>(&dom)))
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);

			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEPoroElasticMaterialPoint* pt = (mp.ExtractData<FEPoroElasticMaterialPoint>());

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
	else if (dynamic_cast<FEBiphasicSoluteDomain*>(&dom))
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average flux
			ew = vec3d(0,0,0);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FESolutePoroElasticMaterialPoint* pt = (mp.ExtractData<FESolutePoroElasticMaterialPoint>());
				
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
				FESolutePoroElasticMaterialPoint* pt = (mp.ExtractData<FESolutePoroElasticMaterialPoint>());
				
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
				FESolutePoroElasticMaterialPoint* pt = (mp.ExtractData<FESolutePoroElasticMaterialPoint>());
				
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
			dom.UnpackElement(e);
			int n = e.Nodes();
			vec3d* D = e.Dt();
			for (int j=0; j<n; ++j)
			{
				double h = e.m_h0[j] * D[j].norm();
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
	FEPoroSolidDomain* pe = dynamic_cast<FEPoroSolidDomain*>(&dom);
	FEBiphasicDomain* pd = dynamic_cast<FEBiphasicDomain*>(&dom);
	FEBiphasicSoluteDomain* psd = dynamic_cast<FEBiphasicSoluteDomain*>(&dom);
	if (pe)
	{
		int N = pe->Nodes();
		for (int i=0; i<N; ++i)
		{
			FENode& node = pe->Node(i);
			a.push_back((float) node.m_pt);
		}
		return true;
	}
	else if (pd)
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
			a.push_back((float) node.m_ct);
		}
		return true;
	}
	return false;
}
