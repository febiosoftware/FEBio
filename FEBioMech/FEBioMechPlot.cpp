#include "stdafx.h"
#include "FEBioMechPlot.h"
#include "FEBioMech/FEDamageNeoHookean.h"
#include "FEBioMech/FEDamageTransIsoMooneyRivlin.h"
#include "FEBioMech/FERigidSolidDomain.h"
#include "FEBioMech/FERigidShellDomain.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEUT4Domain.h"

//=============================================================================
//                            N O D E   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotNodeDisplacement::Save(FEMesh& m, vector<float>& a)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) (node.m_rt.x - node.m_r0.x);
		xf[1] = (float) (node.m_rt.y - node.m_r0.y);
		xf[2] = (float) (node.m_rt.z - node.m_r0.z);

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeVelocity::Save(FEMesh& m, vector<float>& a)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_vt.x;
		xf[1] = (float) node.m_vt.y;
		xf[2] = (float) node.m_vt.z;

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeAcceleration::Save(FEMesh& m, vector<float>& a)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_at.x;
		xf[1] = (float) node.m_at.y;
		xf[2] = (float) node.m_at.z;

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Store nodal reaction forces
bool FEPlotNodeReactionForces::Save(FEMesh& m, vector<float>& a)
{
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = m.Node(i);
		a.push_back((float) node.m_Fr.x);
		a.push_back((float) node.m_Fr.y);
		a.push_back((float) node.m_Fr.z);
	}
	return true;
}

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

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
				s[0] += (float) (f*pt.m_s.xx());
				s[1] += (float) (f*pt.m_s.yy());
				s[2] += (float) (f*pt.m_s.zz());
				s[3] += (float) (f*pt.m_s.xy());
				s[4] += (float) (f*pt.m_s.yz());
				s[5] += (float) (f*pt.m_s.xz());
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
				s[0] += (float) (f*pt.m_s.xx());
				s[1] += (float) (f*pt.m_s.yy());
				s[2] += (float) (f*pt.m_s.zz());
				s[3] += (float) (f*pt.m_s.xy());
				s[4] += (float) (f*pt.m_s.yz());
				s[5] += (float) (f*pt.m_s.xz());
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
				mat3ds& es = pt.m_s;
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
bool FEPlotStrainEnergyDensity::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average strain energy
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEElasticMaterialPoint* pt = (mp.ExtractData<FEElasticMaterialPoint>());
				
				if (pt) ew += pt->m_sed;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSpecificStrainEnergy::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average strain energy
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEElasticMaterialPoint* pt = (mp.ExtractData<FEElasticMaterialPoint>());
				
				if (pt) ew += pt->m_sed/pt->m_rhor;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotDensity::Save(FEDomain &dom, vector<float>& a)
{
	int i, j;
	double ew;
	FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&dom);
	if (pbd)
	{
		for (i=0; i<pbd->Elements(); ++i)
		{
			FESolidElement& el = pbd->Element(i);
			
			// calculate average mass density
			ew = 0;
			for (j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.m_State[j];
				FEElasticMaterialPoint* pt = (mp.ExtractData<FEElasticMaterialPoint>());
				
				if (pt) ew += pt->m_rhor;
			}
			
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
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
				
				if (pt) ew += pt->m_J;
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
				r.x += pt.m_Q[0][0];
				r.y += pt.m_Q[1][0];
				r.z += pt.m_Q[2][0];
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
