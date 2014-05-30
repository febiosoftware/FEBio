#include "stdafx.h"
#include "FEBioMechPlot.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FERemodelingElasticMaterial.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEElasticMixture.h"
#include "FEUT4Domain.h"
#include "FEPreStrainTransIsoMR.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEContactSurface.h"
#include "FECore/FERigidBody.h"
#include "FESPRProjection.h"
#include "FEUncoupledElasticMixture.h"
#include "FERigidMaterial.h"

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

//-----------------------------------------------------------------------------
bool FEPlotRigidReactionTorque::Save(FEMesh& m, vector<float>& a)
{
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = m.Node(i);
		if (node.m_rid >= 0)
		{
			FERigidBody& rb = static_cast<FERigidBody&>(*m_pfem->Object(node.m_rid));
			a.push_back((float)rb.m_Mr.x);
			a.push_back((float)rb.m_Mr.y);
			a.push_back((float)rb.m_Mr.z);
		}
		else
		{
			a.push_back(0.f);
			a.push_back(0.f);
			a.push_back(0.f);
		}
	}
	return true;
}

//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
// Plot contact gap
bool FEPlotContactGap::Save(FESurface& surf, vector<float>& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	double gn[MFN];
	a.assign(MFN*NF, 0.f);
	for (int i=0; i<NF; ++i) 
	{
		FESurfaceElement& f = pcs->Element(i);
		pcs->GetNodalContactGap(i, gn);
		int ne = f.m_lnode.size();
		for (int j = 0; j< ne; ++j) a[MFN*i + j] = (float) gn[j];
	}
	return true;
}

//-----------------------------------------------------------------------------
// Plot contact pressure
bool FEPlotContactPressure::Save(FESurface &surf, vector<float>& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	double tn[MFN];
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = pcs->Element(i);
		pcs->GetNodalContactPressure(i, tn);
		int ne = el.Nodes();
		for (int k=0; k<ne; ++k) a[MFN*i + k] = (float) tn[k];
	}
	return true;
}

//-----------------------------------------------------------------------------
// Plot contact traction
bool FEPlotContactTraction::Save(FESurface &surf, std::vector<float> &a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;

	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(3*MFN*NF, 0.f);
	vec3d tn[MFN];
	for (int j=0; j<NF; ++j)
	{
		FESurfaceElement& el = pcs->Element(j);
		pcs->GetNodalContactTraction(j, tn);

		// store in archive
		int ne = el.Nodes();
		for (int k=0; k<ne; ++k)
		{
			a[3*MFN*j +3*k   ] = (float) tn[k].x;
			a[3*MFN*j +3*k +1] = (float) tn[k].y;
			a[3*MFN*j +3*k +2] = (float) tn[k].z;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotContactForce::Save(FESurface &surf, std::vector<float> &a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(3*MFN*NF, 0.f);
	vec3d fn = pcs->GetContactForce();
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

//-----------------------------------------------------------------------------
// Plot contact area
bool FEPlotContactArea::Save(FESurface &surf, vector<float>& a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	int NF = pcs->Elements();
	const int MFN = FEBioPlotFile::PLT_MAX_FACET_NODES;
	a.assign(MFN*NF, 0.f);
	double area;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = pcs->Element(i);
		area = pcs->GetContactArea();
		int ne = el.Nodes();
		for (int k=0; k<ne; ++k) a[MFN*i + k] = (float) area;
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
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		float s[6] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
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
//! Store the average elasticity for each element.
bool FEPlotElementElasticity::Save(FEDomain& dom, vector<float>& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
    if ((pme == 0) || pme->IsRigid()) return false;
    
    tens4ds c;
    
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
        
		float s[21] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
        
		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& pt = *el.GetMaterialPoint(j);
            c = pme->Tangent(pt);
            
            for (int k=0; k<21; ++k) s[k] += (float) (f*c.d[k]);
		}
        
        for (int k=0; k<21; ++k) a.push_back(s[k]);
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotStrainEnergyDensity::Save(FEDomain &dom, vector<float>& a)
{
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		int N = bd.Elements();
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average strain energy
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FERemodelingMaterialPoint* pt = (mp.ExtractData<FERemodelingMaterialPoint>());
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
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		int N = bd.Elements();
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average strain energy
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FERemodelingMaterialPoint* rpt = (mp.ExtractData<FERemodelingMaterialPoint>());
				
				if (rpt) ew += rpt->m_sed/rpt->m_rhor;
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
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		int N = bd.Elements();
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average mass density
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FERemodelingMaterialPoint* pt = (mp.ExtractData<FERemodelingMaterialPoint>());
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
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		int N = bd.Elements();
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average flux
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
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
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if (pme == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		int BE = bd.Elements();
		for (int i=0; i<BE; ++i)
		{
			FESolidElement& el = bd.Element(i);
			int n = el.GaussPoints();
			vec3d r = vec3d(0,0,0);
			for (int j=0; j<n; ++j)
			{
				FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
				vec3d ri;
				ri.x = pt.m_Q[0][0];
				ri.y = pt.m_Q[1][0];
				ri.z = pt.m_Q[2][0];

				r += pt.m_F*ri;
			}
			r /= (double) n;

			float f[3];
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
	if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
		int NS = sd.Elements();
		for (int i=0; i<NS; ++i)
		{
			FEShellElement& e = sd.Element(i);
			int n = e.Nodes();
			for (int j=0; j<n; ++j)
			{
				vec3d D = sd.GetMesh()->Node(e.m_node[j]).m_Dt;
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
	if (m.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& d = static_cast<FESolidDomain&>(m);
		int N = d.Elements();
		for (int i=0; i<N; ++i)
		{
			FESolidElement& el = d.Element(i);

			float D = 0.f;
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEDamageMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEDamageMaterialPoint>());
				if (ppt)
				{
					FEDamageMaterialPoint& pt = *ppt;
					D += (float) pt.m_D;
				}

				FETIMRDamageMaterialPoint* pt2 = (el.GetMaterialPoint(j)->ExtractData<FETIMRDamageMaterialPoint>());
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
			FEElasticMixtureMaterialPoint& pt = *e.GetMaterialPoint(n)->ExtractData<FEElasticMixtureMaterialPoint>();
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


//-----------------------------------------------------------------------------
bool FEPlotShellStrain::Save(FEDomain &dom, std::vector<float> &a)
{
	if (dom.Class() != FE_DOMAIN_SHELL) return false;

	FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
	int NE = sd.Elements();
	for (int i=0; i<NE; ++i)
	{
		FEShellElement& el = sd.Element(i);
		int ni = el.Nodes();
		mat3ds E; E.zero();
		for (int j=0; j<ni; ++j)
		{
			FEElasticMaterialPoint& ptm = *(el.GetMaterialPoint(j + ni)->ExtractData<FEElasticMaterialPoint>());
			FEElasticMaterialPoint& pti = *(el.GetMaterialPoint(j     )->ExtractData<FEElasticMaterialPoint>());
			FEElasticMaterialPoint& pto = *(el.GetMaterialPoint(j+2*ni)->ExtractData<FEElasticMaterialPoint>());

			E += ptm.Strain();
			E += pto.Strain();
			E += pti.Strain();
		}
		E /= (3.0*ni);

		a.push_back((float) E.xx());
		a.push_back((float) E.yy());
		a.push_back((float) E.zz());
		a.push_back((float) E.xy());
		a.push_back((float) E.yz());
		a.push_back((float) E.xz());
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFiberPreStretch::Save(FEDomain& dom, vector<float>& a)
{
	FEPreStrainTransIsoMR* pm = dynamic_cast<FEPreStrainTransIsoMR*>(dom.GetMaterial());
	if (pm == 0) 
	{
		// see if we have a mixture
		FEUncoupledElasticMixture* pmix = dynamic_cast<FEUncoupledElasticMixture*>(dom.GetMaterial());
		if (pmix == 0) return false;

		// For now, we just grab the first match we find
		int N = pmix->Materials();
		for (int i=0; i<N; ++i)
		{
			pm = dynamic_cast<FEPreStrainTransIsoMR*>(pmix->GetMaterial(i));
			if (pm) break;
		}
		if (pm == 0) return false;
	}

	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int nint = e.GaussPoints();
		double lam = 0.0;
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(j)->GetPointData(0);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEFiberPreStretchMaterialPoint& psp = *mp.ExtractData<FEFiberPreStretchMaterialPoint>();
			mat3d& F = pt.m_F;

			vec3d a0(pt.m_Q[0][0], pt.m_Q[1][0], pt.m_Q[2][0]);
			vec3d a = F*a0;
			double lRtor = a.norm();

			lam += psp.m_lam*lRtor;
//			lam += psp.m_lam;
		}
		lam /= (double) nint;
		a.push_back((float)lam);
	}
	return true;
}


//-----------------------------------------------------------------------------
bool FEPlotSPRStresses::Save(FEDomain& dom, vector<float>& a)
{
	const int LUT[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};

	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	int NN = sd.Nodes();
	int NE = sd.Elements();

	// build the element data array
	vector< vector<double> > ED;
	ED.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		ED[i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	vector<double> val[6];

	// loop over stress components
	for (int n=0; n<6; ++n)
	{
		// fill the ED array
		for (int i=0; i<NE; ++i)
		{
			FESolidElement& el = sd.Element(i);
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEElasticMaterialPoint& ep = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
				mat3ds& s = ep.m_s;
				ED[i][j] = s(LUT[n][0], LUT[n][1]);
			}
		}

		// project to nodes
		map.Project(sd, ED, val[n]);
	}

	// copy results to archive
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float)val[0][i]);
		a.push_back((float)val[1][i]);
		a.push_back((float)val[2][i]);
		a.push_back((float)val[3][i]);
		a.push_back((float)val[4][i]);
		a.push_back((float)val[5][i]);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRPrincStresses::Save(FEDomain& dom, vector<float>& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	int NN = sd.Nodes();
	int NE = sd.Elements();

	// build the element data array
	vector< vector<double> > ED;
	ED.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		ED[i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	vector<double> val[3];

	// loop over stress components
	for (int n=0; n<3; ++n)
	{
		// fill the ED array
		for (int i=0; i<NE; ++i)
		{
			FESolidElement& el = sd.Element(i);
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEElasticMaterialPoint& ep = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
				mat3ds& s = ep.m_s;
				double l[3];
				s.exact_eigen(l);
				ED[i][j] = l[n];
			}
		}

		// project to nodes
		map.Project(sd, ED, val[n]);
	}

	// copy results to archive
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float)val[0][i]);
		a.push_back((float)val[1][i]);
		a.push_back((float)val[2][i]);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRTestLinear::Save(FEDomain& dom, vector<float>& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	int NN = sd.Nodes();
	int NE = sd.Elements();

	// build the element data array
	vector< vector<double> > ED;
	ED.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		ED[i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	vector<double> val[3];

	// loop over stress components
	for (int n=0; n<3; ++n)
	{
		// fill the ED array
		for (int i=0; i<NE; ++i)
		{
			FESolidElement& el = sd.Element(i);
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEElasticMaterialPoint& ep = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
				vec3d r = ep.m_rt;
				double l[3] = {r.x, r.y, r.z};
				ED[i][j] = l[n];
			}
		}

		// project to nodes
		map.Project(sd, ED, val[n]);
	}

	// copy results to archive
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float)val[0][i]);
		a.push_back((float)val[1][i]);
		a.push_back((float)val[2][i]);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRTestQuadratic::Save(FEDomain& dom, vector<float>& a)
{
	// For now, this is only available for solid domains
	if (dom.Class() != FE_DOMAIN_SOLID) return false;

	// get the domain
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	int NN = sd.Nodes();
	int NE = sd.Elements();

	// build the element data array
	vector< vector<double> > ED;
	ED.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		ED[i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	vector<double> val[6];

	// loop over stress components
	for (int n=0; n<6; ++n)
	{
		// fill the ED array
		for (int i=0; i<NE; ++i)
		{
			FESolidElement& el = sd.Element(i);
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEElasticMaterialPoint& ep = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
				vec3d r = ep.m_rt;
				double l[6] = {r.x*r.x, r.y*r.y, r.z*r.z, r.x*r.y, r.y*r.z, r.x*r.z};
				ED[i][j] = l[n];
			}
		}

		// project to nodes
		map.Project(sd, ED, val[n]);
	}

	// copy results to archive
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float)val[0][i]);
		a.push_back((float)val[1][i]);
		a.push_back((float)val[2][i]);
		a.push_back((float)val[3][i]);
		a.push_back((float)val[4][i]);
		a.push_back((float)val[5][i]);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidPotentialEnergy::Save(FEDomain& dom, vector<float>& a)
{
	// only store this on rigid bodies
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;

	// we evaluate the potential energy of all rigid bodies
	int nrb = m_pfem->Objects();
	double Ep = 0.0;
	for (int i=0; i<nrb; ++i)
	{
		// get the corresponding rigid body
		FERigidBody& rb = static_cast<FERigidBody&>(*m_pfem->Object(i));

		// evaluate potential energy
		// NOTE: we assume g = 1 and body force applied in z-direction
		Ep += rb.m_mass*rb.m_rt.z;
	}

	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i) a.push_back((float) Ep);

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidKineticEnergy::Save(FEDomain& dom, vector<float>& a)
{
	// only store this on rigid bodies
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;

	// we evaluate the potential energy of all rigid bodies
	int nrb = m_pfem->Objects();
	double K = 0.0;
	for (int i=0; i<nrb; ++i)
	{
		// get the corresponding rigid body
		FERigidBody& rb = static_cast<FERigidBody&>(*m_pfem->Object(i));

		// get translational kinetic energy
		double KT = 0.5*rb.m_mass*(rb.m_vt*rb.m_vt);

		// get rotational matrix
		mat3d Q = rb.m_qt.RotationMatrix();

		// calculate moment of inertia
		mat3d It = Q*rb.m_moi*Q.transpose();

		// get rotational energy
		double KR = 0.5*(rb.m_wt*(It*rb.m_wt));

		// total kinetic energy
		K += KT + KR;
	}

	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i) a.push_back((float) K);

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidTotalEnergy::Save(FEDomain& dom, vector<float>& a)
{
	// only store this on rigid bodies
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;

	// we evaluate the potential energy of all rigid bodies
	int nrb = m_pfem->Objects();
	double E = 0.0;
	for (int i=0; i<nrb; ++i)
	{
		// get the corresponding rigid body
		FERigidBody& rb = static_cast<FERigidBody&>(*m_pfem->Object(i));

		// get translational kinetic energy
		double KT = 0.5*rb.m_mass*(rb.m_vt*rb.m_vt);

		// get rotational matrix
		mat3d Q = rb.m_qt.RotationMatrix();

		// calculate moment of inertia
		mat3d It = Q*rb.m_moi*Q.transpose();

		// get rotational energy
		double KR = 0.5*(rb.m_wt*(It*rb.m_wt));

		// total kinetic energy
		double K = KT + KR;

		// evaluate potential energy
		// NOTE: we assume g = 1 and body force applied in z-direction
		double H = rb.m_mass*rb.m_rt.z;

		// total energy
		E += K + H;
	}

	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i) a.push_back((float) E);

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidVelocity::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);

	// get the rigid body
	FERigidBody& rb = static_cast<FERigidBody&>(*m_pfem->Object(prm->GetRigidBodyID()));

	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) rb.m_vt.x);
		a.push_back((float) rb.m_vt.y);
		a.push_back((float) rb.m_vt.z);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularVelocity::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);

	// get the rigid body
	FERigidBody& rb = static_cast<FERigidBody&>(*m_pfem->Object(prm->GetRigidBodyID()));

	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) rb.m_wt.x);
		a.push_back((float) rb.m_wt.y);
		a.push_back((float) rb.m_wt.z);
	}

	return true;
}
