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
			FERigidBody& rb = dynamic_cast<FERigidBody&>(*m_pfem->Object(node.m_rid));
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
//! Store the average elasticity for each element.
bool FEPlotElementElasticity::Save(FEDomain& dom, vector<float>& a)
{
	// write solid elasticity
	FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&dom);
	if (pbd) return WriteSolidElasticity(*pbd, a);
    
	FELinearSolidDomain* pbl = dynamic_cast<FELinearSolidDomain*>(&dom);
	if (pbl) return WriteLinearSolidElasticity(*pbl, a);
    
	// write shell elasticity
	FEElasticShellDomain* pbs = dynamic_cast<FEElasticShellDomain*>(&dom);
	if (pbs) return WriteShellElasticity(*pbs, a);
    
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElementElasticity::WriteSolidElasticity(FEElasticSolidDomain& d, vector<float>& a)
{
    FEMaterial* pm = dynamic_cast<FEMaterial*> (d.GetMaterial());
    FEElasticMaterial* pme = pm->GetElasticMaterial();
    if (pme == 0) return false;
    
    tens4ds c;
    
	// write solid element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FESolidElement& el = d.Element(i);
        
		float s[21] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
        
		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& pt = (*el.m_State[j]->ExtractData<FEMaterialPoint>());
            c = pme->Tangent(pt);
            
            for (int k=0; k<21; ++k) s[k] += (float) (f*c.d[k]);
		}
        
        for (int k=0; k<21; ++k) a.push_back(s[k]);
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotElementElasticity::WriteShellElasticity(FEElasticShellDomain& d, vector<float>& a)
{
    FEMaterial* pm = dynamic_cast<FEMaterial*> (d.GetMaterial());
    FEElasticMaterial* pme = pm->GetElasticMaterial();
    if (pme == 0) return false;
    
    tens4ds c;
    
	// write shell element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FEShellElement& el = d.Element(i);
        
		float s[21] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
        
		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& pt = (*el.m_State[j]->ExtractData<FEMaterialPoint>());
            c = pme->Tangent(pt);
            
            for (int k=0; k<21; ++k) s[k] += (float) (f*c.d[k]);
		}
        
        for (int k=0; k<21; ++k) a.push_back(s[k]);
	}
    
	return true;
}


//-----------------------------------------------------------------------------
bool FEPlotElementElasticity::WriteLinearSolidElasticity(FELinearSolidDomain& d, vector<float>& a)
{
    FEMaterial* pm = dynamic_cast<FEMaterial*> (d.GetMaterial());
    FEElasticMaterial* pme = pm->GetElasticMaterial();
    if (pme == 0) return false;
    
    tens4ds c;
    
	// write solid element data
	for (int i=0; i<d.Elements(); ++i)
	{
		FESolidElement& el = d.Element(i);
        
		float s[21] = {0};
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
        
		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& pt = (*el.m_State[j]->ExtractData<FEMaterialPoint>());
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
				vec3d ri;
				ri.x = pt.m_Q[0][0];
				ri.y = pt.m_Q[1][0];
				ri.z = pt.m_Q[2][0];

				r += pt.m_F*ri;
			}
			r /= (double) n;
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


//-----------------------------------------------------------------------------
bool FEPlotShellStrain::Save(FEDomain &dom, std::vector<float> &a)
{
	FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&dom);
	if (psd == 0) return false;
	int NE = psd->Elements();
	for (int i=0; i<NE; ++i)
	{
		FEShellElement& el = psd->Element(i);
		int ni = el.Nodes();
		mat3ds E; E.zero();
		for (int j=0; j<ni; ++j)
		{
			FEElasticMaterialPoint& ptm = *(el.m_State[j + ni]->ExtractData<FEElasticMaterialPoint>());
			FEElasticMaterialPoint& pti = *(el.m_State[j     ]->ExtractData<FEElasticMaterialPoint>());
			FEElasticMaterialPoint& pto = *(el.m_State[j+2*ni]->ExtractData<FEElasticMaterialPoint>());

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

	FESolidDomain& d = dynamic_cast<FESolidDomain&>(dom);
	int NE = d.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& e = d.Element(i);
		int nint = e.GaussPoints();
		double lam = 0.0;
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *e.m_State[j]->GetPointData(0);
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
	if (dynamic_cast<FESolidDomain*>(&dom) == 0) return false;

	// get the domain
	FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
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
				FEElasticMaterialPoint& ep = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
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
	if (dynamic_cast<FESolidDomain*>(&dom) == 0) return false;

	// get the domain
	FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
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
				FEElasticMaterialPoint& ep = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
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
	if (dynamic_cast<FESolidDomain*>(&dom) == 0) return false;

	// get the domain
	FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
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
				FEElasticMaterialPoint& ep = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
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
	if (dynamic_cast<FESolidDomain*>(&dom) == 0) return false;

	// get the domain
	FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);
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
				FEElasticMaterialPoint& ep = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
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
