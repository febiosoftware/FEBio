#include "stdafx.h"
#include "FEBioMechPlot.h"
#include "FEDamageNeoHookean.h"
#include "FEDamageTransIsoMooneyRivlin.h"
#include "FERemodelingElasticMaterial.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FEUT4Domain.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEContactSurface.h"
#include "FECore/FERigidBody.h"
#include "FESPRProjection.h"
#include "FEUncoupledElasticMixture.h"
#include "FERigidMaterial.h"
#include "FEVolumeConstraint.h"
#include "FEMicroMaterial.h"
#include "FEMicroMaterial2O.h"
#include "FEFacet2FacetSliding.h"

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
bool FEPlotRigidReactionForce::Save(FEMesh& m, vector<float>& a)
{
    int N = m.Nodes();
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
    for (int i=0; i<N; ++i)
    {
        FENode& node = m.Node(i);
        if (node.m_rid >= 0)
        {
            FERigidBody& rb = *rigid.Object(node.m_rid);
            a.push_back((float)rb.m_Fr.x);
            a.push_back((float)rb.m_Fr.y);
            a.push_back((float)rb.m_Fr.z);
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

//-----------------------------------------------------------------------------
bool FEPlotRigidReactionTorque::Save(FEMesh& m, vector<float>& a)
{
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = m.Node(i);
		if (node.m_rid >= 0)
		{
			FERigidBody& rb = *rigid.Object(node.m_rid);
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
		int ne = f.Nodes();
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

//-----------------------------------------------------------------------------
// Plot contact penalty parameter
bool FEPlotContactPenalty::Save(FESurface& surf, vector<float>& a)
{
	FEFacetSlidingSurface* ps = dynamic_cast<FEFacetSlidingSurface*>(&surf);
	if (ps)
	{
		int NF = ps->Elements();
		for (int i=0; i<NF; ++i)
		{
			FESurfaceElement& el = ps->Element(i);
			int ni = el.GaussPoints();
			double p = 0.0;
			for (int n=0; n<ni; ++n)
			{
				FEFacetSlidingSurface::Data& pt = ps->m_Data[i][n];
				p += pt.m_eps;
			}
			if (ni > 0) p /= (double) ni;

			a.push_back((float) p);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotVolumePressure::Save(FESurface& S, vector<float>& a)
{
	FEVolumeSurface* pvs = dynamic_cast<FEVolumeSurface*>(&S);
	if (pvs == 0) return false;

	FEVolumeSurface& s = *pvs;
	double p = s.m_p;

	int NE = s.Elements();
	for (int i=0; i<NE; ++i) a.push_back((float)p);

	return true;
}

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Store the average deformation Hessian (G) for each element. 
bool FEPlotElementGnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3drs Gavg; Gavg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		Gavg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				Gavg += (pt2O.m_G)*f;
				
			}
		}

		L2_norm = (float) sqrt(Gavg.tripledot3rs(Gavg));

		a.push_back(L2_norm);
	}
	
	return true;
}

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

				pt.m_F_prev = pt.m_F;

				FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());

				if (mmppt)
				{
					FEMicroMaterialPoint& mmpt = *mmppt;
					mmpt.m_macro_energy += mmpt.m_macro_energy_inc;
					mmpt.m_micro_energy += mmpt.m_micro_energy_inc;
					mmpt.m_energy_diff = fabs(mmpt.m_macro_energy - mmpt.m_micro_energy); 
				}
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
//! Store the norm of the average Cauchy stress for each element. 
bool FEPlotElementsnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	mat3ds s_avg; s_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		s_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
			if (ppt)
			{
				FEElasticMaterialPoint& pt = *ppt;
				s_avg += (pt.m_s)*f;
			}
		}

		L2_norm = (float) sqrt(s_avg.dotdot(s_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average Cauchy stress moment for each element. 
bool FEPlotElementtaunorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3ds tau_avg; tau_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		tau_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				tau_avg += (pt2O.m_tau)*f;
				
				FEElasticMaterialPoint& pt = *ppt;
				double norm = pt2O.m_G.tripledot3rs(pt2O.m_G);

				pt2O.m_G_prev = pt2O.m_G;

				pt2O.m_macro_energy += pt2O.m_macro_energy_inc;
				pt2O.m_micro_energy += pt2O.m_micro_energy_inc;
				pt2O.m_energy_diff = fabs(pt2O.m_macro_energy - pt2O.m_micro_energy); 
			}
		}

		L2_norm = (float) sqrt(tau_avg.tripledot3s(tau_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK1 stress for each element.
bool FEPlotElementPK1norm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	mat3d PK1_avg; PK1_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		PK1_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			
			if (mmppt)
			{
				FEMicroMaterialPoint& mmpt = *mmppt;
				PK1_avg += (mmpt.m_PK1)*f;	
			}
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
				if (mmppt2O)
				{
					FEMicroMaterialPoint2O& mmpt2O = *mmppt2O;
					PK1_avg += (mmpt2O.m_PK1)*f;	
				}
			}
		}

		L2_norm = (float) sqrt(PK1_avg.dotdot(PK1_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK1 stress moment for each element. 
bool FEPlotElementQK1norm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3drs QK1_avg; QK1_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		QK1_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				QK1_avg += (pt2O.m_QK1)*f;
				
			}
		}

		L2_norm = (float) sqrt(QK1_avg.tripledot3rs(QK1_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK2 stress for each element.
bool FEPlotElementSnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	mat3ds S_avg; S_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		S_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			
			if (mmppt)
			{
				FEMicroMaterialPoint& mmpt = *mmppt;
				S_avg += (mmpt.m_S)*f;	
			}
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
				if (mmppt2O)
				{
					FEMicroMaterialPoint2O& mmpt2O = *mmppt2O;
					S_avg += (mmpt2O.m_S)*f;	
				}
			}
		}

		L2_norm = (float) sqrt(S_avg.dotdot(S_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK2 stress moment for each element. 
bool FEPlotElementTnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3ds T_avg; T_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		T_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				T_avg += (pt2O.m_T)*f;
				
			}
		}

		L2_norm = (float) sqrt(T_avg.tripledot3s(T_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average infinitesimal strain gradient for each element. 
bool FEPlotElementinfstrnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3ds inf_strain_avg; inf_strain_avg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		inf_strain_avg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				inf_strain_avg += (pt2O.m_inf_str_grad)*f;
				
			}
		}

		L2_norm = (float) sqrt(inf_strain_avg.tripledot3s(inf_strain_avg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average Green-Lagrange strain gradient for each element. 
bool FEPlotElementGLstrnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3ds Havg; Havg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		Havg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				Havg += (pt2O.m_H)*f;
				
			}
		}

		L2_norm = (float) sqrt(Havg.tripledot3s(Havg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average Euler-Almansi strain gradient for each element. 
bool FEPlotElementEAstrnorm::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float L2_norm; L2_norm = 0.;
	tens3ds havg; havg.zero();

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
		havg.zero();

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
			
			if (ppt2O)
			{
				FEMicroMaterialPoint2O& pt2O = *ppt2O;
				havg += (pt2O.m_h)*f;
				
			}
		}

		L2_norm = (float) sqrt(havg.tripledot3s(havg));

		a.push_back(L2_norm);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Element macro-micro energy difference
bool FEPlotElementenergydiff::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float energy_diff = 0.;

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			
			if (mmppt)
			{
				FEMicroMaterialPoint& mmpt = *mmppt;
				energy_diff += (mmpt.m_energy_diff)*f;	
			}
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
				if (mmppt2O)
				{
					FEMicroMaterialPoint2O& mmpt2O = *mmppt2O;
					energy_diff += (mmpt2O.m_energy_diff)*f;	
				}
			}
		}

		a.push_back(energy_diff);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Element macro energy
bool FEPlotElementMacroEnergy::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float macro_energy = 0.;

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			
			if (mmppt)
			{
				FEMicroMaterialPoint& mmpt = *mmppt;
				macro_energy += (mmpt.m_macro_energy)*f;	
			}
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
				if (mmppt2O)
				{
					FEMicroMaterialPoint2O& mmpt2O = *mmppt2O;
					macro_energy += (mmpt2O.m_macro_energy)*f;	
				}
			}
		}

		a.push_back(macro_energy);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Element macro energy
bool FEPlotElementMicroEnergy::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	float micro_energy = 0.;

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			
			if (mmppt)
			{
				FEMicroMaterialPoint& mmpt = *mmppt;
				micro_energy += (mmpt.m_micro_energy)*f;	
			}
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
				if (mmppt2O)
				{
					FEMicroMaterialPoint2O& mmpt2O = *mmppt2O;
					micro_energy += (mmpt2O.m_micro_energy)*f;	
				}
			}
		}

		a.push_back(micro_energy);
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
    FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
    if ((pme == 0) || pme->IsRigid()) return false;
    
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average strain energy
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double sed = pme->StrainEnergyDensity(mp);
                ew += sed;
			}
			ew /= el.GaussPoints();
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotDevStrainEnergyDensity::Save(FEDomain &dom, vector<float>& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
    FEUncoupledMaterial* pmu = dynamic_cast<FEUncoupledMaterial*>(pme);
    if ((pme == 0) || pme->IsRigid() || (pmu == 0)) return false;
    
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average strain energy
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                double sed = pmu->DevStrainEnergyDensity(mp);
                ew += sed;
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
//			r /= (double) n;
			r.unit();

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
bool FEPlotFiberStretch::Save(FEDomain &dom, vector<float>& a)
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
			double l = 0.0;
			for (int j=0; j<n; ++j)
			{
				FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
				vec3d ri;
				ri.x = pt.m_Q[0][0];
				ri.y = pt.m_Q[1][0];
				ri.z = pt.m_Q[2][0];

				vec3d r = pt.m_F*ri;

				l += r.norm();
			}
			l /= (double) n;
			a.push_back((float) l);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotDevFiberStretch::Save(FEDomain &dom, vector<float>& a)
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
			double lamd = 0.0;
			for (int j=0; j<n; ++j)
			{
				FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();

				// get the deformation gradient
				mat3d& F = pt.m_F;
				double J = pt.m_J;
				double Jm13 = pow(J, -1.0/3.0);

				// get the material fiber axis
				vec3d ri;
				ri.x = pt.m_Q[0][0];
				ri.y = pt.m_Q[1][0];
				ri.z = pt.m_Q[2][0];

				// apply deformation
				vec3d r = pt.m_F*ri;

				// calculate the deviatoric fiber stretch
				double lam = r.norm();
				lamd += lam*Jm13;
			}
			lamd /= (double) n;
			a.push_back((float) lamd);
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
bool FEPlotDamage::Save(FEDomain &dom, vector<float>& a)
{
	int N = dom.Elements();
	FEElasticMaterial* pmat = dom.GetMaterial()->GetElasticMaterial();
	if (dynamic_cast<FEElasticMixture*>(pmat)||dynamic_cast<FEUncoupledElasticMixture*>(pmat))
	{
		int NC = pmat->Properties();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);

			float D = 0.f;
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEElasticMixtureMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMixtureMaterialPoint>();
				for (int k=0; k<NC; ++k)
				{
					FEDamageMaterialPoint* ppd = pt.GetPointData(k)->ExtractData<FEDamageMaterialPoint>();
					if (ppd) D += (float) ppd->m_D;
				}
			}
			D /= (float) nint;
			a.push_back(D);
		}
	}
    else if (dynamic_cast<FEElasticMultigeneration*>(pmat))
    {
        FEElasticMultigeneration* pmg = dynamic_cast<FEElasticMultigeneration*>(pmat);
        int NC = pmg->Properties();
        for (int i=0; i<N; ++i)
        {
            FEElement& el = dom.ElementRef(i);
            
            float D = 0.f;
            int nint = el.GaussPoints();
            for (int j=0; j<nint; ++j)
            {
                FEMultigenerationMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEMultigenerationMaterialPoint>();
                for (int k=0; k<NC; ++k)
                {
                    FEDamageMaterialPoint* ppd = pt.GetPointData(k)->ExtractData<FEDamageMaterialPoint>();
                    FEElasticMixtureMaterialPoint* pem = pt.GetPointData(k)->ExtractData<FEElasticMixtureMaterialPoint>();
                    if (ppd) D += (float) ppd->m_D;
                    else if (pem)
                    {
                        int NE = (int)pem->m_w.size();
                        for (int l=0; l<NE; ++l)
                        {
                            FEDamageMaterialPoint* ppd = pem->GetPointData(l)->ExtractData<FEDamageMaterialPoint>();
                            if (ppd) D += (float) ppd->m_D;
                        }
                    }
                }
            }
            D /= (float) nint;
            a.push_back(D);
        }
    }
	else
	{
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);

			float D = 0.f;
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j)
			{
				FEMaterialPoint& pt = *el.GetMaterialPoint(j);
				FEDamageMaterialPoint* ppd = pt.ExtractData<FEDamageMaterialPoint>();
				if (ppd) D += (float) ppd->m_D;
			}
			D /= (float) nint;
			a.push_back(D);
		}
	}
	return true;
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
bool FEPlotRigidDisplacement::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) rb.m_rt.x);
		a.push_back((float) rb.m_rt.y);
		a.push_back((float) rb.m_rt.z);
	}
    
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
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
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
bool FEPlotRigidAcceleration::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) rb.m_at.x);
		a.push_back((float) rb.m_at.y);
		a.push_back((float) rb.m_at.z);
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidRotation::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    vec3d q = rb.m_qt.GetVector()*rb.m_qt.GetAngle();
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) q.x);
		a.push_back((float) q.y);
		a.push_back((float) q.z);
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
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
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

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularAcceleration::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) rb.m_alt.x);
		a.push_back((float) rb.m_alt.y);
		a.push_back((float) rb.m_alt.z);
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidKineticEnergy::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    vec3d v = rb.m_vt;
    double m = rb.m_mass;
    vec3d w = rb.m_wt;
    mat3d Rt = rb.m_qt.RotationMatrix();
    mat3ds Jt = (Rt*rb.m_moi*Rt.transpose()).sym();
    double ke = ((v*v)*m + w*(Jt*w))/2;
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i) a.push_back((float) ke);
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidEuler::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());

	// get the Euler angles
	double E[3];
	quat2euler(rb.m_qt, E);
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) E[0]);
		a.push_back((float) E[1]);
		a.push_back((float) E[2]);
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidRotationVector::Save(FEDomain& dom, vector<float>& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());

	// get the rotation vector and angle
	double w = rb.m_qt.GetAngle();
	vec3d r = rb.m_qt.GetVector()*w;
    
	// copy results to archive
	int NN = dom.Nodes();
	for (int i=0; i<NN; ++i)
	{
		a.push_back((float) r.x);
		a.push_back((float) r.y);
		a.push_back((float) r.z);
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodalStresses::Save(FEDomain& dom, vector<float>& a)
{
	// make sure this is a solid-domain class
	FESolidDomain* pd = dynamic_cast<FESolidDomain*>(&dom);
	if (pd == 0) return false;

	// stress component look-up table
	int LUT[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};

	// temp storage 
	mat3ds s[FEElement::MAX_NODES];
	double si[27];	// 27 = max nr of integration points for now.
	double sn[FEElement::MAX_NODES];

	// loop over all elements
	int NE = pd->Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& e = pd->Element(i);
		int ne = e.Nodes();
		int ni = e.GaussPoints();

		// loop over stress-components
		for (int j=0; j<6; ++j)
		{
			// get the integration point values
			int j0 = LUT[j][0];
			int j1 = LUT[j][1];
			for (int k=0; k<ni; ++k) 
			{
				FEMaterialPoint& mp = *e.GetMaterialPoint(k);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				si[k] = pt.m_s(j0, j1);
			}

			// project to nodes
			e.project_to_nodes(si, sn);

			// store stress component
			for (int k=0; k<ne; ++k) s[k](j0, j1) = sn[k];
		}

		// push data to archive
		for (int j=0; j<ne; ++j)
		{
			a.push_back((float)s[j].xx());
			a.push_back((float)s[j].yy());
			a.push_back((float)s[j].zz());
			a.push_back((float)s[j].xy());
			a.push_back((float)s[j].yz());
			a.push_back((float)s[j].xz());
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average Euler-lagrange strain
bool FEPlotLagrangeStrain::Save(FEDomain& dom, vector<float>& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;

	// write solid element data
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		float s[6] = { 0 };
		int nint = el.GaussPoints();
		double f = 1.0 / (double)nint;
		mat3dd I(1.0); // identity tensor

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (int j = 0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
			if (ppt)
			{
				mat3d C = ppt->RightCauchyGreen();
				mat3ds E = ((C - I)*0.5).sym();
				s[0] += (float)(f*E.xx());
				s[1] += (float)(f*E.yy());
				s[2] += (float)(f*E.zz());
				s[3] += (float)(f*E.xy());
				s[4] += (float)(f*E.yz());
				s[5] += (float)(f*E.xz());
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
