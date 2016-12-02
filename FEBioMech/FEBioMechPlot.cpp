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
#include "FEMortarSlidingContact.h"
#include <FECore/FERigidSystem.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotNodeVelocity::Save(FEMesh& m, FEDataStream& a)
{
	const int dof_VX = GetFEModel()->GetDOFIndex("vx");
	const int dof_VY = GetFEModel()->GetDOFIndex("vy");
	const int dof_VZ = GetFEModel()->GetDOFIndex("vz");
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a << node.get_vec3d(dof_VX, dof_VY, dof_VZ);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeAcceleration::Save(FEMesh& m, FEDataStream& a)
{
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a << node.m_at;
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Store nodal reaction forces
bool FEPlotNodeReactionForces::Save(FEMesh& m, FEDataStream& a)
{
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = m.Node(i);
		a << node.m_Fr;
	}
	return true;
}

//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
// Plot contact gap
bool FEPlotContactGap::Save(FESurface& surf, FEDataStream& a)
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
bool FEPlotContactPressure::Save(FESurface &surf, FEDataStream& a)
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
bool FEPlotContactTraction::Save(FESurface &surf, FEDataStream& a)
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
bool FEPlotContactForce::Save(FESurface &surf, FEDataStream &a)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surf);
	if (pcs == 0) return false;
    
	vec3d fn = pcs->GetContactForce();
	a << fn;
    
	return true;
}

//-----------------------------------------------------------------------------
// Plot contact area
bool FEPlotContactArea::Save(FESurface &surf, FEDataStream& a)
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
bool FEPlotContactPenalty::Save(FESurface& surf, FEDataStream& a)
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
bool FEPlotMortarContactGap::Save(FESurface& S, FEDataStream& a)
{
	FEMortarSlidingSurface* ps = dynamic_cast<FEMortarSlidingSurface*>(&S);
	if (ps)
	{
		int N = ps->Nodes();
		for (int i=0; i<N; ++i)
		{
			vec3d vA = ps->m_nu[i];
			vec3d gA = ps->m_gap[i];
			double gap = gA*vA;
			a << gap;
		}
		return true;
	}
	else return false;
}

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Store the average deformation Hessian (G) for each element. 
bool FEPlotElementGnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial2O* pme = dynamic_cast<FEElasticMaterial2O*>(dom.GetMaterial()->GetElasticMaterial());
	if (pme == 0) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();

		tens3drs Gavg; Gavg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint2O& pt2O = *(el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint2O>());
			Gavg += pt2O.m_G;
		}
		Gavg /= (double) nint;

		double L2_norm = sqrt(Gavg.tripledot(Gavg));
		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element. 
bool FEPlotElementStress::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;

	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		mat3ds s; s.zero();
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
				s += pt.m_s;
			}
		}
		s *= f;

		a << s;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average Cauchy stress for each element. 
bool FEPlotElementsnorm::Save(FEDomain& dom, FEDataStream& a)
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
bool FEPlotElementtaunorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEMicroMaterial2O* pme = dynamic_cast<FEMicroMaterial2O*>(dom.GetMaterial()->GetElasticMaterial());
	if (pme == 0) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		tens3ds tau_avg;
		tau_avg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O& pt2O = *(el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//			tau_avg += pt2O.m_tau;
		}
		tau_avg *= f;

		double L2_norm = sqrt(tau_avg.tripledot(tau_avg));

		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK1 stress for each element.
bool FEPlotElementPK1norm::Save(FEDomain& dom, FEDataStream& a)
{
	FEMicroMaterial* pm1O = dynamic_cast<FEMicroMaterial*>(dom.GetMaterial()->GetElasticMaterial());
	if (pm1O)
	{
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int nint = el.GaussPoints();
			double f = 1.0 / (double) nint;

			mat3d PK1_avg; PK1_avg.zero();
			for (int j=0; j<nint; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEMicroMaterialPoint* mmppt = mp.ExtractData<FEMicroMaterialPoint>();
				PK1_avg += pm1O->AveragedStressPK1(mmppt->m_rve, mp);
			}
			PK1_avg *= f;

			double L2_norm = sqrt(PK1_avg.dotdot(PK1_avg));
			a << L2_norm;
		}
		return true;
	}

	FEMicroMaterial2O* pm2O = dynamic_cast<FEMicroMaterial2O*>(dom.GetMaterial()->GetElasticMaterial());
	if (pm2O == 0)
	{
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int nint = el.GaussPoints();
			double f = 1.0 / (double) nint;

			mat3d PK1_avg; PK1_avg.zero();
			for (int j=0; j<nint; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEMicroMaterialPoint2O* mmppt = mp.ExtractData<FEMicroMaterialPoint2O>();
				PK1_avg += mmppt->m_rve.AveragedStressPK1(mp);
			}
			PK1_avg *= f;

			double L2_norm = sqrt(PK1_avg.dotdot(PK1_avg));
			a << L2_norm;
		}
		return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK1 stress moment for each element. 
bool FEPlotElementQK1norm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial2O* pme = dynamic_cast<FEElasticMaterial2O*>(dom.GetMaterial()->GetElasticMaterial());
	if (pme == 0) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
	
		tens3drs QK1_avg; QK1_avg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint2O& pt2O = *(el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint2O>());
			QK1_avg += pt2O.m_Q;
		}
		QK1_avg *= f;

		double L2_norm = sqrt(QK1_avg.tripledot(QK1_avg));
		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK2 stress for each element.
bool FEPlotElementSnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		mat3ds S_avg; S_avg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			
			if (mmppt) S_avg += mmppt->m_S;
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//				if (mmppt2O) S_avg += mmppt2O->m_S;
			}
		}
		S_avg *= f;

		double L2_norm = sqrt(S_avg.dotdot(S_avg));

		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK2 stress moment for each element. 
bool FEPlotElementTnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		tens3ds T_avg; T_avg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//			if (ppt2O) T_avg += ppt2O->m_T;
		}
		T_avg *= f;

		double L2_norm = sqrt(T_avg.tripledot(T_avg));
		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average infinitesimal strain gradient for each element. 
bool FEPlotElementinfstrnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		tens3ds inf_strain_avg; inf_strain_avg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//			if (ppt2O) inf_strain_avg += ppt2O->m_inf_str_grad;
		}
		inf_strain_avg *= f;

		double L2_norm = sqrt(inf_strain_avg.tripledot(inf_strain_avg));
		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average Green-Lagrange strain gradient for each element. 
bool FEPlotElementGLstrnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		tens3ds Havg; Havg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//			if (ppt2O) Havg += ppt2O->m_H;
		}
		Havg *= f;

		double L2_norm = sqrt(Havg.tripledot(Havg));
		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average Euler-Almansi strain gradient for each element. 
bool FEPlotElementEAstrnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		tens3ds havg; havg.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint2O* ppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//			if (ppt2O) havg += ppt2O->m_h;
		}
		havg *= f;

		double L2_norm = sqrt(havg.tripledot(havg));
		a << L2_norm;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Element macro-micro energy difference
bool FEPlotElementenergydiff::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		double energy_diff = 0.;
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			if (mmppt) energy_diff += mmppt->m_energy_diff;	
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//				if (mmppt2O) energy_diff += mmppt2O->m_energy_diff;	
			}
		}
		energy_diff *= f;

		a << energy_diff;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Element macro energy
bool FEPlotElementMacroEnergy::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;
	
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;

		double macro_energy = 0.;
		for (int j=0; j<nint; ++j)
		{
			FEMicroMaterialPoint* mmppt = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
			if (mmppt) macro_energy += mmppt->m_macro_energy;
			else
			{
				FEMicroMaterialPoint2O* mmppt2O = (el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint2O>());
//				if (mmppt2O) macro_energy += mmppt2O->m_macro_energy;
			}
		}
		macro_energy *= f;

		a << macro_energy;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! Element macro energy
bool FEPlotElementMicroEnergy::Save(FEDomain& dom, FEDataStream& a)
{
	FEMicroMaterial* pm1O = dynamic_cast<FEMicroMaterial*>(dom.GetMaterial()->GetElasticMaterial());
	if (pm1O)
	{
		// write solid element data
		int N = dom.Elements();
		for (int i=0; i<N; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int nint = el.GaussPoints();
			double f = 1.0 / (double) nint;

			double micro_energy = 0.;
			for (int j=0; j<nint; ++j)
			{
				FEMicroMaterialPoint& mmpt = *(el.GetMaterialPoint(j)->ExtractData<FEMicroMaterialPoint>());
				micro_energy += mmpt.m_micro_energy;
			}
			micro_energy *= f;
			a << micro_energy;
		}

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store the average elasticity for each element.
bool FEPlotElementElasticity::Save(FEDomain& dom, FEDataStream& a)
{
    FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
    if ((pme == 0) || pme->IsRigid()) return false;
    
	// write solid element data
	int N = dom.Elements();
	for (int i=0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);
        
		tens4ds s(0.0);
		int nint = el.GaussPoints();
		double f = 1.0 / (double) nint;
        
		// we output the average tangent values of the gauss points
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& pt = *el.GetMaterialPoint(j);
            tens4ds c = pme->Tangent(pt);
            s  += c;
		}
		s *= f;
        
		a << s;
	}
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotStrainEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEMaterial* pmm = dom.GetMaterial();
    FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmm);
    if (pme == nullptr) pme = dom.GetMaterial()->GetElasticMaterial();

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
bool FEPlotDevStrainEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEMaterial* pmm = dom.GetMaterial();
    FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmm);
    if (pme == nullptr) pme = dom.GetMaterial()->GetElasticMaterial();
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
bool FEPlotSpecificStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
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
// TODO: Should I call the density for remodeling materials something else? 
//       Or maybe the FEElasticMaterialPoint should define a density parameter
//       that will be updated by the materials to define the current density?
bool FEPlotDensity::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		FEElasticMaterial* em = dynamic_cast<FEElasticMaterial*>(bd.GetMaterial());
		if (em == 0) return false;

		FERemodelingElasticMaterial* rm = dynamic_cast<FERemodelingElasticMaterial*>(em);
		if (rm)
		{
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
			
				a << ew;
			}
			return true;
		}
		else
		{
			double rho0 = em->Density();
			int N = bd.Elements();
			for (int i = 0; i<bd.Elements(); ++i)
			{
				FESolidElement& el = bd.Element(i);
				int nint = el.GaussPoints();

				// calculate average mass density
				double ew = 0;
				for (int j = 0; j<nint; ++j)
				{
					FEElasticMaterialPoint& mp = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
					double J = mp.m_F.det();
					ew += rho0 / J;
				}
				ew /= nint;

				a << ew;
			}
			return true;
		}
	}
	return false;
}


//-----------------------------------------------------------------------------
bool FEPlotRelativeVolume::Save(FEDomain &dom, FEDataStream& a)
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
bool FEPlotFiberVector::Save(FEDomain &dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if (pme == 0) return false;

	int BE = dom.Elements();
	for (int i=0; i<BE; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int n = el.GaussPoints();
		vec3d r = vec3d(0,0,0);
		for (int j=0; j<n; ++j)
		{
			FEElasticMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>();
			vec3d ri = pt.m_Q.col(0);

			r += pt.m_F*ri;
		}
//		r /= (double) n;
		r.unit();

		a << r;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFiberStretch::Save(FEDomain &dom, FEDataStream& a)
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
bool FEPlotDevFiberStretch::Save(FEDomain &dom, FEDataStream& a)
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
bool FEPlotShellThickness::Save(FEDomain &dom, FEDataStream &a)
{
    const int dof_X = GetFEModel()->GetDOFIndex("x");
    const int dof_Y = GetFEModel()->GetDOFIndex("y");
    const int dof_Z = GetFEModel()->GetDOFIndex("z");
	const int dof_U = GetFEModel()->GetDOFIndex("u");
	const int dof_V = GetFEModel()->GetDOFIndex("v");
	const int dof_W = GetFEModel()->GetDOFIndex("w");
	if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
		int NS = sd.Elements();
		FEMesh& mesh = *sd.GetMesh();
		for (int i=0; i<NS; ++i)
		{
			FEShellElement& e = sd.Element(i);
			int n = e.Nodes();
			for (int j=0; j<n; ++j)
			{
				FENode& nj = mesh.Node(e.m_node[j]);
				vec3d D = nj.m_d0 + nj.get_vec3d(dof_X, dof_Y, dof_Z) - nj.get_vec3d(dof_U, dof_V, dof_W);
				double h = D.norm();
				a << h;
			}
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Store shell directors
bool FEPlotShellDirector::Save(FEDomain &dom, FEDataStream &a)
{
    const int dof_X = GetFEModel()->GetDOFIndex("x");
    const int dof_Y = GetFEModel()->GetDOFIndex("y");
    const int dof_Z = GetFEModel()->GetDOFIndex("z");
	const int dof_U = GetFEModel()->GetDOFIndex("u");
	const int dof_V = GetFEModel()->GetDOFIndex("v");
	const int dof_W = GetFEModel()->GetDOFIndex("w");
	if (dom.Class() == FE_DOMAIN_SHELL)
	{
		FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
		int NS = sd.Elements();
		FEMesh& mesh = *sd.GetMesh();
		for (int i=0; i<NS; ++i)
		{
			FEShellElement& e = sd.Element(i);
			int n = e.Nodes();
			for (int j=0; j<n; ++j)
			{
				FENode& nj = mesh.Node(e.m_node[j]);
				vec3d D = nj.m_d0 + nj.get_vec3d(dof_X, dof_Y, dof_Z) - nj.get_vec3d(dof_U, dof_V, dof_W);
				a << D;
			}
		}
		return true;
	}
	return false;
}


//-----------------------------------------------------------------------------
bool FEPlotDamage::Save(FEDomain &dom, FEDataStream& a)
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
FEPlotNestedDamage::FEPlotNestedDamage(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
    m_pfem = pfem;
    m_nmat = -1;
}

//-----------------------------------------------------------------------------
// Resolve nested damage material by number
bool FEPlotNestedDamage::SetFilter(int nmat)
{
    m_nmat = nmat-1;
    return (m_nmat != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotNestedDamage::Save(FEDomain &dom, FEDataStream& a)
{
    int N = dom.Elements();
    FEElasticMaterial* pmat = dom.GetMaterial()->GetElasticMaterial();
    if (dynamic_cast<FEElasticMixture*>(pmat)||dynamic_cast<FEUncoupledElasticMixture*>(pmat))
    {
        int NC = pmat->Properties();
        if ((m_nmat > -1) && (m_nmat < NC))
        {
            for (int i=0; i<N; ++i)
            {
                FEElement& el = dom.ElementRef(i);
                
                float D = 0.f;
                int nint = el.GaussPoints();
                for (int j=0; j<nint; ++j)
                {
                    FEElasticMixtureMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEElasticMixtureMaterialPoint>();
                    FEDamageMaterialPoint* ppd = pt.GetPointData(m_nmat)->ExtractData<FEDamageMaterialPoint>();
                    if (ppd) D += (float) ppd->m_D;
                }
                D /= (float) nint;
                a.push_back(D);
            }
        }
    }
    else if (dynamic_cast<FEElasticMultigeneration*>(pmat))
    {
        FEElasticMultigeneration* pmg = dynamic_cast<FEElasticMultigeneration*>(pmat);
        int NC = pmg->Properties();
        if ((m_nmat > -1) && (m_nmat < NC))
        {
            for (int i=0; i<N; ++i)
            {
                FEElement& el = dom.ElementRef(i);
                
                float D = 0.f;
                int nint = el.GaussPoints();
                for (int j=0; j<nint; ++j)
                {
                    FEMultigenerationMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEMultigenerationMaterialPoint>();
                    FEDamageMaterialPoint* ppd = pt.GetPointData(m_nmat)->ExtractData<FEDamageMaterialPoint>();
                    FEElasticMixtureMaterialPoint* pem = pt.GetPointData(m_nmat)->ExtractData<FEElasticMixtureMaterialPoint>();
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
                D /= (float) nint;
                a.push_back(D);
            }
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
bool FEPlotMixtureVolumeFraction::Save(FEDomain &m, FEDataStream &a)
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
bool FEPlotUT4NodalStresses::Save(FEDomain& dom, FEDataStream& a)
{
	FEUT4Domain* pd = dynamic_cast<FEUT4Domain*>(&dom);
	if (pd == 0) return false;

	int N = pd->Nodes();
	for (int i=0; i<N; ++i)
	{
		FEUT4Domain::UT4NODE& n = pd->UT4Node(i);
		a << n.si;
	}
	return true;
}


//-----------------------------------------------------------------------------
bool FEPlotShellStrain::Save(FEDomain &dom, FEDataStream &a)
{
    if (dom.Class() == FE_DOMAIN_SHELL) {
        FEShellDomain& sd = static_cast<FEShellDomain&>(dom);
        int NE = sd.Elements();
        for (int i=0; i<NE; ++i)
        {
            FEShellElement& el = sd.Element(i);
            int nint = el.GaussPoints();
            mat3ds E; E.zero();
            for (int j=0; j<nint; ++j)
            {
                FEElasticMaterialPoint& pt = *(el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
                E += pt.Strain();
            }
            E /= nint;
            
            a << E;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotSPRStresses::Save(FEDomain& dom, FEDataStream& a)
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
bool FEPlotSPRLinearStresses::Save(FEDomain& dom, FEDataStream& a)
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
	map.SetInterpolationOrder(1);
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
bool FEPlotSPRPrincStresses::Save(FEDomain& dom, FEDataStream& a)
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
bool FEPlotSPRTestLinear::Save(FEDomain& dom, FEDataStream& a)
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
bool FEPlotSPRTestQuadratic::Save(FEDomain& dom, FEDataStream& a)
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
bool FEPlotRigidDisplacement::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());

	// store the rigid body position
	// TODO: why do we not store displacement?
	a << rb.m_rt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidVelocity::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// store the rigid velocity
	a << rb.m_vt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAcceleration::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// store rigid body acceleration
	a << rb.m_at;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidRotation::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    vec3d q = rb.m_qt.GetVector()*rb.m_qt.GetAngle();
    
	// store rotation vector
	a << q;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularVelocity::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// store rigid angular velocity
	a << rb.m_wt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidAngularAcceleration::Save(FEDomain& dom, FEDataStream& a)
{
	// get the rigid material
	FEMaterial* pm = dom.GetMaterial();
	if (pm->IsRigid() == false) return false;
	FERigidMaterial* prm = static_cast<FERigidMaterial*>(pm);
    
	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
	FERigidBody& rb = *rigid.Object(prm->GetRigidBodyID());
    
	// store angular acceleration
	a << rb.m_alt;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidKineticEnergy::Save(FEDomain& dom, FEDataStream& a)
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
    
	// store kinetic energy
	a << ke;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidEuler::Save(FEDomain& dom, FEDataStream& a)
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
    
	// store Euler
	a << E[0] << E[1] << E[2];
    
	return true;
}

//-----------------------------------------------------------------------------
// TODO: I think this already gets stored somewhere
bool FEPlotRigidRotationVector::Save(FEDomain& dom, FEDataStream& a)
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
    
	// store rotation vector
	a << r;
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodalStresses::Save(FEDomain& dom, FEDataStream& a)
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
		for (int j=0; j<ne; ++j) a << s[j];
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Store the average Euler-lagrange strain
bool FEPlotLagrangeStrain::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial* pme = dom.GetMaterial()->GetElasticMaterial();
	if ((pme == 0) || pme->IsRigid()) return false;

	// write solid element data
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		int nint = el.GaussPoints();
		double f = 1.0 / (double)nint;
		mat3dd I(1.0); // identity tensor

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		mat3ds s; s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
			if (ppt)
			{
				mat3d C = ppt->RightCauchyGreen();
				mat3ds E = ((C - I)*0.5).sym();
				s += E;
			}
		}
		s *= f;

		// write average strain
		a << s;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidReactionForce::Save(FEDomain& dom, FEDataStream& a)
{
	// get the material
	FEMaterial* pmat = dom.GetMaterial();
	if ((pmat==0) || (pmat->IsRigid() == false)) return false;

	// get the rigid body ID
	int nrid = pmat->GetRigidBodyID();
	if (nrid < 0) return false;

	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
    FERigidBody& rb = *rigid.Object(nrid);

	a << rb.m_Fr;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRigidReactionTorque::Save(FEDomain& dom, FEDataStream& a)
{
	// get the material
	FEMaterial* pmat = dom.GetMaterial();
	if ((pmat==0) || (pmat->IsRigid() == false)) return false;

	// get the rigid body ID
	int nrid = pmat->GetRigidBodyID();
	if (nrid < 0) return false;

	// get the rigid body
	FERigidSystem& rigid = *m_pfem->GetRigidSystem();
    FERigidBody& rb = *rigid.Object(nrid);

	a << rb.m_Mr;

	return true;
}
