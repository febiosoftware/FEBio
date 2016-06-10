#include "stdafx.h"
#include "FEUncoupledViscoElasticMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEUncoupledViscoElasticMaterial, FEUncoupledMaterial)
	ADD_PARAMETER(m_t[0], FE_PARAM_DOUBLE, "t1");
	ADD_PARAMETER(m_t[1], FE_PARAM_DOUBLE, "t2");
	ADD_PARAMETER(m_t[2], FE_PARAM_DOUBLE, "t3");
	ADD_PARAMETER(m_t[3], FE_PARAM_DOUBLE, "t4");
	ADD_PARAMETER(m_t[4], FE_PARAM_DOUBLE, "t5");
	ADD_PARAMETER(m_t[5], FE_PARAM_DOUBLE, "t6");
	ADD_PARAMETER(m_g0  , FE_PARAM_DOUBLE, "g0");
	ADD_PARAMETER(m_g[0], FE_PARAM_DOUBLE, "g1");
	ADD_PARAMETER(m_g[1], FE_PARAM_DOUBLE, "g2");
	ADD_PARAMETER(m_g[2], FE_PARAM_DOUBLE, "g3");
	ADD_PARAMETER(m_g[3], FE_PARAM_DOUBLE, "g4");
	ADD_PARAMETER(m_g[4], FE_PARAM_DOUBLE, "g5");
	ADD_PARAMETER(m_g[5], FE_PARAM_DOUBLE, "g6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEUncoupledViscoElasticMaterial::FEUncoupledViscoElasticMaterial(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_g0 = 1;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		m_t[i] = 1;
		m_g[i] = 0;
	}
	m_binit = false;

	AddProperty(&m_pBase, "elastic");
}

//-----------------------------------------------------------------------------
void FEUncoupledViscoElasticMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
	FEElasticMaterial* pme2 = GetBaseMaterial();
	pme2->SetLocalCoordinateSystem(el, n, mp);
}

//-----------------------------------------------------------------------------
//! data initialization and checking
bool FEUncoupledViscoElasticMaterial::Init()
{
	// combine bulk modulus from base material and uncoupled viscoelastic material
	if (m_binit == false) m_K += m_pBase->m_K;

	if (FEUncoupledMaterial::Init() == false) return false;

	m_binit = true;

	return true;
}

//-----------------------------------------------------------------------------
//! Create material point data
FEMaterialPoint* FEUncoupledViscoElasticMaterial::CreateMaterialPointData()
{ 
	return new FEViscoElasticMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEUncoupledViscoElasticMaterial::DevStress(FEMaterialPoint& mp)
{
	double dt = GetFEModel()->GetTime().timeIncrement;
	if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();
	
	// Calculate the new elastic Cauchy stress
	pt.m_se = m_pBase->DevStress(mp);
	
	// pull-back to get PK2 stress
	mat3ds Se = ep.pull_back(pt.m_se);
	
	// get elastic PK2 stress of previous timestep
	mat3ds Sep = pt.m_Sep;
	
	// calculate new history variables
	// terms are accumulated in S, the total PK2-stress
	mat3ds S = Se*m_g0;
	double g, h;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = (1 - g)/(dt/m_t[i]);
		
		pt.m_H[i] = pt.m_Hp[i]*g + (Se - Sep)*h;
		S += pt.m_H[i]*m_g[i];
	}
	
	// return the total Cauchy stress,
	// which is the push-forward of S
	return ep.push_forward(S);
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4ds FEUncoupledViscoElasticMaterial::DevTangent(FEMaterialPoint& pt)
{
	double dt = GetFEModel()->GetTime().timeIncrement;

	// calculate the spatial elastic tangent
	tens4ds C = m_pBase->DevTangent(pt);
	if (dt == 0.0) return C;
	
	// calculate the visco scale factor
	double f = m_g0, g, h;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = ( 1 - exp(-dt/m_t[i]) )/( dt/m_t[i] );
		f += m_g[i]*h; 
	}
	
	// multiply tangent with visco-factor
	return C*f;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FEUncoupledViscoElasticMaterial::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
/*    if (mp.dt == 0) return 0;
    
	// get the viscoelastic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();
	
	// Calculate the new elastic strain energy density
	pt.m_sed = m_pBase->DevStrainEnergyDensity(mp);
    double sed = pt.m_sed;
	
	// get elastic strain energy density of previous timestep
	double sedp = pt.m_sedp;
	
	// calculate new history variables
	// terms are accumulated in sedt, the total strain energy density
	double sedt = sed*m_g0;
	double dt = mp.dt, g, h;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = (1 - g)/(dt/m_t[i]);
		
		pt.m_Hsed[i] = pt.m_Hsedp[i]*g + (sed - sedp)*h;
		sedt += pt.m_Hsed[i]*m_g[i];
	}
	
	// return the total strain energy density
	return sedt;*/
    return 0;
}
