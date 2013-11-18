#include "stdafx.h"
#include "FEUncoupledViscoElasticMaterial.h"

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
	m_pBase = 0;
	m_binit = false;
}

//-----------------------------------------------------------------------------
//! data initialization and checking
void FEUncoupledViscoElasticMaterial::Init()
{
	FEUncoupledMaterial::Init();
	if (m_pBase == 0) throw MaterialError("This material needs a base material.");
	m_pBase->Init();
	
	// combine bulk modulus from base material and uncoupled viscoelastic material
	if (m_binit == false) m_K += m_pBase->m_K;

	m_binit = true;
}

//-----------------------------------------------------------------------------
//! This material only has one property
int FEUncoupledViscoElasticMaterial::Properties()
{
	return 1;
}

//-----------------------------------------------------------------------------
FEMaterial* FEUncoupledViscoElasticMaterial::GetProperty(int i)
{
	if (i == 0) return m_pBase;
	assert(false);
	return 0;
}

//-----------------------------------------------------------------------------
//! Create material point data
FEMaterialPoint* FEUncoupledViscoElasticMaterial::CreateMaterialPointData()
{ 
	return new FEViscoElasticMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEUncoupledViscoElasticMaterial::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "elastic") == 0) return 0; else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEUncoupledViscoElasticMaterial::SetProperty(int i, FEMaterial* pm)
{
	if (i==0)
	{
		FEUncoupledMaterial* pme = dynamic_cast<FEUncoupledMaterial*>(pm);
		if (pme)
		{ 
			SetBaseMaterial(pme);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEUncoupledViscoElasticMaterial::DevStress(FEMaterialPoint& mp)
{
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
	double dt = mp.dt, g, h;
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
	// calculate the spatial elastic tangent
	tens4ds C = m_pBase->DevTangent(pt);
	
	// calculate the visco scale factor
	double dt = pt.dt;
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
//! Get a material parameter
FEParam* FEUncoupledViscoElasticMaterial::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEUncoupledMaterial::GetParameter(s);
	if (s == "elastic") return m_pBase->GetParameter(s.next());
	return 0;
}
