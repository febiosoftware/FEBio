// FEBiphasic.cpp: implementation of the FEBiphasic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBiphasic.h"
#include "fem.h"

// register the material with the framework
REGISTER_MATERIAL(FEBiphasic, "biphasic");

// Material parameters for the FEBiphasic material
BEGIN_PARAMETER_LIST(FEBiphasic, FEMaterial)
	ADD_PARAMETER(m_rhoTw, FE_PARAM_DOUBLE, "fluid_density");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEBiphasic::FEBiphasic()
{ m_pSolid = 0; m_pPerm = 0; m_rhoTw = 0; }

//-----------------------------------------------------------------------------
void FEBiphasic::Init()
{
	FEMaterial::Init();
	m_pSolid->Init();
	m_pPerm->Init();
	
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	
	// relative volume
	double J = et.J;
	// porosity
	double phiw = 1 - m_pPerm->m_phi0/J;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! The stress of a poro-elastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasic::Stress(FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// add fluid pressure
	s.xx() -= pt.m_p;
	s.yy() -= pt.m_p;
	s.zz() -= pt.m_p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the sum of the elastic tangent plus the fluid tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEBiphasic::Tangent(FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
	
	// call solid tangent routine
	tens4ds c = m_pSolid->Tangent(mp);
	
	// fluid pressure
	double p = pt.m_p;
	
	// adjust tangent for pressures
	double D[6][6] = {0};
	c.extract(D);
	
	D[0][0] -= -p;
	D[1][1] -= -p;
	D[2][2] -= -p;
	
	D[0][1] -= p; D[1][0] -= p;
	D[1][2] -= p; D[2][1] -= p;
	D[0][2] -= p; D[2][0] -= p;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! Calculate fluid flux from the hydraulic permeability and the fluid pressure
//! gradient

vec3d FEBiphasic::Flux(FEMaterialPoint& pt)
{
	FEPoroElasticMaterialPoint& ppt = *pt.ExtractData<FEPoroElasticMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	mat3ds kt = m_pPerm->Permeability(pt);
	
	return -(kt*gradp);
}

//-----------------------------------------------------------------------------
//! actual fluid pressure (same as effective pressure here)

double FEBiphasic::Pressure(FEMaterialPoint& pt)
{
	FEPoroElasticMaterialPoint& ppt = *pt.ExtractData<FEPoroElasticMaterialPoint>();
	
	return ppt.m_p;
}

//-----------------------------------------------------------------------------
//! Return the permeability tensor as a double array

void FEBiphasic::Permeability(double k[3][3], FEMaterialPoint& pt)

{
	mat3ds kt = m_pPerm->Permeability(pt);
	
	k[0][0] = kt.xx();
	k[1][1] = kt.yy();
	k[2][2] = kt.zz();
	k[0][1] = k[1][0] = kt.xy();
	k[1][2] = k[2][1] = kt.yz();
	k[2][0] = k[0][2] = kt.xz();
	
}

//-----------------------------------------------------------------------------
//! tangent of permeability

tens4ds FEBiphasic::Tangent_Permeability_Strain(FEMaterialPoint& pt)
{
	return m_pPerm->Tangent_Permeability_Strain(pt);
}


//-----------------------------------------------------------------------------
//! serialization
void FEBiphasic::Serialize(DumpFile &ar)
{
	// serialize material parameters
	FEMaterial::Serialize(ar);

	// serialize sub-materials
	if (ar.IsSaving())
	{
		ar << m_pSolid->GetTypeString();
		m_pSolid->Serialize(ar);

		ar << m_pPerm->GetTypeString();
		m_pPerm->Serialize(ar);
	}
	else
	{
		char sz[256] = {0};

		FEBioKernel& febio = FEBioKernel::GetInstance();

		ar >> sz;
		m_pSolid = dynamic_cast<FEElasticMaterial*>(febio.CreateMaterial(sz, ar.GetFEM()));
		assert(m_pSolid);
		m_pSolid->Serialize(ar);
		m_pSolid->Init();

		ar >> sz;
		m_pPerm = dynamic_cast<FEHydraulicPermeability*>(febio.CreateMaterial(sz, ar.GetFEM()));
		assert(m_pPerm);
		m_pPerm->Serialize(ar);
		m_pPerm->Init();
	}
}


//-----------------------------------------------------------------------------
// Material parameters for FEHydraulicPermeability
BEGIN_PARAMETER_LIST(FEHydraulicPermeability, FEMaterial)
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
END_PARAMETER_LIST();

void FEHydraulicPermeability::Init()
{
	FEMaterial::Init();
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
}

//-----------------------------------------------------------------------------
// Derivative of permeability w.r.t. solute concentration at material point
// Set this to zero by default because poroelasticity problems do not require it
mat3ds FEHydraulicPermeability::Tangent_Permeability_Concentration(FEMaterialPoint& pt)
{
	return mat3ds(0,0,0,0,0,0);
}
