// FEBiphasic.cpp: implementation of the FEBiphasic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBiphasic.h"

// Material parameters for the FEBiphasic material
BEGIN_PARAMETER_LIST(FEBiphasic, FEMaterial)
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
	ADD_PARAMETER(m_rhoTw, FE_PARAM_DOUBLE, "fluid_density");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FEBiphasic constructor

FEBiphasic::FEBiphasic()
{ 
	m_rhoTw = 0; 
	m_phi0 = 0;

	AddComponent<FEElasticMaterial      >(&m_pSolid, "solid"         );
	AddComponent<FEHydraulicPermeability>(&m_pPerm , "permeability"  );
	AddComponent<FESolventSupply        >(&m_pSupp , "solvent_supply");
}

//-----------------------------------------------------------------------------
void FEBiphasic::Init()
{
	FEMaterial::Init();
	m_pSolid->Init();
	m_pPerm->Init();
	if (m_pSupp) m_pSupp->Init();
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 <= phi0 <= 1");
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.J;
	// porosity
//	double phiw = 1 - m_phi0/J;
	double phi0 = pet.m_phi0;
	double phiw = 1 - phi0/J;
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
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
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
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
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
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
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
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
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
//! serialization
void FEBiphasic::Serialize(DumpFile &ar)
{
	// serialize material parameters
	FEMaterial::Serialize(ar);

	FEBioKernel& febio = FEBioKernel::GetInstance();

	// serialize sub-materials
	if (ar.IsSaving())
	{
		ar << febio.GetTypeStr<FEMaterial>(m_pSolid);
		m_pSolid->Serialize(ar);

		ar << febio.GetTypeStr<FEMaterial>(m_pPerm);
		m_pPerm->Serialize(ar);

		ar << febio.GetTypeStr<FEMaterial>(m_pSupp);
		m_pSupp->Serialize(ar);
	}
	else
	{
		char sz[256] = {0};

		ar >> sz;
		m_pSolid = dynamic_cast<FEElasticMaterial*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolid);
		m_pSolid->Serialize(ar);
		m_pSolid->Init();

		ar >> sz;
		m_pPerm = dynamic_cast<FEHydraulicPermeability*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pPerm);
		m_pPerm->Serialize(ar);
		m_pPerm->Init();

		ar >> sz;
		m_pSupp = dynamic_cast<FESolventSupply*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSupp);
		m_pSupp->Serialize(ar);
		m_pSupp->Init();
	}
}

//-----------------------------------------------------------------------------
FEParam* FEBiphasic::GetParameter(const ParamString& s)
{
	// see if this is a composite parameter
	if (s.count() == 1) return FEMultiMaterial::GetParameter(s);

	// else find the component's parameter
	if      (s == "solid"       ) return m_pSolid->GetParameter(s.next());
	else if (s == "permeability") return m_pPerm ->GetParameter(s.next());
	else return 0;
}

//-----------------------------------------------------------------------------
// Material parameters for FEHydraulicPermeability
void FEHydraulicPermeability::Init()
{
	FEMaterial::Init();
}

//-----------------------------------------------------------------------------
// Derivative of permeability w.r.t. solute concentration at material point
// Set this to zero by default because poroelasticity problems do not require it
mat3ds FEHydraulicPermeability::Tangent_Permeability_Concentration(FEMaterialPoint& pt, const int isol)
{
	return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
// Material parameters for FESolventSupply
void FESolventSupply::Init()
{
	FEMaterial::Init();
}

//-----------------------------------------------------------------------------
// Derivative of supply w.r.t. solute concentration at material point
// Set this to zero by default because biphasic problems do not require it
double FESolventSupply::Tangent_Supply_Concentration(FEMaterialPoint& pt, const int isol)
{
	return 0;
}
