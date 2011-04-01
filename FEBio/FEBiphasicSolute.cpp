// FESolutePoroElastic.cpp: implementation of the FESolutePoroElastic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMaterial.h"
#include "fem.h"

// register the material with the framework
REGISTER_MATERIAL(FEBiphasicSolute, "biphasic-solute");

// Material parameters for the FESolutePoroElastic material
BEGIN_PARAMETER_LIST(FEBiphasicSolute, FEMaterial)
	ADD_PARAMETER(m_rhoTw, FE_PARAM_DOUBLE, "fluid_density");
	ADD_PARAMETER(m_rhoTu, FE_PARAM_DOUBLE, "solute_density");
	ADD_PARAMETER(m_Mu, FE_PARAM_DOUBLE, "solute_molecular_weight");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FEBiphasicSolute constructor

FEBiphasicSolute::FEBiphasicSolute()
{	m_pPerm = 0; m_pDiff = 0; m_pSolub = 0; m_pOsmC= 0; 
	m_rhoTw = 0; m_rhoTu = 0; m_Mu = 0; m_Rgas = 0; m_Tabs = 0; }

//-----------------------------------------------------------------------------
void FEBiphasicSolute::Init()
{
	FEMaterial::Init();
	m_pSolid->Init();
	m_pPerm->Init();
	m_pDiff->Init();
	m_pSolub->Init();
	m_pOsmC->Init();
	
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
	if (m_rhoTu < 0) throw MaterialError("solute_density must be positive");
	if (m_Mu < 0) throw MaterialError("solute_molecular_weight must be positive");
	
	m_Rgas = FEM::GetGlobalConstant("R");
	m_Tabs = FEM::GetGlobalConstant("T");
	
	if (m_Rgas <= 0) throw MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) throw MaterialError("A positive absolute temperature T must be defined in Globals section");
	
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasicSolute::Porosity(FEMaterialPoint& pt)
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
//! The stress of a solute-poroelastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasicSolute::Stress(FEMaterialPoint& mp)
{
	FESolutePoroElasticMaterialPoint& pt = *mp.ExtractData<FESolutePoroElasticMaterialPoint>();
	
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// add fluid pressure
	s.xx() -= pt.m_pa;
	s.yy() -= pt.m_pa;
	s.zz() -= pt.m_pa;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the elastic tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEBiphasicSolute::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ept = *mp.ExtractData<FEElasticMaterialPoint>();
	FESolutePoroElasticMaterialPoint& pt = *mp.ExtractData<FESolutePoroElasticMaterialPoint>();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	
	// relative volume
	double J = ept.J;
	
	// fluid pressure and solute concentration
	double p = pt.m_pa;
	double c = pt.m_c;
	
	// solubility and its derivative w.r.t. strain
	double kappa = m_pSolub->Solubility(mp);
	double dkdJ = m_pSolub->Tangent_Solubility_Strain(mp);
	
	// osmotic coefficient and its derivative w.r.t. strain
	double osmc = m_pOsmC->OsmoticCoefficient(mp);
	double dodJ = m_pOsmC->Tangent_OsmoticCoefficient_Strain(mp);
	
	double dp = m_Rgas*m_Tabs*c*J*(dodJ*kappa+osmc*dkdJ);
	
	// adjust tangent for pressures
	double D[6][6] = {0};
	C.extract(D);
	
	D[0][0] -= -p + dp;
	D[1][1] -= -p + dp;
	D[2][2] -= -p + dp;
	
	D[0][1] -= p + dp; D[1][0] -= p + dp;
	D[1][2] -= p + dp; D[2][1] -= p + dp;
	D[0][2] -= p + dp; D[2][0] -= p + dp;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! Calculate fluid flux

vec3d FEBiphasicSolute::FluidFlux(FEMaterialPoint& pt)
{
	FESolutePoroElasticMaterialPoint& ppt = *pt.ExtractData<FESolutePoroElasticMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c = ppt.m_c;
	
	// concentration gradient
	vec3d gradc = ppt.m_gradc;
	
	// hydraulic permeability
	mat3ds kt = m_pPerm->Permeability(pt);
	
	// solute diffusivity in mixture
	mat3ds D = m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double kappa = m_pSolub->Solubility(pt);
	
	// identity matrix
	mat3dd I(1);
	
	// effective hydraulic permeability
	mat3ds ke = kt.inverse() + (I-D/D0)*(m_Rgas*m_Tabs*kappa*c/phiw/D0);
	ke = ke.inverse();
	
	// fluid flux w
	vec3d w = -(ke*(gradp + (D*gradc)*(m_Rgas*m_Tabs*kappa/D0)));
	
	return w;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEBiphasicSolute::SoluteFlux(FEMaterialPoint& pt)
{
	FESolutePoroElasticMaterialPoint& ppt = *pt.ExtractData<FESolutePoroElasticMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c = ppt.m_c;
	
	// concentration gradient
	vec3d gradc = ppt.m_gradc;
	
	// solute diffusivity in mixture
	mat3ds D = m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double kappa = m_pSolub->Solubility(pt);
	
	// fluid flux w
	vec3d w = FluidFlux(pt);
	
	// solute flux j
	vec3d j = D*(w*(c/D0) - gradc*phiw)*kappa;
	
	return j;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEBiphasicSolute::Pressure(FEMaterialPoint& pt)
{
	FESolutePoroElasticMaterialPoint& ppt = *pt.ExtractData<FESolutePoroElasticMaterialPoint>();
	
	// effective pressure
	double p = ppt.m_p;
	
	// effective concentration
	double c = ppt.m_c;
	
	// osmotic coefficient
	double osmc = m_pOsmC->OsmoticCoefficient(pt);
	
	// solubility
	double kappa = m_pSolub->Solubility(pt);
	
	// actual pressure
	double pa = p + m_Rgas*m_Tabs*osmc*kappa*c;
	
	return pa;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEBiphasicSolute::Concentration(FEMaterialPoint& pt)
{
	FESolutePoroElasticMaterialPoint& ppt = *pt.ExtractData<FESolutePoroElasticMaterialPoint>();
	
	// solubility
	double kappa = m_pSolub->Solubility(pt);
	
	// actual concentration = solubility * effective concentration
	double ca = kappa*ppt.m_c;
	
	return ca;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEBiphasicSolute::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_Rgas << m_Tabs;

		ar << m_pSolid->GetTypeString(); m_pSolid->Serialize(ar);
		ar << m_pPerm ->GetTypeString(); m_pPerm ->Serialize(ar);
		ar << m_pDiff ->GetTypeString(); m_pDiff ->Serialize(ar);
		ar << m_pSolub->GetTypeString(); m_pSolub->Serialize(ar);
		ar << m_pOsmC ->GetTypeString(); m_pOsmC ->Serialize(ar);

	}
	else
	{
		ar >> m_Rgas >> m_Tabs;

		char sz[256] = {0};
		ar >> sz;
		m_pSolid = dynamic_cast<FEElasticMaterial*>(FEMaterialFactory::CreateMaterial(sz));
		assert(m_pSolid); m_pSolid->Serialize(ar);
		m_pSolid->Init();

		ar >> sz;
		m_pPerm = dynamic_cast<FEHydraulicPermeability*>(FEMaterialFactory::CreateMaterial(sz));
		assert(m_pPerm); m_pPerm->Serialize(ar);
		m_pPerm->Init();

		ar >> sz;
		m_pDiff = dynamic_cast<FESoluteDiffusivity*>(FEMaterialFactory::CreateMaterial(sz));
		assert(m_pDiff); m_pDiff->Serialize(ar);
		m_pDiff->Init();

		ar >> sz;
		m_pSolub = dynamic_cast<FESoluteSolubility*>(FEMaterialFactory::CreateMaterial(sz));
		assert(m_pSolub); m_pSolub->Serialize(ar);
		m_pSolub->Init();

		ar >> sz;
		m_pOsmC = dynamic_cast<FEOsmoticCoefficient*>(FEMaterialFactory::CreateMaterial(sz));
		assert(m_pOsmC); m_pOsmC->Serialize(ar);
		m_pOsmC->Init();
	}
}
