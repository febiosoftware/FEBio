// FEBiphasicSolute.cpp: implementation of the FEBiphasicSolute class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBiphasicSolute.h"
#include "FECore/FEModel.h"

// Material parameters for the FEBiphasicSolute material
BEGIN_PARAMETER_LIST(FEBiphasicSolute, FEMaterial)
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
	ADD_PARAMETER(m_rhoTw, FE_PARAM_DOUBLE, "fluid_density");
	ADD_PARAMETER(m_rhoTu, FE_PARAM_DOUBLE, "solute_density");
	ADD_PARAMETER(m_Mu, FE_PARAM_DOUBLE, "solute_molar_mass");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FEBiphasicSolute constructor

FEBiphasicSolute::FEBiphasicSolute()
{
	m_phi0 = 0;
	m_rhoTw = 0;
	m_rhoTu = 0;
	m_Mu = 0;
	m_Rgas = 0;
	m_Tabs = 0; 

	AddComponent<FEElasticMaterial      >(&m_pSolid, "solid"              );
	AddComponent<FEHydraulicPermeability>(&m_pPerm , "permeability"       );
	AddComponent<FESoluteDiffusivity    >(&m_pDiff , "diffusivity"        );
	AddComponent<FESoluteSolubility     >(&m_pSolub, "solubility"         );
	AddComponent<FEOsmoticCoefficient   >(&m_pOsmC , "osmotic_coefficient");
	AddComponent<FESoluteSupply         >(&m_pSupp , "supply"             ); 
}

//-----------------------------------------------------------------------------
void FEBiphasicSolute::Init()
{
	FEMaterial::Init();
	m_pSolid->Init();
	m_pPerm->Init();
	m_pDiff->Init();
	m_pSolub->Init();
	m_pOsmC->Init();
	if (m_pSupp) m_pSupp->Init();
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 <= phi0 <= 1");
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
	if (m_rhoTu < 0) throw MaterialError("solute_density must be positive");
	if (m_Mu < 0) throw MaterialError("solute_molar_mass must be positive");
	
	m_Rgas = FEModel::GetGlobalConstant("R");
	m_Tabs = FEModel::GetGlobalConstant("T");
	
	if (m_Rgas <= 0) throw MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) throw MaterialError("A positive absolute temperature T must be defined in Globals section");
	
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEBiphasicSolute::Porosity(FEMaterialPoint& pt)
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
//! The stress of a solute-poroelastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEBiphasicSolute::Stress(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
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
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	
	// relative volume
	double J = ept.J;
	
	// fluid pressure and solute concentration
	double p = ppt.m_pa;
	double c = spt.m_c;
	
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
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *pt.ExtractData<FESoluteMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c = spt.m_c;
	
	// concentration gradient
	vec3d gradc = spt.m_gradc;
	
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
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *pt.ExtractData<FESoluteMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c = spt.m_c;
	
	// concentration gradient
	vec3d gradc = spt.m_gradc;
	
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
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *pt.ExtractData<FESoluteMaterialPoint>();
	
	// effective pressure
	double p = ppt.m_p;
	
	// effective concentration
	double c = spt.m_c;
	
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
	FESoluteMaterialPoint& spt = *pt.ExtractData<FESoluteMaterialPoint>();
	
	// solubility
	double kappa = m_pSolub->Solubility(pt);
	
	// actual concentration = solubility * effective concentration
	double ca = kappa*spt.m_c;
	
	return ca;
}

//-----------------------------------------------------------------------------
//! referential solute concentration
double FEBiphasicSolute::ReferentialConcentration(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();

	double J = ept.J;
	double phiw = Porosity(pt);
	double cr = J*phiw*Concentration(pt);
	
	return cr;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEBiphasicSolute::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);
	FEBioKernel& febio = FEBioKernel::GetInstance();

	if (ar.IsSaving())
	{
		ar << m_Rgas << m_Tabs;

		ar << febio.GetTypeStr<FEMaterial>(m_pSolid); m_pSolid->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pPerm ); m_pPerm ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pDiff ); m_pDiff ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pSolub); m_pSolub->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pOsmC ); m_pOsmC ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pSupp ); m_pSupp ->Serialize(ar);
	}
	else
	{
		ar >> m_Rgas >> m_Tabs;

		char sz[256] = {0};
		ar >> sz;
		m_pSolid = dynamic_cast<FEElasticMaterial*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolid); m_pSolid->Serialize(ar);
		m_pSolid->Init();

		ar >> sz;
		m_pPerm = dynamic_cast<FEHydraulicPermeability*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pPerm); m_pPerm->Serialize(ar);
		m_pPerm->Init();

		ar >> sz;
		m_pDiff = dynamic_cast<FESoluteDiffusivity*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pDiff); m_pDiff->Serialize(ar);
		m_pDiff->Init();

		ar >> sz;
		m_pSolub = dynamic_cast<FESoluteSolubility*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolub); m_pSolub->Serialize(ar);
		m_pSolub->Init();

		ar >> sz;
		m_pOsmC = dynamic_cast<FEOsmoticCoefficient*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pOsmC); m_pOsmC->Serialize(ar);
		m_pOsmC->Init();

		ar >> sz;
		m_pSupp = dynamic_cast<FESoluteSupply*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSupp); m_pSupp->Serialize(ar);
		m_pSupp->Init();

	}
}
