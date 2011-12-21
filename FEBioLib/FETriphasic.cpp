// FETriphasic.cpp: implementation of the FETriphasic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETriphasic.h"
#include "FECore/FEModel.h"

// Material parameters for the FETriphasic material
BEGIN_PARAMETER_LIST(FETriphasic, FEMaterial)
ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
ADD_PARAMETER(m_rhoTw, FE_PARAM_DOUBLE, "fluid_density");
ADD_PARAMETER(m_cFr, FE_PARAM_DOUBLE, "fixed_charge_density");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! FETriphasic constructor

FETriphasic::FETriphasic()
{	m_pPerm = 0;
	m_pOsmC = 0; 
	m_phi0 = 0;
	m_rhoTw = 0;
	m_pSolute[0] = m_pSolute[1] = 0;
	m_cFr = 0;
	m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
}

//-----------------------------------------------------------------------------
void FETriphasic::Init()
{
	FEMaterial::Init();
	m_pSolid->Init();
	m_pPerm->Init();
	m_pOsmC->Init();
	m_pSolute[0]->Init(); m_pSolute[1]->Init();
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 <= phi0 <= 1");
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
	if ((m_pSolute[0]->ChargeNumber() != 1) && (m_pSolute[0]->ChargeNumber() != -1))
		throw MaterialError("charge_number for solute id=1 must be +1 or -1");
	if ((m_pSolute[1]->ChargeNumber() != 1) && (m_pSolute[1]->ChargeNumber() != -1))
		throw MaterialError("charge_number for solute id=2 must be +1 or -1");
	if (m_pSolute[0]->ChargeNumber() != -m_pSolute[1]->ChargeNumber())
		throw MaterialError("charge_number of solutes must have opposite signs");
	
	m_Rgas = FEModel::GetGlobalConstant("R");
	m_Tabs = FEModel::GetGlobalConstant("T");
	m_Fc = FEModel::GetGlobalConstant("Fc");
	
	if (m_Rgas <= 0) throw MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) throw MaterialError("A positive absolute temperature T must be defined in Globals section");
	if (m_Fc <= 0) throw MaterialError("A positive Faraday constant Fc must be defined in Globals section");
	
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FETriphasic::Porosity(FEMaterialPoint& pt)
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
//! Fixed charge density in current configuration
double FETriphasic::FixedChargeDensity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.J;
	double phi0 = pet.m_phi0;
	double cF = m_cFr*(1-phi0)/(J-phi0);
	
	return cF;
}

//-----------------------------------------------------------------------------
//! Electric potential
double FETriphasic::ElectricPotential(FEMaterialPoint& pt)
{
	FESaltMaterialPoint& set = *pt.ExtractData<FESaltMaterialPoint>();
	
	// Fixed charge density
	double cF = FixedChargeDensity(pt);
	// effective concentration
	double c[2] = {set.m_c[0], set.m_c[1]};
	// solubility
	double kappa[2] = {m_pSolute[0]->m_pSolub->Solubility(pt),
		m_pSolute[1]->m_pSolub->Solubility(pt)};
	double z[2] = {m_pSolute[0]->ChargeNumber(),m_pSolute[1]->ChargeNumber()};
	double psi = -m_Rgas*m_Tabs/m_Fc/z[0]
	*log((-cF+sqrt(cF*cF+4*kappa[0]*c[0]*kappa[1]*c[1]))/(2*z[0]*kappa[0]*c[0]));
	
	return psi;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FETriphasic::Concentration(FEMaterialPoint& pt, const int ion)
{
	FESaltMaterialPoint& spt = *pt.ExtractData<FESaltMaterialPoint>();
	
	// electric potential
	double psi = ElectricPotential(pt);
	// effective concentration
	double c = spt.m_c[ion];
	// solubility
	double kappa = m_pSolute[ion]->m_pSolub->Solubility(pt);
	// charge number
	int z = (int) m_pSolute[ion]->ChargeNumber();
	
	// actual concentration
	double ca = kappa*c*exp(-z*m_Fc*psi/m_Rgas/m_Tabs);
	
	return ca;
}

//-----------------------------------------------------------------------------
//! The stress of a triphasic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FETriphasic::Stress(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// fluid pressure
	double p = Pressure(mp);
	
	// add fluid pressure
	s.xx() -= p;
	s.yy() -= p;
	s.zz() -= p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the elastic tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FETriphasic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ept = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESaltMaterialPoint& spt = *mp.ExtractData<FESaltMaterialPoint>();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	
	// relative volume
	double J = ept.J;
	
	// fluid pressure and solute concentration
	double p = Pressure(mp);
	double ca[2] = {Concentration(mp, 0), Concentration(mp, 1)};
	double z[2] = {m_pSolute[0]->m_z, m_pSolute[1]->m_z};
	double cF = FixedChargeDensity(mp);
	double dcFdJ = -cF/(ept.J - ppt.m_phi0);
	
	// solubility and its derivative w.r.t. strain
	double kappa[2] = {m_pSolute[0]->m_pSolub->Solubility(mp),
		m_pSolute[1]->m_pSolub->Solubility(mp)};
	double dkdJ[2] = {m_pSolute[0]->m_pSolub->Tangent_Solubility_Strain(mp),
		m_pSolute[1]->m_pSolub->Tangent_Solubility_Strain(mp)};
	
	// osmotic coefficient and its derivative w.r.t. strain
	double osmc = m_pOsmC->OsmoticCoefficient(mp);
	double dodJ = m_pOsmC->Tangent_OsmoticCoefficient_Strain(mp);
	
	double dp = m_Rgas*m_Tabs*J*(cF*osmc/(ca[0]+ca[1])*
								 (dcFdJ+z[0]*ca[0]*dkdJ[0]/kappa[0]
								  +z[1]*ca[1]*dkdJ[1]/kappa[1])
								 +ca[0]*(dodJ+osmc*dkdJ[0]/kappa[0])
								 +ca[1]*(dodJ+osmc*dkdJ[1]/kappa[1]));
	
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

vec3d FETriphasic::FluidFlux(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESaltMaterialPoint& spt = *pt.ExtractData<FESaltMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c[2] = {spt.m_c[0], spt.m_c[1]};
	
	// concentration gradient
	vec3d gradc[2] = {spt.m_gradc[0], spt.m_gradc[1]};
	
	// hydraulic permeability
	mat3ds kt = m_pPerm->Permeability(pt);
	
	// solute diffusivity in mixture
	mat3ds D[2] = {m_pSolute[0]->m_pDiff->Diffusivity(pt),
		m_pSolute[1]->m_pDiff->Diffusivity(pt)};
	
	// solute free diffusivity
	double D0[2] = {m_pSolute[0]->m_pDiff->Free_Diffusivity(pt),
		m_pSolute[1]->m_pDiff->Free_Diffusivity(pt)};
	
	// solubility
	double kappa[2] = {m_pSolute[0]->m_pSolub->Solubility(pt),
		m_pSolute[1]->m_pSolub->Solubility(pt)};
	
	// identity matrix
	mat3dd I(1);
	
	// effective hydraulic permeability
	mat3ds ke = kt.inverse() + ((I-D[0]/D0[0])*(kappa[0]*c[0]/D0[0])
								+ (I-D[1]/D0[1])*(kappa[1]*c[1]/D0[1]))*(m_Rgas*m_Tabs/phiw);
	ke = ke.inverse();
	
	// fluid flux w
	vec3d w = -(ke*(gradp + ((D[0]*gradc[0])*(kappa[0]/D0[0])
					+ (D[1]*gradc[1])*(kappa[1]/D0[1]))*m_Rgas*m_Tabs));
	
	return w;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FETriphasic::SoluteFlux(FEMaterialPoint& pt, const int ion)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESaltMaterialPoint& spt = *pt.ExtractData<FESaltMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// concentration
	double c = spt.m_c[ion];
	
	// concentration gradient
	vec3d gradc = spt.m_gradc[ion];
	
	// solute diffusivity in mixture
	mat3ds D = m_pSolute[ion]->m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pSolute[ion]->m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double kappa = m_pSolute[ion]->m_pSolub->Solubility(pt);
	
	// fluid flux w
	vec3d w = FluidFlux(pt);
	
	// solute flux j
	vec3d j = D*(w*(c/D0) - gradc*phiw)*kappa;
	
	return j;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FETriphasic::Pressure(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESaltMaterialPoint& spt = *pt.ExtractData<FESaltMaterialPoint>();
	
	// effective pressure
	double p = ppt.m_p;
	
	// effective concentration
	double ca[2] = {Concentration(pt, 0), Concentration(pt, 1)};
	
	// osmotic coefficient
	double osmc = m_pOsmC->OsmoticCoefficient(pt);
	
	// actual pressure
	double pa = p + m_Rgas*m_Tabs*osmc*(ca[0]+ca[1]);
	
	return pa;
}

//-----------------------------------------------------------------------------
//! Current density
vec3d FETriphasic::CurrentDensity(FEMaterialPoint& pt)
{
	vec3d j[2];
	j[0] = SoluteFlux(pt, 0);
	j[1] = SoluteFlux(pt, 1);
	double z[2] = {m_pSolute[0]->ChargeNumber(), m_pSolute[1]->ChargeNumber()};
	
	vec3d Ie = (j[0]*z[0] + j[1]*z[1])*m_Fc;
	
	return Ie;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FETriphasic::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	if (ar.IsSaving())
	{
		ar << m_phi0 << m_rhoTw << m_cFr << m_Rgas << m_Tabs << m_Fc;
		
		ar << febio.GetTypeStr<FEMaterial>(m_pSolid    ); m_pSolid->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pPerm     ); m_pPerm ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pSolute[0]); m_pSolute[0] ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pSolute[1]); m_pSolute[1] ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pOsmC     ); m_pOsmC ->Serialize(ar);
	}
	else
	{
		ar >> m_phi0 >> m_rhoTw >> m_cFr >> m_Rgas >> m_Tabs >> m_Fc;
		
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
		m_pSolute[0] = dynamic_cast<FESolute*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolute[0]); m_pSolute[0]->Serialize(ar);
		m_pSolute[0]->Init();
		
		ar >> sz;
		m_pSolute[1] = dynamic_cast<FESolute*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolute[1]); m_pSolute[1]->Serialize(ar);
		m_pSolute[1]->Init();
		
		ar >> sz;
		m_pOsmC = dynamic_cast<FEOsmoticCoefficient*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pOsmC); m_pOsmC->Serialize(ar);
		m_pOsmC->Init();
		
	}
}
