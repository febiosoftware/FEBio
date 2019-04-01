#include "stdafx.h"
#include "FESupplyMichaelisMenten.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// define the material parameters
BEGIN_FECORE_CLASS(FESupplyMichaelisMenten, FESoluteSupply)
	ADD_PARAMETER(m_Vmax, FE_RANGE_GREATER_OR_EQUAL(0.0), "Vmax");
	ADD_PARAMETER(m_Km  , FE_RANGE_GREATER         (0.0), "Km"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplyMichaelisMenten::FESupplyMichaelisMenten(FEModel* pfem) : FESoluteSupply(pfem)
{
	m_Vmax = m_Km = 0;
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESupplyMichaelisMenten::Supply(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca[0];
	double phi0 = ppt.m_phi0;
	double cr = (J-phi0)*ca;
	double crhat = -m_Vmax*cr/(m_Km+cr);
	
	return crhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESupplyMichaelisMenten::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to referential concentration
double FESupplyMichaelisMenten::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca[0];
	double phi0 = ppt.m_phi0;
	double cr = (J-phi0)*ca;
	double dcrhatdcr = -m_Vmax*m_Km/SQR(m_Km+cr);
	
	return dcrhatdcr;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplyMichaelisMenten::ReceptorLigandSupply(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Solute supply at steady state
double FESupplyMichaelisMenten::SupplySS(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand concentration at steady-state
double FESupplyMichaelisMenten::ReceptorLigandConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Referential solid supply (moles of solid/referential volume/time)
double FESupplyMichaelisMenten::SolidSupply(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Referential solid concentration (moles of solid/referential volume)
//! at steady-state
double FESupplyMichaelisMenten::SolidConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}

