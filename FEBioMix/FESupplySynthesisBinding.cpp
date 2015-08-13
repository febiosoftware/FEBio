/*
 *  FESupplySynthesisBinding.cpp
 *
 *  Created by Gerard Ateshian on 8/18/11.
 *
 */

#include "FESupplySynthesisBinding.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FESupplySynthesisBinding, FESoluteSupply)
	ADD_PARAMETER(m_supp, FE_PARAM_DOUBLE, "supp");
	ADD_PARAMETER2(m_kf , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "kf"  );
	ADD_PARAMETER2(m_kr , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "kr"  );
	ADD_PARAMETER2(m_crt, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Rtot");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplySynthesisBinding::FESupplySynthesisBinding(FEModel* pfem) : FESoluteSupply(pfem)
{
	m_supp = m_kf = m_kr = m_crt = 0;
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESupplySynthesisBinding::Supply(FEMaterialPoint& mp)
{
	double crhat = m_supp-ReceptorLigandSupply(mp);
	return crhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESupplySynthesisBinding::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to referential concentration
double FESupplySynthesisBinding::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	FESolutesMaterialPoint& pt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double crc = pt.m_sbmr[0];
	double dcrhatdcr = -m_kf*(m_crt - crc);
	
	return dcrhatdcr;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplySynthesisBinding::ReceptorLigandSupply(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca[0];
	double phi0 = ppt.m_phi0;
	double cr = (J-phi0)*ca;
	double crc = spt.m_sbmr[0];
	double crchat = m_kf*cr*(m_crt-crc) - m_kr*crc;
	return crchat;
}

//-----------------------------------------------------------------------------
//! Solute supply at steady state
double FESupplySynthesisBinding::SupplySS(FEMaterialPoint& mp)
{
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand concentration at steady-state
double FESupplySynthesisBinding::ReceptorLigandConcentrationSS(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca[0];
	double phi0 = ppt.m_phi0;
	double cr = (J-phi0)*ca;
	double Kd = m_kr/m_kf;	// dissociation constant
	double crc = m_crt*cr/(Kd+cr);
	return crc;
}

//-----------------------------------------------------------------------------
//! Referential solid supply (moles of solid/referential volume/time)
double FESupplySynthesisBinding::SolidSupply(FEMaterialPoint& mp)
{
	return ReceptorLigandSupply(mp);
}

//-----------------------------------------------------------------------------
//! Referential solid concentration (moles of solid/referential volume)
//! at steady-state
double FESupplySynthesisBinding::SolidConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}
