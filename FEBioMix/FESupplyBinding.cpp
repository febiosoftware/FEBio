/*
 *  FESupplyBinding.cpp
 *
 *  Created by Gerard Ateshian on 8/8/11.
 *
 */

#include "FESupplyBinding.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FESupplyBinding, FESoluteSupply)
	ADD_PARAMETER(m_kf, FE_PARAM_DOUBLE, "kf");
	ADD_PARAMETER(m_kr, FE_PARAM_DOUBLE, "kr");
	ADD_PARAMETER(m_crt, FE_PARAM_DOUBLE, "Rtot");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplyBinding::FESupplyBinding()
{
	m_kf = m_kr = m_crt = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESupplyBinding::Init()
{
	if (m_kf <= 0) throw MaterialError("kf must be > 0");
	if (m_kr < 0) throw MaterialError("kr must be >= 0");
	if (m_crt < 0) throw MaterialError("Rtot must be >= 0");
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESupplyBinding::Supply(FEMaterialPoint& mp)
{
	double crhat = -ReceptorLigandSupply(mp);
	return crhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESupplyBinding::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to referential concentration
double FESupplyBinding::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	FESoluteMaterialPoint& pt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	double crc = pt.m_crc;
	double dcrhatdcr = -m_kf*(m_crt - crc);
	
	return dcrhatdcr;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplyBinding::ReceptorLigandSupply(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca;
	double phi0 = ppt.m_phi0;
	double cr = (J-phi0)*ca;
	double crc = spt.m_crc;
	double crchat = m_kf*cr*(m_crt-crc) - m_kr*crc;
	return crchat;
}

//-----------------------------------------------------------------------------
//! Solute supply at steady state
double FESupplyBinding::SupplySS(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand concentration at steady-state
double FESupplyBinding::ReceptorLigandConcentrationSS(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	double J = et.m_J;
	double ca = spt.m_ca;
	double phi0 = ppt.m_phi0;
	double cr = (J-phi0)*ca;
	double Kd = m_kr/m_kf;	// dissociation constant
	double crc = m_crt*cr/(Kd+cr);
	return crc;
}

//-----------------------------------------------------------------------------
//! Referential solid supply (moles of solid/referential volume/time)
double FESupplyBinding::SolidSupply(FEMaterialPoint& mp)
{
	return ReceptorLigandSupply(mp);
}

//-----------------------------------------------------------------------------
//! Referential solid concentration (moles of solid/referential volume)
//! at steady-state
double FESupplyBinding::SolidConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}

