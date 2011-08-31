/*
 *  FESupplySynthesisBinding.cpp
 *
 *  Created by Gerard Ateshian on 8/18/11.
 *
 */

#include "stdafx.h"
#include "FESupplySynthesisBinding.h"

// register the material with the framework
REGISTER_MATERIAL(FESupplySynthesisBinding, "supply-synthesis-binding");

// define the material parameters
BEGIN_PARAMETER_LIST(FESupplySynthesisBinding, FESoluteSupply)
	ADD_PARAMETER(m_supp, FE_PARAM_DOUBLE, "supp");
	ADD_PARAMETER(m_kf, FE_PARAM_DOUBLE, "kf");
	ADD_PARAMETER(m_kr, FE_PARAM_DOUBLE, "kr");
	ADD_PARAMETER(m_crt, FE_PARAM_DOUBLE, "Rtot");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplySynthesisBinding::FESupplySynthesisBinding()
{
	m_supp = m_kf = m_kr = m_crt = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESupplySynthesisBinding::Init()
{
	if (m_kf <= 0) throw MaterialError("kf must be > 0");
	if (m_kr < 0) throw MaterialError("kr must be >= 0");
	if (m_crt < 0) throw MaterialError("Rtot must be >= 0");
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
	FESoluteMaterialPoint& pt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	double crc = pt.m_crc;
	double dcrhatdcr = -m_kf*(m_crt - crc);
	
	return dcrhatdcr;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplySynthesisBinding::ReceptorLigandSupply(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESoluteMaterialPoint& spt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	double J = et.J;
	double ca = spt.m_ca;
	double phiw = ppt.m_phiw;
	double cr = phiw*J*ca;
	double crc = spt.m_crc;
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
	FESoluteMaterialPoint& spt = *mp.ExtractData<FESoluteMaterialPoint>();
	
	double J = et.J;
	double ca = spt.m_ca;
	double phiw = ppt.m_phiw;
	double cr = phiw*J*ca;
	double Kd = m_kr/m_kf;	// dissociation constant
	double crc = m_crt*cr/(Kd+cr);
	return crc;
}
