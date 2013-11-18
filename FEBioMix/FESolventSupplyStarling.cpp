#include <sstream>
#include <iostream>
#include "FESolventSupplyStarling.h"
#include "FEMultiphasic.h" // need this for FESolutesMaterialPoint. can we move elsewhere?

// define the material parameters
BEGIN_PARAMETER_LIST(FESolventSupplyStarling, FESolventSupply)
	ADD_PARAMETER(m_kp, FE_PARAM_DOUBLE, "kp");
	ADD_PARAMETER(m_pv, FE_PARAM_DOUBLE, "pv");
	ADD_PARAMETER(m_qc[0], FE_PARAM_DOUBLE, "q1");
	ADD_PARAMETER(m_qc[1], FE_PARAM_DOUBLE, "q2");
	ADD_PARAMETER(m_qc[2], FE_PARAM_DOUBLE, "q3");
	ADD_PARAMETER(m_qc[3], FE_PARAM_DOUBLE, "q4");
	ADD_PARAMETER(m_qc[4], FE_PARAM_DOUBLE, "q5");
	ADD_PARAMETER(m_qc[5], FE_PARAM_DOUBLE, "q6");
	ADD_PARAMETER(m_qc[6], FE_PARAM_DOUBLE, "q7");
	ADD_PARAMETER(m_qc[7], FE_PARAM_DOUBLE, "q8");
	ADD_PARAMETER(m_cv[0], FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_cv[1], FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_cv[2], FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_cv[3], FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_cv[4], FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_cv[5], FE_PARAM_DOUBLE, "c6");
	ADD_PARAMETER(m_cv[6], FE_PARAM_DOUBLE, "c7");
	ADD_PARAMETER(m_cv[7], FE_PARAM_DOUBLE, "c8");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESolventSupplyStarling::FESolventSupplyStarling(FEModel* pfem) : FESolventSupply(pfem)
{
	m_kp = 0;
	m_pv = 0;
	for (int i=0; i<MAX_CDOFS; ++i) {
		m_qc[i] = 0;
		m_cv[i] = 0;
	}
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESolventSupplyStarling::Init()
{
}

//-----------------------------------------------------------------------------
//! Solvent supply
double FESolventSupplyStarling::Supply(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& mpt = *mp.ExtractData<FESolutesMaterialPoint>();

	// evaluate solvent supply from pressure drop
	double phiwhat = m_kp*(m_pv - ppt.m_p);
	
	// evaluate solvent supply from concentration drop
	if (&mpt) {
		int nsol = mpt.m_nsol;
		for (int isol=0; isol<nsol; ++isol) {
			phiwhat += m_qc[isol]*(m_cv[isol] - mpt.m_c[isol]);
		}
	}
	
	return phiwhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to strain
mat3ds FESolventSupplyStarling::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	mat3dd Phie(Supply(mp));
	
	return Phie;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to pressure
double FESolventSupplyStarling::Tangent_Supply_Pressure(FEMaterialPoint &mp)
{
	return -m_kp;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to concentration
double FESolventSupplyStarling::Tangent_Supply_Concentration(FEMaterialPoint &mp, const int isol)
{
	FESolutesMaterialPoint& mpt = *mp.ExtractData<FESolutesMaterialPoint>();
	if (isol < mpt.m_nsol) {
		return -m_qc[isol];
	}
	
	return 0;
}

