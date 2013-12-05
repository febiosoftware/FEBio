#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a solvent supply following
// Starling's equation

class FESolventSupplyStarling :	public FESolventSupply
{
public:
	//! constructor
	FESolventSupplyStarling(FEModel* pfem);
	
	//! Solute supply
	double Supply(FEMaterialPoint& pt);
	
	//! Tangent of supply with respect to strain
	mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp);
	
	//! Tangent of supply with respect to pressure
	double Tangent_Supply_Pressure(FEMaterialPoint& mp);
	
	//! Tangent of supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp, const int isol);
	
	//! data initialization and checking
	void Init();
	
public:
	double		m_kp;				//!< coefficient of pressure drop
	double		m_pv;				//!< prescribed (e.g., vascular) pressure
	double		m_qc[MAX_CDOFS];	//!< coefficients of concentration drops
	double		m_cv[MAX_CDOFS];	//!< prescribed (e.g., vascular) concentrations
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
