/*
 *  FECellGrowth.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 7/8/11.
 *  Copyright 2011 Columbia University. All rights reserved.
 *
 */

#include "stdafx.h"
#include "FECellGrowth.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FECellGrowth, FEElasticMaterial)
ADD_PARAMETER(m_phir, FE_PARAM_DOUBLE, "phir");
ADD_PARAMETER(m_cr, FE_PARAM_DOUBLE, "cr");
ADD_PARAMETER(m_ce, FE_PARAM_DOUBLE, "ce");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FECellGrowth::Init()
{
	FEElasticMaterial::Init();

	if (m_phir < 0) throw MaterialError("phir must be positive.");
	if (m_cr < 0) throw MaterialError("cr must be positive.");
	if (m_ce < 0) throw MaterialError("ce must be positive.");
	
	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	
	if (m_Rgas <= 0) throw MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) throw MaterialError("A positive absolute temperature T must be defined in Globals section");
	
}

//-----------------------------------------------------------------------------
mat3ds FECellGrowth::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate intracellular osmolarity relative to mixture volume in reference configuration
	double c = m_cr/(J-m_phir);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(c - m_ce);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FECellGrowth::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate intracellular osmolarity relative to mixture volume in reference configuration
	double c = m_cr/(J-m_phir);
		
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(c - m_ce);
	
	mat3dd I(1.0);	// Identity
	
	tens4ds I1 = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds C = I4*(2*p) - I1*(p-m_Rgas*m_Tabs*c*J/(J-m_phir));
	return C;
}

