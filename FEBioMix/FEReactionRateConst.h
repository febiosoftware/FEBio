#pragma once
#include "FEMultiphasic.h"

class FEBIOMIX_API FEReactionRateConst : public FEReactionRate
{
public:
	//! constructor
	FEReactionRateConst(FEModel* pfem) : FEReactionRate(pfem) { m_k = 0; }
	
	//! reaction rate at material point
	double ReactionRate(FEMaterialPoint& pt) override { return m_k; }
	
	//! tangent of reaction rate with strain at material point
	mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) override { return mat3ds(0); }
	
	//! tangent of reaction rate with effective fluid pressure at material point
	double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) override {return 0; }

public:
	double	m_k;		//!< reaction rate
	
	DECLARE_FECORE_CLASS();	
};
