/*
 *  FEReactionRateConst.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/6/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEMultiphasic.h"

class FEReactionRateConst : public FEReactionRate
{
public:
	//! constructor
	FEReactionRateConst() { m_k = 0; }
	
	//! data initialization and checking
	void Init();
	
	//! reaction rate at material point
	double ReactionRate(FEMaterialPoint& pt) { return m_k; }
	
	//! tangent of reaction rate with strain at material point
	mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) { return mat3dd(0); }
	
	//! tangent of reaction rate with effective fluid pressure at material point
	double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) {return 0; }
	
public:
	double	m_k;		//!< reaction rate
	
	// declare as registered
	DECLARE_REGISTERED(FEReactionRateConst);
	
	DECLARE_PARAMETER_LIST();	
};