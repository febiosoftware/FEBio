/*
 *  FEMichaelisMenten.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/8/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

//-----------------------------------------------------------------------------
//! Forward chemical reaction following Michaelis-Menten kinetics.
//! The maximum uptake rate is given by the forward reaction rate
//! which is defined in the parent class FEChemicalReaction.
//! Optionally, a minimum concentration may be prescribed for the
//! reactant to trigger the reaction.

#include "FEMultiphasic.h"

class FEMichaelisMenten : public FEChemicalReaction
{
public:
	//! constructor
	FEMichaelisMenten(FEModel* pfem) : FEChemicalReaction(pfem) {m_Rid = m_Pid = -1; m_Km = m_c0 = 0; m_Rtype = false; }
	
	//! data initialization and checking
	void Init();
	
	//! molar supply at material point
	double ReactionSupply(FEMaterialPoint& pt);
	
	//! tangent of molar supply with strain at material point
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective pressure at material point
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective concentration at material point
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);

public:
	double	m_Km;			//!< concentration at which half-maximum rate occurs
	int		m_Rid;			//!< local id of reactant
	int		m_Pid;			//!< local id of product
	bool	m_Rtype;		//!< flag for reactant type (solute = false, sbm = true)
	double	m_c0;			//!< minimum reactant concentration to trigger reaction
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();	
};

