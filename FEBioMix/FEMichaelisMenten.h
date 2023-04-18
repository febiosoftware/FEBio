/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEChemicalReaction.h"

//-----------------------------------------------------------------------------
//! Forward chemical reaction following Michaelis-Menten kinetics.
//! The maximum uptake rate is given by the forward reaction rate
//! which is defined in the parent class FEChemicalReaction.
//! Optionally, a minimum concentration may be prescribed for the
//! reactant to trigger the reaction.

class FEBIOMIX_API FEMichaelisMenten : public FEChemicalReaction
{
public:
	//! constructor
	FEMichaelisMenten(FEModel* pfem);
	
	//! data initialization and checking
	bool Init() override;
	
	//! molar supply at material point
	double ReactionSupply(FEMaterialPoint& pt) override;
	
	//! tangent of molar supply with strain at material point
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) override;
	
	//! tangent of molar supply with effective pressure at material point
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) override;
	
	//! tangent of molar supply with effective concentration at material point
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) override;

public:
	double	m_Km;			//!< concentration at which half-maximum rate occurs
	int		m_Rid;			//!< local id of reactant
	int		m_Pid;			//!< local id of product
	bool	m_Rtype;		//!< flag for reactant type (solute = false, sbm = true)
	double	m_c0;			//!< minimum reactant concentration to trigger reaction
	
	// declare parameter list
	DECLARE_FECORE_CLASS();	
};

