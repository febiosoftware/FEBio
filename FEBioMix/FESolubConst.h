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
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute solubility

class FEBIOMIX_API FESolubConst : public FESoluteSolubility
{
public:
	//! constructor
	FESolubConst(FEModel* pfem);
	
	//! solubility
	double Solubility(FEMaterialPoint& pt) override;
	
	//! Tangent of solubility with respect to strain
	double Tangent_Solubility_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of solubility with respect to concentration
	double Tangent_Solubility_Concentration(FEMaterialPoint& mp, const int isol) override;
	
	//! Cross derivative of solubility with respect to strain and concentration
	double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp, const int isol) override;
	
	//! Second derivative of solubility with respect to strain
	double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp) override;

	//! Second derivative of solubility with respect to concentration
	double Tangent_Solubility_Concentration_Concentration(FEMaterialPoint& mp, const int isol, const int jsol) override;

public:
	double	m_solub;			//!< solubility
	
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
