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
#include "FEFluidSupply.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a fluid supply following
// Starling's equation
class FEBIOFLUID_API FEFluidSupplyStarling :	public FEFluidSupply
{
public:
	//! constructor
    FEFluidSupplyStarling(FEModel* pfem);
	
	//! Solute supply
	double Supply(FEMaterialPoint& pt) override;
	
	//! Tangent of supply with respect to solid strain
	mat3d Tangent_Supply_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of supply with respect to pressure
	double Tangent_Supply_Dilatation(FEMaterialPoint& mp) override;
	
    //! tangent of fluid supply with respect to rate of deformation
    mat3ds Tangent_Supply_RateOfDeformation(FEMaterialPoint& mp) override { return mat3ds(0); }

    //! Initialization
    bool Init() override { return FEFluidSupply::Init(); }
	
   
public:
	FEParamDouble		m_kp;				//!< coefficient of pressure drop
    FEParamDouble		m_pv;				//!< prescribed (e.g., vascular) pressure
	
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
