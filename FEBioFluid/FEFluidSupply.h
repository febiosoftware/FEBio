/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2023 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FEMaterial.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Base class for solvent supply.
//! These materials need to define the supply and tangent supply functions.
//!
class FEBIOFLUID_API FEFluidSupply : public FEMaterialProperty
{
public:
    FEFluidSupply(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~FEFluidSupply(){}
	
	//! fluid supply
	virtual double Supply(FEMaterialPoint& pt) = 0;
	
	//! tangent of fluid supply with respect to strain
	virtual mat3d Tangent_Supply_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of fluid supply with respect to pressure
	virtual double Tangent_Supply_Dilatation(FEMaterialPoint& mp) = 0;
	
    //! tangent of fluid supply with respect to rate of deformation
    virtual mat3ds Tangent_Supply_RateOfDeformation(FEMaterialPoint& mp) = 0;
    
	//! tangent of fluid supply with respect to concentration
	virtual double Tangent_Supply_Concentration(FEMaterialPoint& mp, const int isol);

	FECORE_BASE_CLASS(FEFluidSupply)
};

