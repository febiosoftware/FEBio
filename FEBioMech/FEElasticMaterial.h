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
#include "FESolidMaterial.h"
#include "FEElasticMaterialPoint.h"

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEBIOMECH_API FEElasticMaterial : public FESolidMaterial
{
public:
	//! constructor 
	FEElasticMaterial(FEModel* pfem);

	//! destructor
	~FEElasticMaterial();

	//! create material point data for this material
	FEMaterialPointData* CreateMaterialPointData() override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
    // get the elastic material
    virtual FEElasticMaterial* GetElasticMaterial() { return this; }

public:
	//! evaluates approximation to Cauchy stress using forward difference
	mat3ds SecantStress(FEMaterialPoint& pt, bool PK2 = false) override;

public:
    virtual double StrongBondSED(FEMaterialPoint& pt) { return StrainEnergyDensity(pt); }
    virtual double WeakBondSED(FEMaterialPoint& pt) { return 0; }

protected:
//	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FEElasticMaterial);
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEElasticStress : public FEDomainParameter
{
public:
	FEElasticStress();
	FEParamValue value(FEMaterialPoint& mp) override;
};
