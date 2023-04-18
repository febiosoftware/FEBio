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
#include "FEModelComponent.h"
#include "FEMaterialPoint.h"
#include "FEModelParam.h"
#include "FEDomainList.h"
#include "FEDomainParameter.h"

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEDomain;
class DumpStream;

//-----------------------------------------------------------------------------
class FECORE_API FEMaterialBase : public FEModelComponent
{
public:
	FEMaterialBase(FEModel* fem);

	//! returns a pointer to a new material point object
	virtual FEMaterialPointData* CreateMaterialPointData();

	//! Update specialized material points at each iteration
	virtual void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp);

	// evaluate local coordinate system at material point
	virtual mat3d GetLocalCS(const FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types
//! From this class all other material classes are derived.

class FECORE_API FEMaterial : public FEMaterialBase
{
	FECORE_SUPER_CLASS(FEMATERIAL_ID)
	FECORE_BASE_CLASS(FEMaterial)

public:
	FEMaterial(FEModel* fem);
	virtual ~FEMaterial();

	//! performs initialization
	bool Init() override;

	//! get a domain parameter
	FEDomainParameter* FindDomainParameter(const std::string& paramName);

	// evaluate local coordinate system at material point
	mat3d GetLocalCS(const FEMaterialPoint& mp) override;

	// set the (local) material axis valuator
	void SetMaterialAxis(FEMat3dValuator* val);

protected:
	FEMat3dValuator*	m_Q;			//!< local material coordinate system

public:
	//! Assign a domain to this material
	void AddDomain(FEDomain* dom);

	//! get the domaint list
	FEDomainList& GetDomainList() { return m_domList; }

protected:
	void AddDomainParameter(FEDomainParameter* p);

private:
	FEDomainList	m_domList;		//!< list of domains that use this material

	std::vector<FEDomainParameter*>	m_param;	//!< list of domain variables

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Material properties are classes that can only be defined as properties of other materials
class FECORE_API FEMaterialProperty : public FEMaterialBase
{
	FECORE_SUPER_CLASS(FEMATERIALPROP_ID)

public:
	FEMaterialProperty(FEModel* fem);

	// evaluate local coordinate system at material point
	mat3d GetLocalCS(const FEMaterialPoint& mp) override;

	DECLARE_FECORE_CLASS();
};
