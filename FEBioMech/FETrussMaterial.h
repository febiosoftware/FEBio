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

#include "FEElasticMaterial.h"
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
// Material point class for truss materials
class FETrussMaterialPoint : public FEElasticMaterialPoint
{
public:
	FEMaterialPointData* Copy()
	{
		FETrussMaterialPoint* pt = new FETrussMaterialPoint(*this);
		if (m_pNext) pt->m_pNext = m_pNext->Copy();
		return pt;
	}

	void Serialize(DumpStream& ar)
	{
		FEElasticMaterialPoint::Serialize(ar);
		ar & m_lam & m_tau;
	}

	void Init()
	{
		FEElasticMaterialPoint::Init();
		m_lam = 1;
		m_tau = 0;
	}

public:
	double	m_lam;	// stretch
	double	m_tau;	// Kirchoff stress
};

//-----------------------------------------------------------------------------
// Base class for truss element materials
class FEBIOMECH_API FETrussMaterial : public FEMaterial
{
public:
	FETrussMaterial(FEModel* pfem);
	~FETrussMaterial();

public:
	double	m_rho;	// density

public:
	//! calculate Kirchhoff stress of truss
	virtual double Stress(FEMaterialPoint& pt) = 0;

	//! calculate elastic tangent
	virtual double Tangent(FEMaterialPoint& pt) = 0;

	//! create material point data
	FEMaterialPointData* CreateMaterialPointData() override { return new FETrussMaterialPoint; }

	//! material density
	double Density();

	// declare the parameter list
	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FETrussMaterial);
};

//-----------------------------------------------------------------------------
class FELinearTrussMaterial : public FETrussMaterial
{
public:
	FELinearTrussMaterial(FEModel* fem);

	//! calculate Kirchhoff stress of truss
	double Stress(FEMaterialPoint& pt) override;

	//! calculate elastic tangent
	double Tangent(FEMaterialPoint& pt) override;

public:
	double	m_E;	// Elastic modulus
	double	m_v;	// Poisson's ratio

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

class FEBIOMECH_API FETrussStress : public FEDomainParameter
{
public:
	FETrussStress();
	FEParamValue value(FEMaterialPoint& mp) override;
};
