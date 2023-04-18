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
#include <FECore/FEMaterial.h>
#include <FECore/DumpStream.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
class FEMembraneMaterialPoint : public FEMaterialPointData
{
public:
	FEMembraneMaterialPoint()
	{
		s[0] = s[1] = s[2] = 0;
	}

	FEMaterialPointData* Copy()
	{
		return new FEMembraneMaterialPoint(*this);
	}

	void Serialize(DumpStream& ar)
	{
		FEMaterialPointData::Serialize(ar);
		ar & g & s;
	}

public:
	// calculate membrane strain
	void strain(double e[3]);

public:
	double	g[6];	// deformation gradient
	double	s[3];	// in-plane PK2 stress
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEMembraneMaterial : public FEMaterial
{
public:
	FEMembraneMaterial(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEMembraneMaterial() {}

	//! create material point data
	FEMaterialPointData* CreateMaterialPointData() { return new FEMembraneMaterialPoint; }

public:
	//! calculate in-plane membrane stress
	virtual void Stress(FEMaterialPoint& mp, double s[3]) = 0;

	//! calculate in-plane membrane tangent
	virtual void Tangent(FEMaterialPoint& mp, double D[3][3]) = 0;

	FECORE_BASE_CLASS(FEMembraneMaterial)
};

//-----------------------------------------------------------------------------
// class for triangular membranes
class FEElasticMembrane : public FEMembraneMaterial
{
public:
	//! constructor
	FEElasticMembrane(FEModel* pfem) : FEMembraneMaterial(pfem) {}

public:
	//! calculate in-plane membrane stress
	void Stress(FEMaterialPoint& mp, double s[3]) override;

	//! calculate in-plane membrane tangent
	void Tangent(FEMaterialPoint& mp, double D[3][3]) override;

public:
	double	m_E;
	double	m_v;

	DECLARE_FECORE_CLASS();
};
