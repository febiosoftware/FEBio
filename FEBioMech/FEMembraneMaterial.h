/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

//-----------------------------------------------------------------------------
class FEMembraneMaterialPoint : public FEMaterialPoint
{
public:
	FEMembraneMaterialPoint()
	{
		s[0] = s[1] = s[2] = 0;
	}

	FEMaterialPoint* Copy()
	{
		return new FEMembraneMaterialPoint();
	}

	void Serialize(DumpStream& ar)
	{
		FEMaterialPoint::Serialize(ar);
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
class FEMembraneMaterial : public FEMaterial
{
public:
	FEMembraneMaterial(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEMembraneMaterial() {}

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() { return new FEMembraneMaterialPoint; }

public:
	//! calculate in-plane membrane stress
	virtual void Stress(FEMaterialPoint& mp, double s[3]) = 0;

	//! calculate in-plane membrane tangent
	virtual void Tangent(FEMaterialPoint& mp, double D[3][3]) = 0;
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
