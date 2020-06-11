/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
// Material point class for truss materials
class FETrussMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FETrussMaterialPoint* pt = new FETrussMaterialPoint(*this);
		if (m_pNext) pt->m_pNext = m_pNext->Copy();
		return pt;
	}

	void Serialize(DumpStream& ar)
	{
		FEMaterialPoint::Serialize(ar);
		ar & m_l & m_tau;
	}

	void Init()
	{
		FEMaterialPoint::Init();
		m_l = 1;
		m_tau = 0;
	}

public:
	double	m_l;	// strech
	double	m_tau;	// Kirchoff stress
};

//-----------------------------------------------------------------------------
// Base class for truss element materials
class FETrussMaterial : public FEMaterial
{
public:
	FETrussMaterial(FEModel* pfem) : FEMaterial(pfem) {}
	~FETrussMaterial(){}

public:
	double	m_E;	// Elastic modulus

public:
	//! calculate Kirchhoff stress of truss
	virtual double Stress(FEMaterialPoint& pt);

	//! calculate elastic tangent
	virtual double Tangent(FEMaterialPoint& pt);

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() override { return new FETrussMaterialPoint; }

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
