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
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
class FEJ2PlasticMaterialPoint : public FEElasticMaterialPoint
{
public:
	FEJ2PlasticMaterialPoint() {}

	FEMaterialPointData* Copy() 
	{
		return new FEJ2PlasticMaterialPoint(*this);
	}

	void Init()
	{
		// intialize data to zero
		e0.zero();
		e1.zero();
		sn.zero();
		b = false;
		Y1 = Y0;

		// don't forget to intialize the nested data
		FEMaterialPointData::Init();
	}

	void Update(const FETimeInfo& timeInfo)
	{
		e0 = e1;
		sn = m_s;
		Y0 = Y1;

		// don't forget to call the base class
		FEElasticMaterialPoint::Update(timeInfo);
	}

	void Serialize(DumpStream& ar)
	{
		FEElasticMaterialPoint::Serialize(ar);
		ar & e0 & e1 & sn;
		ar & Y0 & Y1 & b;
	}

public:
	mat3ds	e0, e1;		// strain at time n and n+1
	mat3ds	sn;			// stress at time n
	double	Y0, Y1;		// yield strenght at time n, n+1
	bool	b;			// plasticity flag
};

//-----------------------------------------------------------------------------
//! This class implements a simple von-Mises plasticity model with isotropic
//! hardening. 
class FEVonMisesPlasticity : public FESolidMaterial
{
public:
	FEVonMisesPlasticity(FEModel* pfem);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

	double	m_K;	//!< bulk modulus
	double	m_G;	//!< shear modulus
	double	m_Y;	//!< initial yield strength
	double	m_H;	//!< hardening modulus 

public:
	FEMaterialPointData* CreateMaterialPointData() override;

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! data initialization and checking
	bool Init() override;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
