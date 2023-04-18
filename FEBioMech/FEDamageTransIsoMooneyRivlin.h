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
#include "FEUncoupledMaterial.h"

#ifdef WIN32
#define max(a,b) ((a)>(b)?(a):(b))
#endif

//-----------------------------------------------------------------------------
// We first define a material point that stores the damage variable.
class FETIMRDamageMaterialPoint : public FEMaterialPointData
{
public:
	FETIMRDamageMaterialPoint(FEMaterialPointData*pt);

	FEMaterialPointData* Copy();

	void Init();
	void Update(const FETimeInfo& timeInfo);

	void Serialize(DumpStream& ar);

public:
	// matrix
	double	m_MEtrial;			//!< trial strain at time t
	double	m_MEmax;			//!< max strain variable up to time t
	double	m_Dm;				//!< damage

	// fiber
	double	m_FEtrial;			//!< trial strain at time t
	double	m_FEmax;			//!< max strain variable up to time t
	double	m_Df;				//!< damage
};

//-----------------------------------------------------------------------------
class FEDamageTransIsoMooneyRivlin : public FEUncoupledMaterial
{
public:
	FEDamageTransIsoMooneyRivlin(FEModel* pfem);

public:
	// Mooney-Rivlin parameters
	double	m_c1;	//!< Mooney-Rivlin coefficient C1
	double	m_c2;	//!< Mooney-Rivlin coefficient C2

	// fiber parameters
	double	m_c3;
	double	m_c4;

	// Matrix damage parameters
	double	m_Mbeta;		//!< damage parameter beta
	double	m_Msmin;		//!< damage parameter psi-min
	double	m_Msmax;		//!< damage parameter psi-max

	// Fiber damage parameters
	double	m_Fbeta;
	double	m_Fsmin;
	double	m_Fsmax;

public:
	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! damage
    double Damage(FEMaterialPoint& pt);
    
protected:
	mat3ds MatrixStress(FEMaterialPoint& mp);
	mat3ds FiberStress (FEMaterialPoint& mp);
	tens4ds MatrixTangent(FEMaterialPoint& pt);
	tens4ds FiberTangent (FEMaterialPoint& pt);
	double MatrixStrainEnergyDensity(FEMaterialPoint& pt);
	double FiberStrainEnergyDensity (FEMaterialPoint& pt);

protected:
	// calculate damage reduction factor for matrix
	double MatrixDamage(FEMaterialPoint& pt);

	// calculate damage reduction factor for fibers
	double FiberDamage(FEMaterialPoint& pt);

	double MatrixDamageDerive(FEMaterialPoint& pt);
	double FiberDamageDerive(FEMaterialPoint& pt);

public:

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
