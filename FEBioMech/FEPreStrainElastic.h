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

//-----------------------------------------------------------------------------
//! Class for stroing pre-strain gradient data
//! This material point class stores the initial pre-strain "guess" and
//! the correction term. The total prestrain gradient is the product of these two.
class FEPrestrainMaterialPoint : public FEMaterialPointData
{
public:
	//! constructor
	FEPrestrainMaterialPoint(FEMaterialPointData* mp);

	//! copy
	FEMaterialPointData* Copy();

	//! initialization
	void Init(bool bflag);

	//! Serialization
	void Serialize(DumpStream& dmp);

public:
	//! get the initial prestrain gradient
	const mat3d& initialPrestrain() const { return F0; }

	//! set the prestrain gradient
	void setInitialPrestrain(const mat3d& F) { F0 = F; }

	//! get the prestrain correction
	const mat3d& PrestrainCorrection() const { return Fc; }

	//! set the prestrain correction
	void setPrestrainCorrection(const mat3d& F) { Fc = F; }

	//! get the total prestrain
	mat3d prestrain() const { return Fc*F0; }

protected:
	mat3d	F0;	//!< initial guess of pre-strain value
	mat3d	Fc;	//!< "correction" used to premultiply the prestrain gradient
};

//-----------------------------------------------------------------------------
//! Base class for algorithms that will be used to calculate a pre-strain gradient.
//! This is used by the FEPrestrainElastic class to calculate the initial pre-strain
//! gradient.
class FEBIOMECH_API FEPrestrainGradient : public FEMaterialProperty
{
public:
	FEPrestrainGradient(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~FEPrestrainGradient(){}

	// evaluate the pre-strain deformation gradient
	virtual mat3d Prestrain(FEMaterialPoint& mp) = 0;

	// initialize the pre-strain gradient based on a deformation gradient
	// This is used by the pre-strain initial condition
	virtual void Initialize(const mat3d& F, FEMaterialPoint& mp) = 0;

	FECORE_BASE_CLASS(FEPrestrainGradient)
};

//-----------------------------------------------------------------------------
class FEPrestrainMaterial
{
public:
	virtual FEPrestrainGradient* PrestrainGradientProperty() = 0;

	virtual FEElasticMaterial* GetElasticMaterial() = 0;
};

//-----------------------------------------------------------------------------
//! This material applies a user-defined prestrain deformation gradient
//! before evaluating the stress and tangents. 
class FEPrestrainElastic : public FEElasticMaterial, public FEPrestrainMaterial
{
public:
	//! constructor
	FEPrestrainElastic(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

	//! return the pre-strain gradient property
	FEPrestrainGradient* PrestrainGradientProperty() override { return m_Fp; }

	//! return the elastic material
	FEElasticMaterial* GetElasticMaterial() override { return m_mat; }

	// evaluate density in (pre-strained) reference configuration
	double Density(FEMaterialPoint& mp) override;

public:
	//! Cauchy stress 
	mat3ds Stress(FEMaterialPoint& mp) override;

	//! spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

protected:
	//! calculate prestrain deformation gradient
	mat3d PrestrainGradient(FEMaterialPoint& mp);

private:
	FEElasticMaterial*		m_mat;	//!< elastic base material
	FEPrestrainGradient*	m_Fp;	//!< pre-strain gradient

	DECLARE_FECORE_CLASS();
};
