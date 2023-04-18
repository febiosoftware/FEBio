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
#include "FEElasticFiberMaterialUC.h"
#include "FEFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Uncoupled formulation of the fiber-exp-linear material for use with uncoupled
//! solid mixtures.
class FEFiberExpLinearUC : public FEFiberMaterialUncoupled
{
public:
	//! Constructor
	FEFiberExpLinearUC(FEModel* pfem);

	//! calculate deviatoric stress at material point
	mat3ds DevFiberStress(FEMaterialPoint& pt, const vec3d& n0) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevFiberTangent(FEMaterialPoint& pt, const vec3d& n0) override;

	//! calculate deviatoric strain energy density at material point
	double DevFiberStrainEnergyDensity(FEMaterialPoint& pt, const vec3d& n0) override;

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers
	double	m_lam1;		//!< fiber stretch for straightened fibers

	DECLARE_FECORE_CLASS();
};

class FEUncoupledFiberExpLinear : public FEElasticFiberMaterialUC_T<FEFiberExpLinearUC>
{
public: 
	FEUncoupledFiberExpLinear(FEModel* fem) : FEElasticFiberMaterialUC_T<FEFiberExpLinearUC>(fem) {}
	DECLARE_FECORE_CLASS();
};
