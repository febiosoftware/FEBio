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
#include "FEUncoupledFiberExpLinear.h"

//-----------------------------------------------------------------------------
class FEMRVonMisesMaterialPoint : public FEMaterialPointData
{
public:
	FEMRVonMisesMaterialPoint(FEMaterialPointData* mp = nullptr);

	FEMaterialPointData* Copy() override;

	void Serialize(DumpStream& ar) override;

public:
	double	m_kf;
	double	m_tp;
};

//-----------------------------------------------------------------------------
//! Transversely Isotropic Multiple material

//! This material has an isotopric Multiple basis and single preferred
//! fiber direction.

class FEMRVonMisesFibers: public FEUncoupledMaterial
{
public:
	FEMRVonMisesFibers (FEModel* pfem);

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2
	double	kf;	//!< Fiber Concentration Factor
	double	tp;	//!< Preferred Fiber Orientation IN RADIANS
	int	gipt;	//!< Gauss Integration Points
	int	vmc;	//!< Choice of von Mises distribution 
	double var_n;	//!< Exponent for the constrained von Mises distribution

public:
	//! calculate stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! create material point data
	FEMaterialPointData* CreateMaterialPointData() override;

	// declare parameter list
	DECLARE_FECORE_CLASS();

protected:
	FEFiberExpLinearUC	m_fib;
};
