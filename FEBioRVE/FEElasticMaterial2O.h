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
#include "FEBioMech/FEElasticMaterial.h"
#include <FECore/tens3d.h>
#include <FECore/tens4d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>

//-----------------------------------------------------------------------------
// second order continuum elastic material point
class FEElasticMaterialPoint2O : public FEMaterialPointData
{
public:
	//! constructor
	FEElasticMaterialPoint2O(FEMaterialPointData* pt);

	//! create a shallow copy
	FEMaterialPointData* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);

public:
	mat3d	m_PK1;			//!< PK1 stress

	tens3drs   m_G;			//!< gradient of deformation gradient
	tens3drs   m_Q;			//!< higher-order stress tensor
};

//-----------------------------------------------------------------------------
//! This is the base class for second-order continuum elastic material
class FEElasticMaterial2O : public FEElasticMaterial
{
public:
	//! constructor
	FEElasticMaterial2O(FEModel* fem);

public: // these functions must be implemented by derived classes

	//! Calculate PK1 stress and higher order stress Q
	virtual void Stress(FEMaterialPoint& mp, mat3d& P, tens3drs& Q) = 0;

	//! Calculate material tangents
	//! C = dP/dF
	//! L = dP/dG
	//! H = dQ/dF
	//! J = dQ/dG
	virtual void Tangent(FEMaterialPoint& mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J) = 0;

	//! create material point data
	FEMaterialPointData* CreateMaterialPointData() override;

public:
	double			m_beta;			//!< beta parameter for DG

	// flags for evaluating effects of stiffness contributions
	// TODO: I'll probably delete this when all bugs are found
	bool	m_bKDG1;
	bool	m_bKDG2;
	bool	m_bKDG3;
	bool	m_buseJ0;

private:
	// We don't need these functions
	// TODO: Perhaps we should derive this class directly from FEMaterial so
	//       that these functions are not inherited. 
	mat3ds Stress(FEMaterialPoint& pt) override;
	tens4ds Tangent(FEMaterialPoint& pt) override;

	DECLARE_FECORE_CLASS();
};
