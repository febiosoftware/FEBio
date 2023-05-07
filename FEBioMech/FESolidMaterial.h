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
#include <FECore/tens4d.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FEBIOMECH_API FESolidMaterial : public FEMaterial
{
public:
	//! constructor
	FESolidMaterial(FEModel* pfem);

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate the 2nd Piola-Kirchhoff stress at material point
	virtual mat3ds PK2Stress(FEMaterialPoint& pt, const mat3ds E);

	//! calculate material tangent stiffness at material point
	virtual tens4dmm MaterialTangent(FEMaterialPoint& pt, const mat3ds E);

    //! calculate secant tangent stiffness at material point
    virtual tens4dmm SecantTangent(FEMaterialPoint& pt, bool mat = false);

	//! return the material density
	void SetDensity(const double d);

	//! evaluate density
	virtual double Density(FEMaterialPoint& pt);

	//! Is this a rigid material or not
	virtual bool IsRigid() const { return false; }

	tens4dmm SolidTangent(FEMaterialPoint& pt);

	virtual mat3ds SecantStress(FEMaterialPoint& pt, bool PK2 = false);
	virtual bool UseSecantTangent() { return false; }

protected:
	FEParamDouble	m_density;	//!< material density

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FESolidMaterial)
};
