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
#include "FEBodyForce.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! This class defines a non-homogeneous body force, i.e. the force can depend
//! on the reference position
class FEGenericBodyForce : public FEBodyForce
{
public:
	//! constructor
	FEGenericBodyForce(FEModel* pfem);

	//! evaluate the body force
	vec3d force(FEMaterialPoint& pt) override;

	//! stiffness 
	mat3ds stiffness(FEMaterialPoint& pt) override;

public:
	FEParamVec3 m_force;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! This class defines a deformation-independent constant force (e.g. gravity)
// NOTE: This class is obsolete!
class FEConstBodyForceOld : public FEBodyForce
{
public:
	FEConstBodyForceOld(FEModel* pfem) : FEBodyForce(pfem) { m_f = vec3d(0, 0, 0); }
	vec3d force(FEMaterialPoint& pt) override { return m_f; }
	mat3ds stiffness(FEMaterialPoint& pt) override { return mat3ds(0, 0, 0, 0, 0, 0); }

protected:
	vec3d	m_f;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Class that implements old body force
// NOTE: This class is obsolete!
class FENonConstBodyForceOld : public FEBodyForce
{
public:
	FENonConstBodyForceOld(FEModel* fem);
	vec3d force(FEMaterialPoint& pt) override;
	mat3ds stiffness(FEMaterialPoint& pt) override { return mat3ds(0, 0, 0, 0, 0, 0); }

private:
	FEParamDouble	m_f[3];

	DECLARE_FECORE_CLASS();
};
