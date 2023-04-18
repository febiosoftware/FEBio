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



#include "stdafx.h"
#include "FEGenericBodyForce.h"
#include "FEElasticMaterial.h"

BEGIN_FECORE_CLASS(FEGenericBodyForce, FEBodyForce);
	ADD_PARAMETER(m_force, "force")->setUnits(UNIT_SPECIFIC_FORCE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEGenericBodyForce::FEGenericBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
vec3d FEGenericBodyForce::force(FEMaterialPoint &mp)
{
	return m_force(mp);
}

//-----------------------------------------------------------------------------
mat3ds FEGenericBodyForce::stiffness(FEMaterialPoint& pt)
{
	return mat3ds(0, 0, 0, 0, 0, 0);
}

//=============================================================================
BEGIN_FECORE_CLASS(FEConstBodyForceOld, FEBodyForce);
	ADD_PARAMETER(m_f.x, "x");
	ADD_PARAMETER(m_f.y, "y");
	ADD_PARAMETER(m_f.z, "z");
END_FECORE_CLASS();


//=============================================================================
BEGIN_FECORE_CLASS(FENonConstBodyForceOld, FEBodyForce);
	ADD_PARAMETER(m_f[0], "x");
	ADD_PARAMETER(m_f[1], "y");
	ADD_PARAMETER(m_f[2], "z");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENonConstBodyForceOld::FENonConstBodyForceOld(FEModel* pfem) : FEBodyForce(pfem)
{
}

//-----------------------------------------------------------------------------
vec3d FENonConstBodyForceOld::force(FEMaterialPoint& pt)
{
	vec3d F;
	F.x = m_f[0](pt);
	F.y = m_f[1](pt);
	F.z = m_f[2](pt);
	return F;
}
