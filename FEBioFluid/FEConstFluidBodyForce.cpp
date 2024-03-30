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
#include "FEConstFluidBodyForce.h"

BEGIN_FECORE_CLASS(FEConstFluidBodyForce, FEBodyForce);
    ADD_PARAMETER(m_force, "force")->setUnits(UNIT_SPECIFIC_FORCE);
END_FECORE_CLASS();

FEConstFluidBodyForce::FEConstFluidBodyForce(FEModel* pfem) : FEBodyForce(pfem) 
{ 
}

vec3d FEConstFluidBodyForce::force(FEMaterialPoint& pt) 
{ 
    return m_force(pt);
}

mat3d FEConstFluidBodyForce::stiffness(FEMaterialPoint& pt) 
{ 
	return mat3ds(0, 0, 0, 0, 0, 0); 
}
