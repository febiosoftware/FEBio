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
#include "FECentrifugalFluidBodyForce.h"

BEGIN_FECORE_CLASS(FECentrifugalFluidBodyForce, FEBodyForce);
	ADD_PARAMETER(w, "angular_speed")->setUnits(UNIT_ANGULAR_VELOCITY);
	ADD_PARAMETER(n, "rotation_axis");
	ADD_PARAMETER(c, "rotation_center")->setUnits(UNIT_LENGTH);
END_FECORE_CLASS();

FECentrifugalFluidBodyForce::FECentrifugalFluidBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	w = 0.0;
	n = vec3d(0,0,1);
	c = vec3d(0,0,0);
}

vec3d FECentrifugalFluidBodyForce::force(FEMaterialPoint& mp)
{
	mat3d K = stiffness(mp);
	return K*(mp.m_rt - c);
}

double FECentrifugalFluidBodyForce::divforce(FEMaterialPoint& mp)
{
    return -2*w*w;
}

mat3d FECentrifugalFluidBodyForce::stiffness(FEMaterialPoint& mp)
{ 
	return (mat3dd(1) - dyad(n))*(-w*w); 
}
