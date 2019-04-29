/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEFiberDensityDistribution.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
// define the ellipsoidal fiber density distributionmaterial parameters
BEGIN_FECORE_CLASS(FEEllipsodialFiberDensityDistribution, FEFiberDensityDistribution)
	ADD_PARAMETER(m_spa , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "spa" );
END_FECORE_CLASS();

double FEEllipsodialFiberDensityDistribution::FiberDensity(const vec3d& n0)
{
    double R = 1.0/sqrt(SQR(n0.x/m_spa[0])+SQR(n0.y/m_spa[1])+SQR(n0.z/m_spa[2]));
    return R;
}

//-----------------------------------------------------------------------------
// define the 3d von Mises fiber density distribution material parameters
BEGIN_FECORE_CLASS(FEVonMises3DFiberDensityDistribution, FEMaterial)
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "b" );
END_FECORE_CLASS();

double FEVonMises3DFiberDensityDistribution::FiberDensity(const vec3d& n0)
{
    // The local x-direction is the principal fiber bundle direction
    // The x-component of n0 is cos(phi)
    double R = exp(m_b*(2*SQR(n0.x)-1));
    return R;
}

//-----------------------------------------------------------------------------
// define the 3d 2-fiber family axisymmetric von Mises fiber density distribution material parameters
BEGIN_FECORE_CLASS(FEVonMises3DTwoFDDAxisymmetric, FEMaterial)
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "b" );
	ADD_PARAMETER(m_c, FE_RANGE_CLOSED(0, 1), "cosg" );
END_FECORE_CLASS();

double FEVonMises3DTwoFDDAxisymmetric::FiberDensity(const vec3d& n0)
{
    // The local x-direction is the principal fiber bundle direction
    // The x-component of n0 is cos(phi)
    double cphi = n0.x; double sphi = sqrt(1-cphi*cphi);
    double sing = sqrt(1-m_c*m_c);
    double cp = cphi*m_c - sphi*sing;
    double cm = cphi*m_c + sphi*sing;
    double R = exp(m_b*(2*SQR(cp)-1)) + exp(m_b*(2*SQR(cm)-1));
    return R;
}

//-----------------------------------------------------------------------------
// define the ellipsoidal fiber density distributionmaterial parameters
BEGIN_FECORE_CLASS(FEEllipticalFiberDensityDistribution, FEMaterial)
	ADD_PARAMETER(m_spa[0], FE_RANGE_GREATER_OR_EQUAL(0.0), "spa1" );
	ADD_PARAMETER(m_spa[1], FE_RANGE_GREATER_OR_EQUAL(0.0), "spa2" );
END_FECORE_CLASS();

double FEEllipticalFiberDensityDistribution::FiberDensity(const vec3d& n0)
{
    // 2d fibers lie in the local x-y plane
    // n0.x = cos(theta) and n0.y = sin(theta)
    double R = 1.0/sqrt(SQR(n0.x/m_spa[0])+SQR(n0.y/m_spa[1]));
    return R;
}

//-----------------------------------------------------------------------------
// define the 2d von Mises fiber density distribution material parameters
BEGIN_FECORE_CLASS(FEVonMises2DFiberDensityDistribution, FEMaterial)
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "b" );
END_FECORE_CLASS();

double FEVonMises2DFiberDensityDistribution::FiberDensity(const vec3d& n0)
{
    // The fiber bundle is in the x-y plane and
    // the local x-direction is the principal fiber bundle direction
    // The x-component of n0 is cos(theta)
    double R = exp(m_b*(2*SQR(n0.x)-1));
    return R;
}
