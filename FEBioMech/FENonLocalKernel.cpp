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
#include "FENonLocalKernel.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <algorithm>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEKernelBell, FENonLocalKernel)
    ADD_PARAMETER(m_R, FE_RANGE_GREATER_OR_EQUAL(0.), "R")->setUnits(UNIT_LENGTH)->setLongName("interaction radius");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
double FEKernelBell::Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt)
{
    double r = (p0.m_r0 - pt.m_r0).unit();
    double k = 0;
    if (r <= m_R*m_mult) 
	{
		double e = 1 - (r*r) / (m_R*m_R);
        k = e*e*m_c;
    }
    
    return k;
}

//! kernel integral over infinite domain
double FEKernelBell::KernelIntegralInfinity()
{
    return 105.0/(32.0*PI*m_R*m_R*m_R);
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEKernelCone, FENonLocalKernel)
    ADD_PARAMETER(m_R, FE_RANGE_GREATER_OR_EQUAL(0.), "R")->setUnits(UNIT_LENGTH)->setLongName("interaction radius");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
double FEKernelCone::Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt)
{
    double r = (p0.m_r0 - pt.m_r0).unit();
    double k = 0;
    if (r <= m_R*m_mult) {
        k = (1-r/m_R)*m_c;
    }
    
    return k;
}

//! kernel integral over infinite domain
double FEKernelCone::KernelIntegralInfinity()
{
    return 3/(PI*m_R*m_R*m_R);
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEKernelGauss, FENonLocalKernel)
    ADD_PARAMETER(m_R, FE_RANGE_GREATER_OR_EQUAL(0.), "R")->setUnits(UNIT_LENGTH)->setLongName("interaction radius");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
double FEKernelGauss::Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt)
{
    double r = (p0.m_r0 - pt.m_r0).unit();
    double k = 0;
    if (r <= m_R*m_mult) {
		double e = r / m_R;
        k = exp(-e*e*0.5)*m_c;
    }
    
    return k;
}
//! kernel integral over infinite domain
double FEKernelGauss::KernelIntegralInfinity()
{
    return pow(m_R*sqrt(2*PI),-3);
}

