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
#include "febiomech_api.h"
#include <FECore/FEMaterial.h>
#include <FECore/FEMeshTopo.h>

//-----------------------------------------------------------------------------
// Virtual base class for nonlocal kernels
class FEBIOMECH_API FENonLocalKernel : public FEMaterialProperty
{
public:
    FENonLocalKernel(FEModel* pfem) : FEMaterialProperty(pfem) { m_R = 0; m_mult = 1; m_c = 1; }
    virtual ~FENonLocalKernel() {}
    
    //! kernel calculation
    virtual double Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt) = 0;
    
    //! kernel integral over infinite domain
    virtual double KernelIntegralInfinity() = 0;

public:
    double  m_R;    //!< interaction radius
    double  m_mult; //!< multiplier of interaction radius such that kernel becomes negligibly small
    double  m_c;    //!< normalization of kernel (depends on problem dimension)

    FECORE_BASE_CLASS(FENonLocalKernel);
};

//-----------------------------------------------------------------------------
// bell-shaped kernel
class FEBIOMECH_API FEKernelBell : public FENonLocalKernel
{
public:
    FEKernelBell(FEModel* pfem) : FENonLocalKernel(pfem) { m_mult = 1; }
    
    //! kernel calculation
    double Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt) override;
    
    //! kernel integral over infinite domain
    double KernelIntegralInfinity() override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// cone-shaped kernel
class FEBIOMECH_API FEKernelCone : public FENonLocalKernel
{
public:
    FEKernelCone(FEModel* pfem) : FENonLocalKernel(pfem) { m_mult = 1; }
    
    double Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt) override;
    
    //! kernel integral over infinite domain
    double KernelIntegralInfinity() override;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Gauss-shaped kernel
class FEBIOMECH_API FEKernelGauss : public FENonLocalKernel
{
public:
    FEKernelGauss(FEModel* pfem) : FENonLocalKernel(pfem) { m_mult = 5; }
    
    double Kernel(FEMaterialPoint& p0, FEMaterialPoint& pt) override;

    //! kernel integral over infinite domain
    double KernelIntegralInfinity() override;
    
    DECLARE_FECORE_CLASS();
};
