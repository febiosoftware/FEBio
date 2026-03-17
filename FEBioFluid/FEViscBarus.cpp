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



#include "FEViscBarus.h"
#include "FEFluidMaterialPoint.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEViscBarus, FEViscosity)
    ADD_PARAMETER(m_alpha, "alpha")->setUnits("1/P");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEViscBarus::FEViscBarus(FEModel* pfem) : FEViscosity(pfem)
{
    m_pFluid = nullptr;
    m_alpha = 0;
}

bool FEViscBarus::Init()
{
    m_pFluid = dynamic_cast<FEFluidMaterial*>(GetAncestor());
    if (m_pFluid == nullptr) return false;
    return true;
}

//-----------------------------------------------------------------------------
//! normalize viscosity
double FEViscBarus::NormalizedViscosity(FEMaterialPoint& mp)
{
    double p = m_pFluid->Pressure(mp);
    return exp(m_alpha*p);
}

//-----------------------------------------------------------------------------
//! tangent of normalized viscosity with respect to temperature
double FEViscBarus::Tangent_NormalizedViscosity_Temperature(FEMaterialPoint& mp)
{
    return 0.;
}

//-----------------------------------------------------------------------------
//! tangent of normalized viscosity with respect to volumetric strain (or J)
double FEViscBarus::Tangent_NormalizedViscosity_Strain(FEMaterialPoint& mp)
{
    double p = m_pFluid->Pressure(mp);
    double dpdJ = m_pFluid->Tangent_Pressure_Strain(mp);
    return m_alpha*dpdJ*exp(m_alpha*p);
}
