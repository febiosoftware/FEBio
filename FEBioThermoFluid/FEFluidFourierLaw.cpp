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
#include "FEFluidFourierLaw.h"
#include "FEThermalCondConst.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFluidFourierLaw, FEFluidHeatFlux)
    ADD_PARAMETER(m_Kr, FE_RANGE_GREATER_OR_EQUAL(0.0), "Kr")->setUnits(UNIT_THERMAL_CONDUCTIVITY)->setLongName("referential thermal conductivity");

// Optionally add strain-dependent normalized viscosity relations
    ADD_PROPERTY(m_Khat, "Khat", FEProperty::Optional)->SetLongName("normalized thermal conductivity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEFluidFourierLaw::FEFluidFourierLaw(FEModel* pfem) : FEFluidHeatFlux(pfem)
{
    m_Khat = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEFluidFourierLaw::Init()
{
    FEModel* pfem = GetFEModel();
    if (m_Khat == nullptr) {
        m_Khat = new FEThermalCondConst(pfem);
    }
    m_Khat->Init();

    return FEFluidHeatFlux::Init();
}

//-----------------------------------------------------------------------------
//! viscous stress
vec3d FEFluidFourierLaw::HeatFlux(FEMaterialPoint& pt)
{
    FEThermoFluidMaterialPoint& tf = *pt.ExtractData<FEThermoFluidMaterialPoint>();
    
    double K = Conductivity(pt);
    
    vec3d q = -tf.m_gradT*K;
        
    return q;
}
