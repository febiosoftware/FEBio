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
#include "FEThermoViscousFluid.h"
#include "FEThermoFluid.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
//! Constructor.
FEThermoViscousFluid::FEThermoViscousFluid(FEModel* pfem) : FEViscousFluid(pfem)
{
    for (int i=0; i<MAX_PROPS; ++i) m_prop[i] = nullptr;
    m_Tr = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEThermoViscousFluid::Init()
{
    m_Tr = GetGlobalConstant("T");
    
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined for thermo-viscous fluids in Globals section"); return false; }
    
    return FEViscousFluid::Init();
}

//-----------------------------------------------------------------------------
void FEThermoViscousFluid::Serialize(DumpStream& ar)
{
    FEViscousFluid::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_Tr;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity
double FEThermoViscousFluid::ThermoProp(FEMaterialPoint& mp, const int iprop)
{
    double scf = 1;
    if (m_prop[iprop]) {
        FEThermoFluidMaterialPoint* tf = mp.ExtractData<FEThermoFluidMaterialPoint>();
        if (tf) {
            double That = (tf->m_T+m_Tr)/m_Tr;
            scf = m_prop[iprop]->value(That);
        }
    }
    return scf;
}

//-----------------------------------------------------------------------------
//! dynamic shear viscosity tangent w.r.t. temperature
double FEThermoViscousFluid::TangentThermoProp(FEMaterialPoint& mp, const int iprop)
{
    double scf = 1;
    if (m_prop[iprop]) {
        FEThermoFluidMaterialPoint* tf = mp.ExtractData<FEThermoFluidMaterialPoint>();
        if (tf) {
            double That = (tf->m_T+m_Tr)/m_Tr;
            scf = m_prop[iprop]->derive(That)/m_Tr;
        }
    }
    return scf;
}
