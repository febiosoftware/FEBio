/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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



#include "FEViscosity.h"
#include "FEFluidMaterialPoint.h"
#include "FEElasticFluid.h"
#include "FEFluid.h"

FEViscosity::FEViscosity(FEModel* pfem) : FEMaterialProperty(pfem)
{
    m_pIGisentropic = nullptr;
}

bool FEViscosity::Init()
{
    FEFluid* pfl = dynamic_cast<FEFluid*>(GetAncestor());
    if (pfl) {
        FEElasticFluid* pef = pfl->GetElastic();
        if (pef) {
            m_pIGisentropic = dynamic_cast<FEIdealGasIsentropic*>(pef);
        }
    }
    return FEMaterialProperty::Init();
}

//! tangent of normalized viscosity with respect to temperature
double FEViscosity::Tangent_NormalizedViscosity_Temperature(FEMaterialPoint& mp)
{
    // if the elastic response of the fluid is that of an isentropic ideal gas
    // then temperature depends on dilatation
    if (m_pIGisentropic) {
        FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
        double J = 1 + fp.m_ef;
        double gamma = m_pIGisentropic->m_gamma;
        double Tr = m_pIGisentropic->m_Tr;
        double T = Tr*pow(J,1-gamma);   // this is the absolute temperature
        double dJdT = J/(1-gamma)/T;
        double detahatdT = Tangent_NormalizedViscosity_Strain(mp)*dJdT;
        return detahatdT;
    }
    return 0.;
}

