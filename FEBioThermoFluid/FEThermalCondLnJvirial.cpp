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



#include "FEThermalCondLnJvirial.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FEBioFluid/FEFluidMaterialPoint.h>
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEThermalCondLnJvirial, FEFluidThermalConductivity)
    ADD_PROPERTY(m_coef[0] , "c0"   )->SetLongName("0th normalized conductivity coefficient");
    ADD_PROPERTY(m_coef[1] , "c1" , FEProperty::Optional)->SetLongName("1st normalized conductivity coefficient");
    ADD_PROPERTY(m_coef[2] , "c2" , FEProperty::Optional)->SetLongName("2nd normalized conductivity coefficient");
    ADD_PROPERTY(m_coef[3] , "c3" , FEProperty::Optional)->SetLongName("3rd normalized conductivity coefficient");
    ADD_PROPERTY(m_coef[4] , "c4" , FEProperty::Optional)->SetLongName("4th normalized conductivity coefficient");
    ADD_PROPERTY(m_coef[5] , "c5" , FEProperty::Optional)->SetLongName("5th normalized conductivity coefficient");
    ADD_PROPERTY(m_coef[6] , "c6" , FEProperty::Optional)->SetLongName("6th normalized conductivity coefficient");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEThermalCondLnJvirial::FEThermalCondLnJvirial(FEModel* pfem) : FEFluidThermalConductivity(pfem)
{
    for (int k=0; k<MAX_COEF; ++k) m_coef[k] = nullptr;
}

bool FEThermalCondLnJvirial::Init()
{
    m_Tr = GetGlobalConstant("T");
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }

    m_ncoef = 0;
    for (int k=0; k<MAX_COEF; ++k) {
        if (m_coef[k]) {
            m_coef[k]->Init();
            ++m_ncoef;
        }
    }
    if (m_ncoef < 1) { feLogError("At least one virial coefficient should be provided for the viscosity"); return false; }
    
    return true;
}

//-----------------------------------------------------------------------------
//! normalize viscosity
double FEThermalCondLnJvirial::NormalizedConductivity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double J = 1 + fp.m_ef;
    double lnJ = log(J);
    double C[MAX_COEF];
    for (int k=0; k<m_ncoef; ++k) C[k] = m_coef[k]->value(T);
    double khat = C[0];
    for (int k=1; k<m_ncoef; ++k) khat += C[k]*pow(lnJ,k);
    
    return khat;
}

//-----------------------------------------------------------------------------
//! tangent of normalized viscosity with respect to temperature
double FEThermalCondLnJvirial::Tangent_NormalizedConductivity_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double J = 1 + fp.m_ef;
    double lnJ = log(J);
    double dC[MAX_COEF];
    for (int k=0; k<m_ncoef; ++k) dC[k] = m_coef[k]->derive(T);
    double dkhatdT = dC[0];
    for (int k=1; k<m_ncoef; ++k) dkhatdT += dC[k]*pow(lnJ,k);
    
    return dkhatdT;
}

//-----------------------------------------------------------------------------
//! tangent of normalized viscosity with respect to volumetric strain (or J)
double FEThermalCondLnJvirial::Tangent_NormalizedConductivity_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T + m_Tr;
    double J = 1 + fp.m_ef;
    double lnJ = log(J);
    double C[MAX_COEF];
    for (int k=1; k<m_ncoef; ++k) C[k] = m_coef[k]->value(T);
    double dkhatdJ = 0;
    for (int k=1; k<m_ncoef; ++k) dkhatdJ += k*C[k]*pow(lnJ,k-1);
    dkhatdJ /= J;
    
    return dkhatdJ;
}
