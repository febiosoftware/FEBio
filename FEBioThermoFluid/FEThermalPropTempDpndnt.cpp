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



#include "FEThermalPropTempDpndnt.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FEBioFluid/FEFluidMaterialPoint.h>
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEThermalPropTempDpndnt, FEThermoFluidProperty)
    ADD_PROPERTY(m_prop , "prophat"   )->SetLongName("normalized property");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEThermalPropTempDpndnt::FEThermalPropTempDpndnt(FEModel* pfem) : FEThermoFluidProperty(pfem)
{
    m_prop = nullptr;
}

bool FEThermalPropTempDpndnt::Init()
{
    m_Tr = GetGlobalConstant("T");
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }

    if (m_prop) return m_prop->Init();
    
    return false;
}

//-----------------------------------------------------------------------------
//! normalize viscosity
double FEThermalPropTempDpndnt::NormalizedProperty(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T/m_Tr+1;
    double prop = m_prop->value(T);
    
    return prop;
}

//-----------------------------------------------------------------------------
//! tangent of normalized property with respect to temperature
double FEThermalPropTempDpndnt::Tangent_NormalizedProperty_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double T = tf.m_T/m_Tr+1;
    double dpropdT = m_prop->derive(T);
    
    return dpropdT;
}

//-----------------------------------------------------------------------------
//! tangent of normalized property with respect to volumetric strain (or J)
double FEThermalPropTempDpndnt::Tangent_NormalizedProperty_Strain(FEMaterialPoint& mp)
{
    return 0.0;
}
