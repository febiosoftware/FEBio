/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
//
//  FEFluidConstantConductivity.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/28/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#include "FETempDependentConductivity.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FETempDependentConductivity, FEFluidThermalConductivity)

    // parameters
    ADD_PARAMETER(m_Kr, "Kr")->setLongName("referential thermal conductivity")->setUnits(UNIT_THERMAL_CONDUCTIVITY);
    // properties
    ADD_PROPERTY(m_Khat, "Khat")->SetLongName("normalized thermal conductivity");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETempDependentConductivity::FETempDependentConductivity(FEModel* pfem) : FEFluidThermalConductivity(pfem)
{
    m_Khat = nullptr;
    m_Kr = 0;
    m_Tr = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FETempDependentConductivity::Init()
{
    m_Tr = GetGlobalConstant("T");
    
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    
    m_Khat->Init();
    
    return FEFluidThermalConductivity::Init();
}

//-----------------------------------------------------------------------------
void FETempDependentConductivity::Serialize(DumpStream& ar)
{
    FEFluidThermalConductivity::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_Kr & m_Tr;
    ar & m_Khat;
    
    if (ar.IsLoading()) {
        m_Khat->Init();
    }
}

//-----------------------------------------------------------------------------
//! calculate thermal conductivity at material point
double FETempDependentConductivity::ThermalConductivity(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (tf.m_T+m_Tr)/m_Tr;
    return m_Khat->value(That)*m_Kr;
}

//-----------------------------------------------------------------------------
//! tangent of thermal conductivity with respect to strain J
double FETempDependentConductivity::Tangent_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of thermal conductivity with respect to temperature T
double FETempDependentConductivity::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (tf.m_T+m_Tr)/m_Tr;
    return m_Khat->derive(That)*m_Kr/m_Tr;
}

