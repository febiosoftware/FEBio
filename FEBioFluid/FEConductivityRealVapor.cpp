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

#include "FEConductivityRealVapor.h"
#include "FEThermoFluid.h"
#include "FERealVapor.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEConductivityRealVapor, FEFluidThermalConductivity)

    // parameters
    ADD_PARAMETER(m_Kr, "Kr")->setLongName("referential thermal conductivity")->setUnits(UNIT_THERMAL_CONDUCTIVITY);
    // properties
    ADD_PROPERTY(m_esat , "esat", FEProperty::Optional)->SetLongName("saturation dilatation");
    ADD_PROPERTY(m_Ksat , "Ksat", FEProperty::Optional)->SetLongName("normalized saturation thermal conductivity");
    ADD_PROPERTY(m_C[0] , "C0"  , FEProperty::Optional)->SetLongName("1st K virial coeff");
    ADD_PROPERTY(m_C[1] , "C1"  , FEProperty::Optional)->SetLongName("2nd K virial coeff");
    ADD_PROPERTY(m_C[2] , "C2"  , FEProperty::Optional)->SetLongName("3rd K virial coeff");
    ADD_PROPERTY(m_C[3] , "C3"  , FEProperty::Optional)->SetLongName("4th K virial coeff");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConductivityRealVapor::FEConductivityRealVapor(FEModel* pfem) : FEFluidThermalConductivity(pfem)
{
    m_Ksat = nullptr;
    m_esat = nullptr;
    for (int k=0; k<MAX_NVC; ++k) m_C[k] = nullptr;
    m_Kr = 0;
    m_Tr = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEConductivityRealVapor::Init()
{
    m_Tr = GetGlobalConstant("T");
    
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    
    if (m_esat) m_esat->Init();
    else {
        FECoreBase* pMat = GetAncestor();
        FEThermoFluid* pFluid = dynamic_cast<FEThermoFluid*>(pMat);
        if (pFluid) {
            FERealVapor* pRV = dynamic_cast<FERealVapor*>(pFluid->GetElastic());
            if (pRV) {
                m_esat = pRV->m_esat;
                m_esat->Init();
                m_Tc = pRV->m_Tc;
                m_alpha = pRV->m_alpha;
            }
            else return false;
        }
        else return false;
    }
    if (m_Ksat) m_Ksat->Init();
    m_nvc = 0;
    for (int k=0; k<MAX_NVC; ++k) {
        if (m_C[k]) {
            m_C[k]->Init();
            ++m_nvc;
        }
    }

    return FEFluidThermalConductivity::Init();
}

//-----------------------------------------------------------------------------
void FEConductivityRealVapor::Serialize(DumpStream& ar)
{
    FEFluidThermalConductivity::Serialize(ar);
    
    if (ar.IsShallow()) return;
    
    ar & m_Kr & m_Tr & m_nvc;
    ar & m_Ksat;
    for (int i=0; i<MAX_NVC; ++i)
        ar & m_C[i];
    
    if (ar.IsLoading()) {
        if (m_Ksat) m_Ksat->Init();
        for (int i=0; i<MAX_NVC; ++i)
            if (m_C[i]) m_C[i]->Init();
    }
}

//-----------------------------------------------------------------------------
//! calculate thermal conductivity at material point
double FEConductivityRealVapor::ThermalConductivity(FEMaterialPoint& mp)
{
    double K = 1.0;
    if (m_Ksat) {
        FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
        FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();
        double T = tf.m_T + m_Tr;
        double That = T/m_Tr;
        double J = 1 + pf.m_ef;
        double y = (That < m_Tc) ? (m_Tc-That)/(m_Tc-1) : 0;
        double q = log(1+pow(y,m_alpha));
        K = m_Ksat->value(q);
        double Jsat = exp(m_esat->value(q));
        double x = 1 - Jsat/J;
        for (int k=0; k<m_nvc; ++k) K += m_C[k]->value(q)*pow(x,k+1);
    }
    return K*m_Kr;
}

//-----------------------------------------------------------------------------
//! tangent of thermal conductivity with respect to strain J
double FEConductivityRealVapor::Tangent_Strain(FEMaterialPoint& mp)
{
    double d = 1e-6;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef+d;
    ft->m_T = tf.m_T;
    FEMaterialPoint tmp(ft);
    double Kp = ThermalConductivity(tmp);
    fp->m_ef = pf.m_ef-d;
    double Km = ThermalConductivity(tmp);
    delete ft;
    double dKJ = (Kp - Km)/(2*d);
    return dKJ;
}

//-----------------------------------------------------------------------------
//! tangent of thermal conductivity with respect to temperature T
double FEConductivityRealVapor::Tangent_Temperature(FEMaterialPoint& mp)
{
    double Tr = GetGlobalConstant("T");
    double d = 1e-6*Tr;
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    FEFluidMaterialPoint& pf = *mp.ExtractData<FEFluidMaterialPoint>();

    FEFluidMaterialPoint* fp = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* ft = new FEThermoFluidMaterialPoint(fp);
    fp->m_ef = pf.m_ef;
    ft->m_T = tf.m_T+d;
    FEMaterialPoint tmp(ft);
    double Kp = ThermalConductivity(tmp);
    ft->m_T = tf.m_T-d;
    double Km = ThermalConductivity(tmp);
    delete ft;
    double dKT = (Kp - Km)/(2*d);
    return dKT;
}

