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


#include "FELinearThermoElasticFluid.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FEBioFluid/FEFluidMaterialPoint.h>
#include "FEThermalPropConst.h"
#include "FEThermoFluid.h"

// define the material parameters
BEGIN_FECORE_CLASS(FELinearThermoElasticFluid, FEThermoElasticFluid)
    ADD_PARAMETER(m_cvr, FE_RANGE_GREATER_OR_EQUAL(0.0), "cvr")->setUnits(UNIT_SPECIFIC_ENTROPY)->setLongName("referential isochoric specific heat capacity");
    ADD_PROPERTY(m_cvhat, "cvhat", FEProperty::Optional)->SetLongName("normalized isochoric specific heat capacity");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
FELinearThermoElasticFluid::FELinearThermoElasticFluid(FEModel* pfem) : FEThermoElasticFluid(pfem)
{
    m_cvhat = nullptr;
    m_Tr = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FELinearThermoElasticFluid::Init()
{
    FEModel* pfem = GetFEModel();
    
    m_Tr = GetGlobalConstant("T");
    if (m_Tr <= 0){
        feLogError("A positive absolute temperature T must be defined in Globals section");
        return false;
    }
    
    if (m_cvhat == nullptr) {
        m_cvhat = new FEThermalPropConst(pfem);
    }
    m_cvhat->Init();
    
    return true;
}

//-----------------------------------------------------------------------------
//! gauge pressure
double FELinearThermoElasticFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
	double k = m_pFluid->m_k;
    double p = -k*fp.m_ef;
    
    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FELinearThermoElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    return -m_pFluid->m_k;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FELinearThermoElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FELinearThermoElasticFluid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
	double k = m_pFluid->m_k;
	double rhor = m_pFluid->m_rhor;

    double a = k/(2*rhor)*pow(fp.m_ef,2);
    return a;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FELinearThermoElasticFluid::Dilatation(const double T, const double p, double& e)
{
    e = -p/m_pFluid->m_k;
    return true;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FELinearThermoElasticFluid::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    double cv = m_cvr;
    if (m_cvhat) cv *= m_cvhat->NormalizedProperty(mp);
    return cv;
}
        
//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FELinearThermoElasticFluid::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double cv = IsochoricSpecificHeatCapacity(mp);
    double p = Pressure(mp);    // current gauge pressure
    double dpT = Tangent_Temperature(mp);
    // evaluate dpT/dpJ in the reference configuration
    double efsafe = fp.m_ef;
    double Tsafe = tf.m_T;
    fp.m_ef = 0;
    tf.m_T = 0;
    double r = Tangent_Temperature(mp)/Tangent_Strain(mp);
    // restore values before continuing
    fp.m_ef = efsafe;
    tf.m_T = Tsafe;
    // evaluate isobaric specific heat capacity
    double cp = cv + r/m_pFluid->m_rhor*(p - (m_Tr + tf.m_T)*dpT);
    return cp;
}
        
//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FELinearThermoElasticFluid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    double dcvdJ = 0;
    if (m_cvhat) dcvdJ = m_cvr*m_cvhat->Tangent_NormalizedProperty_Strain(mp);
    return dcvdJ;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FELinearThermoElasticFluid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    double dcvdT = 0;
    if (m_cvhat) dcvdT = m_cvr*m_cvhat->Tangent_NormalizedProperty_Temperature(mp);
    return dcvdT;
}
