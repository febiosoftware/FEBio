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


#include "FEIdealGas.h"
#include <FECore/log.h>
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluid.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEIdealGas, FEElasticFluid)

    // material parameters
    ADD_PARAMETER(m_M    , FE_RANGE_GREATER(0.0), "M")->setUnits(UNIT_MOLAR_MASS)->setLongName("molar mass");
    ADD_PARAMETER(m_arhat, "ar")->setLongName("normalized referential specific free energy");    // ar normalized by R.Tr/M
    ADD_PARAMETER(m_srhat, "sr")->setLongName("normalized referential specific entropy");        // sr normalized by R/M
    ADD_PROPERTY (m_aohat, "ao")->SetLongName("normalized specific free energy circle");         // a-circle normalized by R.Tr/M
    ADD_PROPERTY (m_cphat, "cp")->SetLongName("normalized isobaric specific heat capacity");     // cp normalized by R/M

END_FECORE_CLASS();

FEIdealGas::FEIdealGas(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_R = m_Pr = m_Tr = m_arhat = m_srhat = 0;
    m_aohat = nullptr;
    m_cphat = nullptr;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGas::Init()
{
    m_R  = GetGlobalConstant("R");
    m_Tr = GetGlobalConstant("T");
    m_Pr = GetGlobalConstant("P");
    
    if (m_R  <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section");    return false; }
    if (m_Tr <= 0) { feLogError("A positive referential absolute temperature T must be defined in Globals section"); return false; }
    if (m_Pr == 0) {
        FEThermoFluid* pMat = dynamic_cast<FEThermoFluid*>(GetParent());
        double rhor = pMat->ReferentialDensity();
        m_Pr = m_R*m_Tr/m_M*rhor;
        feLogWarning("The referential absolute pressure P is calculated internally as %g\n",m_Pr);
    }
    
    m_aohat->Init();
    m_cphat->Init();
    
    return true;
}

//-----------------------------------------------------------------------------
void FEIdealGas::Serialize(DumpStream& ar)
{
    FEElasticFluid::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_R & m_Pr & m_Tr;
}

//-----------------------------------------------------------------------------
//! gauge pressure
double FEIdealGas::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (m_Tr + tf.m_T)/m_Tr;
    double J = 1 + fp.m_ef;

    double phat = That/J-1;

    return phat*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FEIdealGas::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (m_Tr + tf.m_T)/m_Tr;
    double J = 1 + fp.m_ef;

    double dphatJ = -That/pow(J,2);

    return dphatJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FEIdealGas::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (m_Tr + tf.m_T)/m_Tr;
    double J = 1 + fp.m_ef;

    double d2phatJ = 2*That/pow(J,3);

    return d2phatJ*m_Pr;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FEIdealGas::Tangent_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;

    double dphatThat = 1/J;

    return dphatThat*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FEIdealGas::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FEIdealGas::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1 + fp.m_ef;

    double d2phat = -1/pow(J,2);

    return d2phat*m_Pr/m_Tr;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FEIdealGas::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double That = (m_Tr + tf.m_T)/m_Tr;
    
    // referential free energy
    double ahat = J + That*(log(That/J)-1)+ m_arhat - m_srhat*(That-1);
    
    // add a_circle
    ahat += m_aohat->value(That);
    
    return ahat*m_R*m_Tr/m_M;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FEIdealGas::SpecificEntropy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double That = (tf.m_T + m_Tr)/m_Tr;

    double shat = -log(That/J)+m_srhat;
    
    // add strain-dependent contribution
    shat += -m_aohat->derive(That);

    return shat*m_R/m_M;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FEIdealGas::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    
    double J = 1 + fp.m_ef;
    double That = (m_Tr + tf.m_T)/m_Tr;
    
    double what = J + That*(log(That/J)-1);
    
    return what*m_R*m_Tr/m_M;
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FEIdealGas::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (tf.m_T + m_Tr)/m_Tr;

    double cphat = m_cphat->value(That);
    
    return cphat*m_R/m_M;
}
                
//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FEIdealGas::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    double cv = IsobaricSpecificHeatCapacity(mp) - m_R/m_M;

    return cv;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FEIdealGas::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FEIdealGas::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    FEThermoFluidMaterialPoint& tf = *mp.ExtractData<FEThermoFluidMaterialPoint>();
    double That = (tf.m_T + m_Tr)/m_Tr;

    double dcvhat = m_cphat->derive(That);
    
    return dcvhat*m_R*m_Tr/m_M;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FEIdealGas::Dilatation(const double T, const double p, double& e)
{
    double J = (T+m_Tr)/m_Tr/(1+p/m_Pr);
    e = J - 1;
    return true;
}
