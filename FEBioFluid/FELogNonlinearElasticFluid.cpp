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


#include "FELogNonlinearElasticFluid.h"
#include "FEFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
FELogNonlinearElasticFluid::FELogNonlinearElasticFluid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_k = 0;
    m_rhor = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FELogNonlinearElasticFluid::Init()
{
    return true;
}

//-----------------------------------------------------------------------------
// serialization
void FELogNonlinearElasticFluid::Serialize(DumpStream& ar)
{
    if (ar.IsLoading()) return;
    ar & m_k & m_rhor;
}

//-----------------------------------------------------------------------------
//! gage pressure
double FELogNonlinearElasticFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1+fp.m_ef;
    return -m_k*log(J)/J;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FELogNonlinearElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1+fp.m_ef;
    return m_k*(log(J)-1)/pow(J,2);
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FELogNonlinearElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1+fp.m_ef;
    return m_k*(3-2*log(J))/pow(J,3);
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to temperature T
double FELogNonlinearElasticFluid::Tangent_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to temperature T
double FELogNonlinearElasticFluid::Tangent_Temperature_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J and temperature T
double FELogNonlinearElasticFluid::Tangent_Strain_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific free energy
double FELogNonlinearElasticFluid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = 1+fp.m_ef;
    return m_k/2*pow(log(J),2)/m_rhor;
}

//-----------------------------------------------------------------------------
//! specific entropy
double FELogNonlinearElasticFluid::SpecificEntropy(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! specific strain energy
double FELogNonlinearElasticFluid::SpecificStrainEnergy(FEMaterialPoint& mp)
{
    return SpecificFreeEnergy(mp);
}

//-----------------------------------------------------------------------------
//! isobaric specific heat capacity
double FELogNonlinearElasticFluid::IsobaricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! isochoric specific heat capacity
double FELogNonlinearElasticFluid::IsochoricSpecificHeatCapacity(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to strain J
double FELogNonlinearElasticFluid::Tangent_cv_Strain(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of isochoric specific heat capacity with respect to temperature T
double FELogNonlinearElasticFluid::Tangent_cv_Temperature(FEMaterialPoint& mp)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FELogNonlinearElasticFluid::Dilatation(const double T, const double p, double& e)
{
    double errabs = 1e-3;
    double errrel = 1e-3;
    int maxiter = 50;
    bool cnvgd = false;
    // initializations
    e = -p/(m_k+p);                     // initial guess
    double f = log(1+e)+p*(1+e)/m_k;    // function
    double de = 0;
    int iter = 0;
    while (!cnvgd) {
        double df = p/m_k+1.0/(1+e);    // derivative
        double de = -f/df;              // solution increment
        e += de;                        // update solution
        f = log(1+e)+p*(1+e)/m_k;       // function
        // check convergence
        if ((fabs(f) < errabs) || (fabs(de) < errrel*fabs(e))) cnvgd = true;
        else if (++iter > maxiter) return false;
    }
    return true;
}
