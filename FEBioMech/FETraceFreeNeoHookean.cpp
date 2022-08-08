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



#include "stdafx.h"
#include "FETraceFreeNeoHookean.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FETraceFreeNeoHookean, FEElasticMaterial)
    ADD_PARAMETER(m_mu, FE_RANGE_GREATER(0.0), "mu")->setLongName("shear modulus");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FETraceFreeNeoHookean::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double J = pt.m_J;
    
    // get the material parameter
    double mu = m_mu(mp);
    
    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    
    double I1 = b.tr();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = (b*(3/I1) - I)*(mu/J);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FETraceFreeNeoHookean::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double J = pt.m_J;
    
    // get the material parameter
    double mu = m_mu(mp);

    mat3dd I(1);
    
    return dyad4s(I)*(2*mu/J);
}

//-----------------------------------------------------------------------------
double FETraceFreeNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double J = pt.m_J;
    
    // get the material parameter
    double mu = m_mu(mp);

    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    double I1 = b.tr();
    
    double sed = mu*log(pow(I1/3,1.5)/J);
    
    return sed;
}

//-----------------------------------------------------------------------------
mat3ds FETraceFreeNeoHookean::PK2Stress(FEMaterialPoint& mp, const mat3ds ES)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the material parameter
    double mu = m_mu(mp);
    
    mat3ds C = pt.RightCauchyGreen();
    double I1 = C.tr();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds S = (I*(3/I1) - C.inverse())*mu;
    
    return S;
}

//-----------------------------------------------------------------------------
tens4dmm FETraceFreeNeoHookean::MaterialTangent(FEMaterialPoint& mp, const mat3ds ES)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the material parameter
    double mu = m_mu(mp);
    
    mat3ds C = pt.RightCauchyGreen();
    tens4dmm c = dyad4s(C.inverse())*(2*mu);
    
    return c;
}
