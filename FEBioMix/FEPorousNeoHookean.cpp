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



#include "stdafx.h"
#include "FEPorousNeoHookean.h"
#include "FEBiphasicSolute.h"
#include "FEMultiphasic.h"
#include "FETriphasic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPorousNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E    , FE_RANGE_GREATER   (      0.0), "E"       );
	ADD_PARAMETER(m_phisr, FE_RANGE_CLOSED    (0.0 , 1.0), "phi0"    );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FEPorousNeoHookean::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double phisr = ReferentialSolidVolumeFraction(mp);
    double phiwr = 1 - phisr;
    m_mu  = m_E/3*(1+0.5*phiwr*phiwr);
    double J = pt.m_J;
    double Jbar = (J-phisr)/phiwr;
    double lnJbar = log(Jbar);
    double R = pow(Jbar/J, 2./3.);
    double mu1 = m_mu/J;

    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    double I1 = b.tr();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = b*(mu1*R) + I*((mu1*(phisr*R*I1/3. - J) + m_lam*lnJbar)/(J - phisr));
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPorousNeoHookean::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double phisr = ReferentialSolidVolumeFraction(mp);
    double phiwr = 1 - phisr;
    m_mu  = m_E/3*(1+0.5*phiwr*phiwr);
    double J = pt.m_J;
    double Jbar = (J-phisr)/phiwr;
    double lnJbar = log(Jbar);
    double R = pow(Jbar/J, 2./3.);
    double mu1 = m_mu/J;
    double lam1 = m_lam/J;
    
    double g = mu1*R*phisr/(J-phisr);
    double h = (lam1*lnJbar - mu1)*J/(J-phisr);
    double Jdg = mu1*R*(2*phisr - 3*J)*phisr/3./pow(J - phisr, 2);
    double Jdh = J*(mu1*phisr + lam1*(J - phisr*lnJbar))/pow(J - phisr, 2);

    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    double I1 = b.tr();
    
    // Identity
    mat3dd I(1);
    
    tens4ds bIIb = dyad1s(I, b);
    tens4ds II = dyad1s(I);
    tens4ds I4 = dyad4s(I);
    tens4ds c = bIIb*2.*g/3. + II*(Jdg*I1/3 + Jdh) - I4*(2*(g*I1/3 + h));
    
    return c;
}

//-----------------------------------------------------------------------------
double FEPorousNeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double phisr = ReferentialSolidVolumeFraction(mp);
    double phiwr = 1 - phisr;
    m_mu  = m_E/3*(1+0.5*phiwr*phiwr);
    double J = pt.m_J;
    double Jbar = (J-phisr)/phiwr;
    double lnJbar = log(Jbar);
    
    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    double I1bar = b.tr()*pow(Jbar/J, 2./3.);
    
    double sed = m_mu*((I1bar-3)/2.0 - lnJbar)+m_lam*lnJbar*lnJbar/2.0;
    
    return sed;
}

//-----------------------------------------------------------------------------
double FEPorousNeoHookean::ReferentialSolidVolumeFraction(FEMaterialPoint& mp)
{
    if (m_phisr < 1) return m_phisr;
    
    FEBiphasicMaterialPoint& pt = *mp.ExtractData<FEBiphasicMaterialPoint>();
    return pt.m_phi0;
}
