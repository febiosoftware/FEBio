/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
bool FEPorousNeoHookean::Init()
{
    // If porosity is not specified (i.e., if equal to 0) and parent material
    // is multiphasic, extract solid volume fraction from it
    if (m_phisr == 1) {
        FECoreBase* m_pMat = GetParent();
        FEBiphasic* m_pb = dynamic_cast<FEBiphasic*>(m_pMat);
        FEBiphasicSolute* m_pbs = dynamic_cast<FEBiphasicSolute*>(m_pMat);
        FEMultiphasic* m_pm = dynamic_cast<FEMultiphasic*>(m_pMat);
        FETriphasic* m_pt = dynamic_cast<FETriphasic*>(m_pMat);
        if (m_pb) m_phisr = m_pb->m_phi0;
        else if (m_pbs) m_phisr = m_pbs->m_phi0;
        else if (m_pm) m_phisr = m_pm->m_phi0;
        else if (m_pt) m_phisr = m_pt->m_phi0;
    }
    
    m_phiwr = 1. - m_phisr;
    
    // lame parameters (ignore user-specified Poisson ratio)
    m_lam = 0;
    m_mu  = m_E/3*(1+0.5*m_phiwr*m_phiwr);
    
    return true;
}

//-----------------------------------------------------------------------------
mat3ds FEPorousNeoHookean::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double J = pt.m_J;
    double Jbar = (J-m_phisr)/m_phiwr;
    double lnJbar = log(Jbar);
    double R = pow(Jbar/J, 2./3.);
    double mu1 = m_mu/J;

    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    double I1 = b.tr();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = b*(mu1*R) + I*((mu1*(m_phisr*R*I1/3. - J) + m_lam*lnJbar)/(J - m_phisr));
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPorousNeoHookean::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double J = pt.m_J;
    double Jbar = (J-m_phisr)/m_phiwr;
    double lnJbar = log(Jbar);
    double R = pow(Jbar/J, 2./3.);
    double mu1 = m_mu/J;
    double lam1 = m_lam/J;
    
    double g = mu1*R*m_phisr/(J-m_phisr);
    double h = (lam1*lnJbar - mu1)*J/(J-m_phisr);
    double Jdg = mu1*R*(2*m_phisr - 3*J)*m_phisr/3./pow(J - m_phisr, 2);
    double Jdh = J*(mu1*m_phisr + lam1*(J - m_phisr*lnJbar))/pow(J - m_phisr, 2);

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
    
    double J = pt.m_J;
    double Jbar = (J-m_phisr)/m_phiwr;
    double lnJbar = log(Jbar);
    
    // calculate left Cauchy-Green tensor
    mat3ds b = pt.LeftCauchyGreen();
    double I1bar = b.tr()*pow(Jbar/J, 2./3.);
    
    double sed = m_mu*((I1bar-3)/2.0 - lnJbar)+m_lam*lnJbar*lnJbar/2.0;
    
    return sed;
}
