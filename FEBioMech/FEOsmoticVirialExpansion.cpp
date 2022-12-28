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
#include "FEOsmoticVirialExpansion.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOsmoticVirialExpansion, FEElasticMaterial)
    ADD_PARAMETER(m_phiwr, FE_RANGE_CLOSED(0.0, 1.0), "phiw0");
    ADD_PARAMETER(m_cr   , FE_RANGE_GREATER_OR_EQUAL(0.0), "cr");
    ADD_PARAMETER(m_c1, "c1");
    ADD_PARAMETER(m_c2, "c2");
    ADD_PARAMETER(m_c3, "c3");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEOsmoticVirialExpansion::FEOsmoticVirialExpansion(FEModel* pfem) : FEElasticMaterial(pfem)
{ 
	m_c1 = m_c2 = m_c3 = 0; 
}

//-----------------------------------------------------------------------------
mat3ds FEOsmoticVirialExpansion::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // calculate concentration in current configuration
    double c = m_phiwr*m_cr/(J-1+m_phiwr);
    
    // calculate osmotic pressure
    double p = m_c1*c + m_c2*c*c + m_c3*pow(c,3);
    
    // calculate T = -p*I
    mat3dd I(1.0);	// identity tensor
    mat3ds s = -p*I;
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEOsmoticVirialExpansion::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // calculate concentration in current configuration
    double c = m_phiwr*m_cr/(J-1+m_phiwr);
    
    // calculate osmotic pressure
    double p = m_c1*c + m_c2*c*c + m_c3*pow(c,3);
    
    // calculate derivative of osmotic pressure w.r.t. J
    double dpdJ = -(m_c1 + 2*m_c2*c + 3*m_c3*c*c)*c/(J -1 + m_phiwr);
    
    mat3dd I(1.0);	// Identity
    
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    
    // calculate tangent osmotic modulus
    tens4ds C = IxI*(-J*dpdJ) + (2.0*I4 - IxI)*p;
    
    return C;
}
