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
#include "FEReactionRateHuiskes.h"
#include "FEBiphasic.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateHuiskes, FEReactionRate)
	ADD_PARAMETER(m_B, "B");
	ADD_PARAMETER(m_psi0, FE_RANGE_GREATER_OR_EQUAL(0.0), "psi0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateHuiskes::FEReactionRateHuiskes(FEModel* pfem) : FEReactionRate(pfem) 
{ 
    m_B = 0;
    m_psi0 = 0;
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateHuiskes::ReactionRate(FEMaterialPoint& pt)
{
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    double rhor = pbm->SolidReferentialApparentDensity(pt);
    double phir = pbm->SolidReferentialVolumeFraction(pt);
	
    FERemodelingMaterialPoint& rpt = *(pt.ExtractData<FERemodelingMaterialPoint>());
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.m_J;
	double sed = rpt.m_sed;
	double zhat = m_B(pt)*(sed/rhor - m_psi0(pt))/(J-phir);
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateHuiskes::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());

    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    
    double rhor = pbm->SolidReferentialApparentDensity(pt);
    double phir = pbm->SolidReferentialVolumeFraction(pt);
    double p = pbm->GetActualFluidPressure(pt);
	
    double J = et.m_J;
    double zhat = ReactionRate(pt);
    mat3dd I(1);
    mat3ds dzhatde = (I*(-zhat) + (et.m_s+I*p)*(m_B(pt)/rhor))/(J-phir);
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateHuiskes::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

