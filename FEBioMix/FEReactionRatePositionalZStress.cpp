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
#include "FEReactionRatePositionalZStress.h"
#include "FEBiphasic.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include <FECore\log.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRatePositionalZStress, FEReactionRate)
ADD_PARAMETER(stress0, "residual_stress")->setLongName("initial residual stress");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRatePositionalZStress::FEReactionRatePositionalZStress(FEModel* pfem) : FEReactionRate(pfem)
{

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRatePositionalZStress::ReactionRate(FEMaterialPoint& pt)
{
    double z_pos = pt.m_r0.z;
    double xpo = -1.0 * pow(fabs(z_pos), 3) / 3.0;    
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    double phir = pbm->SolidReferentialVolumeFraction(pt);

    double J = et.m_J;
    double m_S = 1.0 + (fabs(et.m_s.tr()) / stress0);
    double zhat = 0.5 * m_S * exp(xpo);
    return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
//! SL! Todo: Figure out what to do with rhor for solutes.
mat3ds FEReactionRatePositionalZStress::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    return mat3ds(0.0);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRatePositionalZStress::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
    return 0;
}

