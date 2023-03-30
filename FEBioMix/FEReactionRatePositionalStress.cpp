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
#include "FEReactionRatePositionalStress.h"
#include "FEBiphasic.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include <FECore\log.h>
#include "FEBioMech/FENeoHookean.h"

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRatePositionalStress, FEReactionRate)
	ADD_PARAMETER(m_k,	"k");
	ADD_PARAMETER(m_b,	"b");
	ADD_PARAMETER(m_n,	"n");
	ADD_PARAMETER(C1,	"E");
	ADD_PARAMETER(stress0, "residual_stress")->setLongName("initial residual stress");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRatePositionalStress::FEReactionRatePositionalStress(FEModel* pfem) : FEReactionRate(pfem) 
{ 

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRatePositionalStress::ReactionRate(FEMaterialPoint& pt)
{
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
    double phir = pbm->SolidReferentialVolumeFraction(pt);
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.m_J;    
    double m_S = 1.0 + (fabs(et.m_s.tr()) / stress0);

	vec3d pos = pt.m_r0;
	double pos_ab = pow(fabs(pos.x) + fabs(pos.y), m_n);
	double r = exp(-pos_ab / m_b);

	double zhat = m_k * m_S * r / (J - phir);
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
//! SL! Todo: Figure out what to do with rhor for solutes.
mat3ds FEReactionRatePositionalStress::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	double phir = pbm->SolidReferentialVolumeFraction(pt);
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	double J = et.m_J;
	vec3d pos = pt.m_r0;

	mat3ds I = mat3ds(1.0);
	double zhat = ReactionRate(pt);
	double trs = et.m_s.tr();
	mat3ds dzdJ = -1.0 * (zhat / (J - phir)) * I;
	
	
	double pos_ab = pow(fabs(pos.x) + fabs(pos.y), m_n);
	double r = exp(-pos_ab / m_b);
	FEMesh& mesh = this->GetMesh();
	FEModel* fem = this->GetFEModel();
	mat3ds B = et.LeftCauchyGreen();
	mat3ds dzdS = ((1.0 / (J - phir)) * (4.0 * C1 * m_k * r / J) * (fabs(trs) / stress0)) * B;
	mat3ds dzde = dzdJ + dzde;
	return dzde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRatePositionalStress::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

