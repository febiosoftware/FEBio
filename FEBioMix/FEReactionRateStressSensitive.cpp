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
#include "FEReactionRateStressSensitive.h"
#include "FEBiphasic.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include "FEBioMech/FENeoHookean.h"
#include <FEBioMech/FEKinematicGrowth.h>
#include <FEBioMech/FEGrowthTensor.h>
#include <FECore\FESolidElement.h>
#include <FECore\FESPRProjection.h>
#include <FECore\FESolidDomain.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateStressSensitive, FEReactionRate)
	ADD_PARAMETER(m_k, "k");
	ADD_PARAMETER(m_a0, "a0");
	ADD_PARAMETER(m_a, "a");
	ADD_PARAMETER(m_b, "b");
	ADD_PARAMETER(stress0, "residual_stress")->setLongName("initial residual stress");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateStressSensitive::FEReactionRateStressSensitive(FEModel* pfem) : FEReactionRate(pfem)
{ 

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateStressSensitive::ReactionRate(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();
	double m_S = m_a0 + m_a / (1.0 + (exp(-(ep.m_s.tr() - m_b) / stress0)));
	double zhat = m_k(pt) * m_S;
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
//! SL! Todo: Figure out what to do with rhor for solutes.
mat3ds FEReactionRateStressSensitive::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	double phir = pbm->SolidReferentialVolumeFraction(pt);
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	double J = et.m_J;
	vec3d pos = pt.m_r0;

	mat3ds I = mat3ds(1.0);
	double zhat = ReactionRate(pt);
	mat3ds dzde = -1.0 * zhat * I;
	return dzde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateStressSensitive::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

