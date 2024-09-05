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
#include "FEReactionRateStressSensitiveERD.h"
#include <FEBioMix/FEBiphasic.h>
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include "FEBioMech/FENeoHookean.h"
#include "FEKinematicGrowthRateDependent.h"
#include "FEGrowthTensorERD.h"
#include <FECore/FESolidElement.h>
#include <FECore/FESPRProjection.h>
#include <FECore/FESolidDomain.h>
#include <iostream>
#include "FEElasticReactionDiffusionSolidDomain.h"
#include "FEElasticReactionDiffusionDomain.h"
#include "FEMassActionForwardERD.h"

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateStressSensitiveERD, FEReactionRateERD)
	ADD_PARAMETER(m_k, "k");
	ADD_PARAMETER(m_a0, "a0");
	ADD_PARAMETER(m_a, "a");
	ADD_PARAMETER(m_b, "b");
	ADD_PARAMETER(m_s0, "s0")->setLongName("initial residual stress");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateStressSensitiveERD::ReactionRate(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();
	double I1 = ep.m_s.tr();
	double m_S = m_a0 + m_a / (1.0 + (exp(-(I1 - m_s0) / m_b)));
	
	return m_k(pt) * m_S;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
//! SL! Todo: Figure out what to do with rhor for solutes.
mat3ds FEReactionRateStressSensitiveERD::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	return mat3ds(0.0);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate k with Cauchy stress (sigma) at material point
mat3ds FEReactionRateStressSensitiveERD::Tangent_ReactionRate_Stress(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();
	FEElement* el = pt.m_elem;

	double trs = ep.m_s.tr();
	double m_S = m_a0 + m_a / (1.0 + (exp(-(trs - m_s0) / m_b)));	
	double dkfds = m_k(pt);
	double temp_exp = exp(-(trs - m_s0) / m_b);
	double dsdtrs = (m_a / m_b) * ((temp_exp) / pow((1.0 + temp_exp),2.0));
	mat3ds dtrsdsigma = mat3ds(1.0);
	return dkfds * dsdtrs * dtrsdsigma;

}