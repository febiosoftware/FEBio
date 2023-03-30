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
#include "FEReactionRatePositional.h"
#include "FEBiphasic.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include "FEReaction.h"
#include "FESoluteInterface.h"
#include <FECore\log.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRatePositional, FEReactionRate)
	ADD_PARAMETER(m_k, "k");
	ADD_PARAMETER(m_b, "b");
	ADD_PARAMETER(m_n, "n");
	//ADD_PARAMETER(m_psi0, FE_RANGE_GREATER_OR_EQUAL(0.0), "psi0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRatePositional::FEReactionRatePositional(FEModel* pfem) : FEReactionRate(pfem) 
{	

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRatePositional::ReactionRate(FEMaterialPoint& pt)
{
	// find whether we have solutes or reactants. Will need to update later.
	int nsol = m_pReact->m_psm->Solutes();
	int nsbm = m_pReact->m_psm->SBMs();
	if (nsol > 0 || nsbm > 0)
	{
		FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
		double phir = pbm->SolidReferentialVolumeFraction(pt);
		FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
		double J = et.m_J;
		vec3d pos = pt.m_r0;
		double pos_ab = pow(fabs(pos.x) + fabs(pos.y), m_n);
		double r = exp(-pos_ab / m_b);
		double zhat = m_k * r / (J - phir);
		return zhat;
	}
	else
	{
		return 0.0;
	}

 
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRatePositional::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	double phir = pbm->SolidReferentialVolumeFraction(pt);
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	double J = et.m_J;
	vec3d pos = pt.m_r0;

	mat3ds I = mat3ds(1.0);
	double zhat = ReactionRate(pt);
	mat3ds dzde = -1.0 * (zhat / (J - phir)) * I;
	return dzde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRatePositional::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

