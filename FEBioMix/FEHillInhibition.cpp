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
#include "FESoluteInterface.h"
#include "FEBiphasic.h"
#include "FEHillInhibition.h"

BEGIN_FECORE_CLASS(FEHillInhibition, FEReactionRate)
	ADD_PARAMETER(m_Kmax, "Kmax");
	ADD_PARAMETER(m_w, "reaction_weight");
	ADD_PARAMETER(m_t, "degradation_rate");
	ADD_PARAMETER(m_E50, "E_50");
	ADD_PARAMETER(m_n, "Hill_coeff");
	ADD_PARAMETER(m_sol_id, "sol_id");
	ADD_PARAMETER(m_sbm_id, "sbm_id");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEHillInhibition::FEHillInhibition(FEModel* pfem) : FEReactionRate(pfem)
{

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEHillInhibition::ReactionRate(FEMaterialPoint& pt)
{
	double c = 0.0;
	if (m_sol_id > 0)
		c = m_pReact->m_psm->GetActualSoluteConcentration(pt, m_sol_id - 1);
	else if (m_sbm_id > 0)
		c = m_pReact->m_psm->SBMConcentration(pt, m_sbm_id - 1);
	double En = pow(m_E50, m_n);
	double B = (En - 1.0) / (2.0 * En - 1.0);
	double K = pow((B - 1.0), (1 / m_n));
	double cn = pow(c, m_n);
	double Kn = pow(K, m_n);
	double zhat = (m_Kmax * m_w / m_t) * (1.0 - ((B * cn) / (Kn + cn)));
	return max(zhat, 0.0);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEHillInhibition::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	// if the reaction supply is insensitive to strain
	if (m_pReact->m_bool_refC)
		return mat3ds(0.0);

	double zhat = ReactionRate(pt);
	mat3dd I(1);
	mat3ds dzhatde = I * (-zhat);
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEHillInhibition::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0.0;
}
