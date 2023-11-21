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
#include "FEHillActivationORActivation.h"

BEGIN_FECORE_CLASS(FEHillActivationORActivation, FEReactionRate)
	ADD_PARAMETER(m_Kmax, "Kmax");
	ADD_PARAMETER(m_w, "reaction_weight");
	ADD_PARAMETER(m_t, "degradation_rate");
	ADD_PARAMETER(m_E50, "E_50");
	ADD_PARAMETER(m_n, "Hill_coeff");
	ADD_PARAMETER(m_sol_id[0], "sol_id_act1");
	ADD_PARAMETER(m_sbm_id[0], "sbm_id_act1");
	ADD_PARAMETER(m_sol_id[1], "sol_id_act2");
	ADD_PARAMETER(m_sbm_id[1], "sbm_id_act2");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEHillActivationORActivation::FEHillActivationORActivation(FEModel* pfem) : FEReactionRate(pfem)
{

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEHillActivationORActivation::ReactionRate(FEMaterialPoint& pt)
{
	double c[2] = { 0.0 };
	double En = pow(m_E50, m_n);
	double B = (En - 1.0) / (2.0 * En - 1.0);
	double K = pow((B - 1.0), (1.0 / m_n));
	double Kn = pow(K, m_n);
	double cn[2] = { 0.0 };
	double F[2] = { 0.0 };
	for (int i = 0; i < 2; ++i)
	{
		if (m_sol_id[i] > 0)
			c[i] = m_pReact->m_psm->GetActualSoluteConcentration(pt, m_sol_id[i] - 1);
		else if (m_sbm_id[i] > 0)
			c[i] = m_pReact->m_psm->SBMConcentration(pt, m_sbm_id[i] - 1);
		cn[i] = pow(c[i], m_n);
		F[i] = (B * cn[i]) / (Kn + cn[i]);
	}
	double zhat = (m_Kmax * m_w / m_t) * (F[0] + F[1] - F[0] * F[1]);
	return max(zhat, 0.0);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEHillActivationORActivation::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
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
double FEHillActivationORActivation::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0.0;
}
