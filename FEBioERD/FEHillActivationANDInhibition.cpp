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
#include "FEElasticReactionDiffusionInterface.h"
#include <FEBioMix/FEBiphasic.h>
#include "FEHillActivationANDInhibition.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEHillActivationANDInhibition, FEChemicalReactionERD)
	ADD_PARAMETER(m_Kmax, "Kmax");
	ADD_PARAMETER(m_w, "reaction_weight");
	ADD_PARAMETER(m_t, "degradation_rate");
	ADD_PARAMETER(m_E50, "E_50");
	ADD_PARAMETER(m_n, "Hill_coeff");
	ADD_PARAMETER(u_sol_id_a, "sol_id_act");
	ADD_PARAMETER(u_sol_id_b, "sol_id_inh");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEHillActivationANDInhibition::FEHillActivationANDInhibition(FEModel* pfem) : FEChemicalReactionERD(pfem)
{
	// set material properties
	ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
}

bool FEHillActivationANDInhibition::Init()
{
	m_B = (pow(m_E50, m_n) - 1.0) / (2.0 * pow(m_E50, m_n) - 1.0);
	m_K = pow((m_B - 1.0), (1.0 / m_n));
	m_Kn = pow(m_K, m_n);
	m_Kb = m_Kmax * m_w / m_t;
	m_sol_id_a = u_sol_id_a - 1; m_sol_id_b = u_sol_id_b - 1;
	if ((m_sol_id_a < 0) || (m_sol_id_b < 0))
	{
		feLogError("sol_id: param not valid");
		return false;
	}
	return FEChemicalReactionERD::Init();
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEHillActivationANDInhibition::ReactionSupply(FEMaterialPoint& pt)
{
	double fa = f_Hill(pt, m_sol_id_a);
	double fb = f_Hill(pt, m_sol_id_b);
	double zhat = m_Kb * fa * (1.0 - fb);
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEHillActivationANDInhibition::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	return mat3ds(0.0);
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEHillActivationANDInhibition::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
	double c_act = m_psm->GetActualSoluteConcentration(pt, sol);
	double c_eff = m_psm->GetEffectiveSoluteConcentration(pt, sol);

	//dzdc_act
	if (sol == m_sol_id_a)
	{
		if (c_act > 0.0 && c_eff > 0.0)
			return dfdc(pt, m_sol_id_a) * m_Kb * (1.0 - f_Hill(pt, m_sol_id_b));
	}

	//dzdc_inh
	else if (sol == m_sol_id_b)
	{
		if (c_act > 0.0 && c_eff > 0.0)
			return - 1.0 * dfdc(pt, m_sol_id_b) * m_Kb * f_Hill(pt, m_sol_id_a);
	}

	//dzdc_gamma
	else
		return 0.0;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with Cauchy stress (sigma) at material point
mat3ds FEHillActivationANDInhibition::Tangent_ReactionSupply_Stress(FEMaterialPoint& pt)
{
	return mat3ds(0.0);
}

double FEHillActivationANDInhibition::f_Hill(FEMaterialPoint& pt, const int sol)
{
	double cn = pow(m_psm->GetActualSoluteConcentration(pt, sol), m_n);
	double f_h = (m_B * cn) / (m_Kn + cn);
	return f_h;
}

double FEHillActivationANDInhibition::dfdc(FEMaterialPoint& pt, const int sol)
{
	double c = m_psm->GetActualSoluteConcentration(pt, sol);
	double cn = pow(c, m_n);
	double dfdc = (m_Kn * m_n * f_Hill(pt, sol)) / (c * (m_Kn + cn));
	return dfdc;
}