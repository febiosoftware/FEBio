#include "stdafx.h"
#include "FESpringMaterial.h"

//-----------------------------------------------------------------------------
// FELinearSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FELinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
END_FECORE_CLASS();

double FELinearSpring::force(double dl)
{
	return m_E*dl;
}

double FELinearSpring::stiffness(double dl)
{
	return m_E;
}

//-----------------------------------------------------------------------------
// FETensionOnlyLinearSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FETensionOnlyLinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
END_FECORE_CLASS();

double FETensionOnlyLinearSpring::force(double dl)
{
	if (dl >= 0) return m_E*dl; else return 0;
}

double FETensionOnlyLinearSpring::stiffness(double dl)
{
	return (dl >= 0 ? m_E : 0);
}

//-----------------------------------------------------------------------------
// FENonLinearSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FENonLinearSpring, FEDiscreteMaterial)
	ADD_PROPERTY(m_F, "force");
END_FECORE_CLASS();

FENonLinearSpring::FENonLinearSpring(FEModel* pfem) : FESpringMaterial(pfem)
{
	m_F = nullptr;
}

double FENonLinearSpring::force(double dl)
{
	return m_F->value(dl);
}

double FENonLinearSpring::stiffness(double dl)
{
	return m_F->derive(dl);
}

//-----------------------------------------------------------------------------
// FEExperimentalSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEExperimentalSpring, FESpringMaterial)
	ADD_PARAMETER(m_E , "E");
	ADD_PARAMETER(m_sM, "sM");
	ADD_PARAMETER(m_sm, "sm");
END_FECORE_CLASS();

FEExperimentalSpring::FEExperimentalSpring(FEModel* pfem) : FESpringMaterial(pfem)
{
	m_E = 0.0;
	m_sM = 0.0;
	m_sm = 0.0;
}

double FEExperimentalSpring::force(double dl)
{
	if (dl >= 0.0)
		return m_sM*(1.0 - exp(-m_E*dl / m_sM));
	else
		return -m_sm*(1.0 - exp(m_E*dl / m_sm));
}

double FEExperimentalSpring::stiffness(double dl)
{
	if (dl >= 0.0)
		return m_E*exp(-m_E*dl / m_sM);
	else
		return m_E*exp(m_E*dl / m_sm);
}
