#include "stdafx.h"
#include "FEDiscreteMaterial.h"

//-----------------------------------------------------------------------------
// FELinearSpring
//-----------------------------------------------------------------------------

double FELinearSpring::force(double dl)
{
	return m_E*dl;
}

double FELinearSpring::stiffness(double dl)
{
	return m_E;
}

void FELinearSpring::Init()
{
	if (m_E < 0) throw MaterialError("Invalid value for E in linear spring material");
}

//-----------------------------------------------------------------------------
// FETensionOnlyLinearSpring
//-----------------------------------------------------------------------------

double FETensionOnlyLinearSpring::force(double dl)
{
	if (dl >= 0) return m_E*dl; else return 0;
}

double FETensionOnlyLinearSpring::stiffness(double dl)
{
	return (dl >= 0 ? m_E : 0);
}

void FETensionOnlyLinearSpring::Init()
{
	if (m_E < 0) throw MaterialError("Invalid value for E in tension only linear spring material");
}

//-----------------------------------------------------------------------------
// FENonLinearSpring
//-----------------------------------------------------------------------------

void FENonLinearSpring::Init()
{
	if (m_nlc < 0) throw MaterialError("Invalid load curve ID for nonlinear spring");
}

double FENonLinearSpring::force(double dl)
{
	assert(m_plc);
	return m_F*m_plc->Value(dl);
}

double FENonLinearSpring::stiffness(double dl)
{
	assert(m_plc);
	return m_F*m_plc->Deriv(dl);
}
