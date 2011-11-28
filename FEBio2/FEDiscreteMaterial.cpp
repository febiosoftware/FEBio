#include "stdafx.h"
#include "FEDiscreteMaterial.h"
#include "fem.h"

//-----------------------------------------------------------------------------
// FELinearSpring
//-----------------------------------------------------------------------------

// register the material with the framework
REGISTER_MATERIAL(FELinearSpring, "linear spring");

// define the material parameters
BEGIN_PARAMETER_LIST(FELinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
END_PARAMETER_LIST();

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

// register the material with the framework
REGISTER_MATERIAL(FETensionOnlyLinearSpring, "tension only linear spring");

// define the material parameters
BEGIN_PARAMETER_LIST(FETensionOnlyLinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
END_PARAMETER_LIST();

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

// register the material with the framework
REGISTER_MATERIAL(FENonLinearSpring, "nonlinear spring");

// define the material parameters
BEGIN_PARAMETER_LIST(FENonLinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER(m_F, FE_PARAM_DOUBLE, "F");
END_PARAMETER_LIST();

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

void FENonLinearSpring::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_F << m_nlc;
	}
	else
	{
		ar >> m_F >> m_nlc;
		m_plc = ar.GetFEModel()->GetLoadCurve(m_nlc);
	}
}
