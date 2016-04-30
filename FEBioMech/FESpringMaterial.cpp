#include "stdafx.h"
#include "FESpringMaterial.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// FELinearSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FELinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER2(m_E, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E");
END_PARAMETER_LIST();

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
BEGIN_PARAMETER_LIST(FETensionOnlyLinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER2(m_E, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E");
END_PARAMETER_LIST();

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
BEGIN_PARAMETER_LIST(FENonLinearSpring, FEDiscreteMaterial)
	ADD_PARAMETER(m_F, FE_PARAM_DOUBLE, "force");
END_PARAMETER_LIST();

FENonLinearSpring::FENonLinearSpring(FEModel* pfem) : FESpringMaterial(pfem)
{
	m_nlc = -1; 
	m_plc = 0;
	m_F = 1; 
}

bool FENonLinearSpring::Init()
{
	if (m_nlc < 0) return MaterialError("Invalid load curve ID for nonlinear spring");
	if (m_plc == 0)
	{
		m_plc = GetFEModel()->GetLoadCurve(m_nlc);
		if (m_plc == 0) return false;
	}
	return FEDiscreteMaterial::Init();
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

void FENonLinearSpring::Serialize(DumpStream& ar)
{
	FESpringMaterial::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_F << m_nlc;
	}
	else
	{
		ar >> m_F >> m_nlc;
		m_plc = ar.GetFEModel().GetLoadCurve(m_nlc);
	}
}

bool FENonLinearSpring::SetParameterAttribute(FEParam& p, const char* szatt, const char* szval)
{
	if (strcmp(p.name(), "force") == 0)
	{
		if (strcmp(szatt, "lc") == 0)
		{
			m_nlc = atoi(szval) - 1;
			return true;
		}
	}
	return false;
}
