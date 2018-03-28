#include "stdafx.h"
#include "FEMathValue.h"
#include "DumpStream.h"

FEMathDouble::FEMathDouble()
{
	m_scale = 1.0;
	setExpression("0.0");
}

void FEMathDouble::setExpression(const std::string& expr)
{
	m_expr = expr;
}

std::string FEMathDouble::getExpression() const
{
	return m_expr;
}

void FEMathDouble::setScale(double s)
{
	m_scale = s;
}

double FEMathDouble::value()
{
	int ierr = 0;
	double v = m_math.eval(m_expr.c_str(), ierr);
	return m_scale*v;
}

void FEMathDouble::setVariable(const char* sz, double v)
{
	m_math.SetVariable(sz, v);
}

void FEMathDouble::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_expr;
		ar << m_scale;
	}
	else
	{
		ar >> m_expr;
		ar >> m_scale;
	}
}
