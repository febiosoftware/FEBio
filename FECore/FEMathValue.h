#pragma once
#include <string>
#include "MathParser.h"

class DumpStream;

class FEMathDouble
{
public:
	FEMathDouble();

	void setExpression(const std::string& expr);
	std::string getExpression() const;

	void setScale(double s);

	double value();

	void setVariable(const char* sz, double v);

public:
	void Serialize(DumpStream& ar);

private:
	std::string	m_expr;		// math expression
	double		m_scale;	// scale value
	MathParser	m_math;		// the math parser
};
