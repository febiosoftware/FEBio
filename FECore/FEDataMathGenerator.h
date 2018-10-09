#pragma once
#include "FEDataGenerator.h"
#include <string>
#include "MathObject.h"

class FENodeSet;
class FEFacetSet;

//-----------------------------------------------------------------------------
class FECORE_API FEDataMathGenerator : public FEDataGenerator
{
public:
	FEDataMathGenerator(FEModel* fem);

	bool Init();

	// set the math expression
	void setExpression(const std::string& math);

private:
	double value(const vec3d& r) override;

private:
	std::string			m_math;
	MSimpleExpression	m_val;

	DECLARE_FECORE_CLASS()
};
