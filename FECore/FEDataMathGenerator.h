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
	void value(const vec3d& r, vector<double>& data) override;

private:
	std::string			m_math;
	MSimpleExpression	m_val[3];

	DECLARE_FECORE_CLASS()
};
