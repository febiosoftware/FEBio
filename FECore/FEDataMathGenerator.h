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

	bool Init() override;

	// set the math expression
	void setExpression(const std::string& math);

private:
	void value(const vec3d& r, double& data) override;
	void value(const vec3d& r, vec3d& data) override;

private:
	std::string			m_math;
	vector<MSimpleExpression>	m_val;

	DECLARE_FECORE_CLASS()
};
