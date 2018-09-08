#pragma once
#include "FEMaterial.h"
#include "MathParser.h"
#include "FESolidDomain.h"

//---------------------------------------------------------------------------------------
// class for evaluating material parameters
class FEMatValuator
{
public:
	FEMatValuator() {}
	virtual ~FEMatValuator() {}

	virtual double eval(const FEMaterialPoint& pt) = 0;
};

//---------------------------------------------------------------------------------------
class FEMatConstValue : public FEMatValuator
{
public:
	FEMatConstValue(double v = 0.0) : m_val(v) {};
	double eval(const FEMaterialPoint& pt) override { return m_val; }

private:
	double	m_val;
};

//---------------------------------------------------------------------------------------
class FEMatExpression : public FEMatValuator
{
public:
	FEMatExpression(const std::string& s);
	double eval(const FEMaterialPoint& pt) override;

private:
	std::string		m_expr;
	MathParser		m_math;
};

//---------------------------------------------------------------------------------------
class FEMatMappedValue : public FEMatValuator
{
public:
	FEMatMappedValue(FEDomain* dom, std::vector<double>& values);

	double eval(const FEMaterialPoint& pt) override;

private:
	FEDomain*		m_dom;
	std::vector<double>	m_val;
};

//---------------------------------------------------------------------------------------
// This class represents a material parameter.
// NOTE: Work in progress!
class FEMaterialParam
{
public:
	FEMaterialParam();

	// set the value
	void setValue(double v);

	// set the valuator
	void setValuator(FEMatValuator* val);

	// evaluate the parameter at a material point
	double eval(const FEMaterialPoint& pt) { return m_val->eval(pt); }

private:
	FEMatValuator*	m_val;
};
