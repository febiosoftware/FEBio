#pragma once
#include "FEValuator.h"
#include "MathObject.h"
#include "FEDataMap.h"
#include "FENodeDataMap.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating scalar parameters
class FECORE_API FEScalarValuator : public FEValuator
{
	FECORE_SUPER_CLASS

public:
	FEScalarValuator(FEModel* fem) : FEValuator(fem) {};

	virtual double operator()(const FEMaterialPoint& pt) = 0;

	virtual FEScalarValuator* copy() = 0;

	virtual bool isConst() { return false; }

	virtual double* constValue() { return nullptr; }
};

//---------------------------------------------------------------------------------------
class FECORE_API FEConstValue : public FEScalarValuator
{
public:
	FEConstValue(FEModel* fem) : FEScalarValuator(fem), m_val(0.0) {};
	double operator()(const FEMaterialPoint& pt) override { return m_val; }

	bool isConst() override { return true; }

	double* constValue() override { return &m_val; }

	FEScalarValuator* copy() override
	{ 
		FEConstValue* val = new FEConstValue(GetFEModel()); 
		val->m_val = m_val;
		return val;
	}

private:
	double	m_val;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
class FECORE_API FEMathValue : public FEScalarValuator
{
public:
	FEMathValue(FEModel* fem) : FEScalarValuator(fem) {}
	~FEMathValue();
	double operator()(const FEMaterialPoint& pt) override;

	bool Init() override;

	FEScalarValuator* copy() override;

	void setMathString(const std::string& s);

	bool create(FECoreBase* pc = 0);

	bool isConst() override { return false; }
	double* constValue() override { return nullptr; }

	void Serialize(DumpStream& ar) override;

private:
	std::string			m_expr;
	MSimpleExpression	m_math;
	std::vector<FEParam*>	m_vars;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
class FECORE_API FEMappedValue : public FEScalarValuator
{
public:
	FEMappedValue(FEModel* fem);
	void setDataMap(FEDataMap* val, double scl = 1.0);

	double operator()(const FEMaterialPoint& pt) override;

	FEScalarValuator* copy() override;

	bool isConst() override { return false; }
	double* constValue() override { return &m_scale; }

private:
	double		m_scale;
	FEDataMap*	m_val;
};

//---------------------------------------------------------------------------------------
class FECORE_API FENodeMappedValue : public FEScalarValuator
{
public:
	FENodeMappedValue(FEModel* fem);
	void setDataMap(FENodeDataMap* val, double scale = 1.0);

	double operator()(const FEMaterialPoint& pt) override;

	FEScalarValuator* copy() override;

	bool isConst() override { return false; }
	double* constValue() override { return &m_scale; }

private:
	double				m_scale;
	FENodeDataMap*		m_val;
};

