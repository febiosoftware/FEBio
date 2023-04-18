/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEValuator.h"
#include "MathObject.h"
#include "FEDataMap.h"
#include "FENodeDataMap.h"

//---------------------------------------------------------------------------------------
// Base class for evaluating scalar parameters
class FECORE_API FEScalarValuator : public FEValuator
{
	FECORE_SUPER_CLASS(FESCALARVALUATOR_ID)
	FECORE_BASE_CLASS(FEScalarValuator)

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

	FEScalarValuator* copy() override;

private:
	double	m_val;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
class FEMathExpression : public MSimpleExpression
{
	struct MathParam
	{
		int			type;	// 0 = param, 1 = map
		FEParam* pp;
		FEDataMap* map;
	};

public:
	bool Init(const std::string& expr, FECoreBase* pc = nullptr);

	void operator = (const FEMathExpression& me);

	double value(FEModel* fem, const FEMaterialPoint& pt);

private:
	std::vector<MathParam>	m_vars;
};

//---------------------------------------------------------------------------------------
class FECORE_API FEMathValue : public FEScalarValuator
{
public:
	FEMathValue(FEModel* fem);
	~FEMathValue();
	double operator()(const FEMaterialPoint& pt) override;

	bool Init() override;

	FEScalarValuator* copy() override;

	void setMathString(const std::string& s);

	bool create(FECoreBase* pc = 0);

	void Serialize(DumpStream& ar) override;

private:
	std::string			m_expr;
	FEMathExpression	m_math;
	FECoreBase*			m_parent;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
class FECORE_API FEMappedValue : public FEScalarValuator
{
public:
	FEMappedValue(FEModel* fem);
	void setDataMap(FEDataMap* val);
	void setScaleFactor(double s);

	FEDataMap* dataMap();

	double operator()(const FEMaterialPoint& pt) override;

	FEScalarValuator* copy() override;

	void Serialize(DumpStream& dmp) override;

private:
	double		m_scale;	// scale factor
	FEDataMap*	m_val;
};

//---------------------------------------------------------------------------------------
class FECORE_API FENodeMappedValue : public FEScalarValuator
{
public:
	FENodeMappedValue(FEModel* fem);
	void setDataMap(FENodeDataMap* val);

	double operator()(const FEMaterialPoint& pt) override;

	FEScalarValuator* copy() override;

private:
	FENodeDataMap*		m_val;
};

