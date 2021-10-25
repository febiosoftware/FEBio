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
#include "MItem.h"
#include <vector>
#include "fecore_api.h"

//-----------------------------------------------------------------------------
typedef std::vector<MVariable*>	MVarList;

//-----------------------------------------------------------------------------
// This class defines the base class for all math objects
// It also stores a list of all the variables
class FECORE_API MathObject
{
public:
	MathObject();
	MathObject(const MathObject& mo);
	void operator = (const MathObject& mo);

	virtual ~MathObject();

	int Dim() { return (int)m_Var.size(); }

	MVariable* AddVariable(const std::string& var, double initVal = 0.0);
	void AddVariables(const std::vector<std::string>& varList);

	void AddVariable(MVariable* pv);
	MVariable* FindVariable(const std::string& s);
	int Variables() const { return (int)m_Var.size(); }

	MVariable* Variable(int i) { return m_Var[i]; }
	const MVariable* Variable(int i) const { return m_Var[i]; }

	virtual void Clear();

	virtual MathObject* copy() = 0;

protected:
	MVarList	m_Var;		// list of variables
};

//-----------------------------------------------------------------------------
// This class defines a simple epxression that can be evaluated by
// setting the values of the variables.
class FECORE_API MSimpleExpression : public MathObject
{
public:
	MSimpleExpression() {}
	MSimpleExpression(const MSimpleExpression& mo);
	void operator = (const MSimpleExpression& mo);

	void SetExpression(MITEM& e) { m_item = e; }
	MITEM& GetExpression() { return m_item; }
	const MITEM& GetExpression() const { return m_item; }

	// Create a simple expression object from a string
	bool Create(const std::string& expr, bool autoVars = false);

	// copy the expression
	MathObject* copy() { return new MSimpleExpression(*this); }

	// These functions are not thread safe since variable values can be overridden by different threads
	// In multithreaded applications, use the thread safe functions below.
	double value() const { return value(m_item.ItemPtr());  }

	// combines Create and value. Not efficient usage! 
	double value(const std::string& s);

	// This is a thread safe function to evaluate the expression
	// The values of the variables are passed as an argument. This function
	// does not call MVariable->value, but uses these passed values insteads.
	// Make sure that the var array has the same size as the variable array of the expression
	double value_s(const std::vector<double>& var) const
	{ 
		assert(var.size() == m_Var.size());
		return value(m_item.ItemPtr(), var); 
	}

	int Items();

protected:
	double value(const MItem* pi) const;
	double value(const MItem* pi, const std::vector<double>& var) const;

protected:
	void fixVariableRefs(MItem* pi);

protected:
	MITEM	m_item;
};
