#pragma once
#include "MItem.h"
#include <vector>

//-----------------------------------------------------------------------------
typedef std::vector<MVariable*>	MVarList;

//-----------------------------------------------------------------------------
// This class defines the base class for all math objects
// It also stores a list of all the variables
class MathObject
{
public:
	MathObject();
	virtual ~MathObject();

	int Dim() { return (int)m_Var.size(); }

	MVariable* AddVariable(const std::string& var);

	void AddVariable(MVariable* pv);
	MVariable* FindVariable(const std::string& s);
	int Variables() { return (int)m_Var.size(); }

	MVariable* Variable(int i) { return m_Var[i]; }

	virtual MathObject* copy() = 0;

protected:
	MVarList	m_Var;		// list of variables
};

//-----------------------------------------------------------------------------
// This class defines a simple epxression that can be evaluated by
// setting the values of the variables.
class MSimpleExpression : public MathObject
{
public:
	void SetExpression(MITEM& e) { m_item = e; }
	MITEM& GetExpression() { return m_item; }
	const MITEM& GetExpression() const { return m_item; }

	// Create a simple expression object from a string
	bool Create(const std::string& expr, bool autoVars = false);

	MathObject* copy()
	{
		MSimpleExpression* po = new MSimpleExpression();
		po->SetExpression(m_item);
		po->m_Var = m_Var;
		return po;
	}

	// These functions are not thread safe since variable values can be overridden by different threads
	// In multithreaded applications, use the thread safe functions below.
	double value() { return value(m_item.ItemPtr());  }

	// This is a thread safe function to evaluate the expression
	// The values of the variables are passed as an argument. This function
	// does not call MVariable->value, but uses these passed values insteads.
	// Make sure that the var array has the same size as the variable array of the expression
	double value_s(const std::vector<double>& var)
	{ 
		assert(var.size() == m_Var.size());
		return value(m_item.ItemPtr(), var); 
	}

	int Items();

protected:
	double value(const MItem* pi);
	double value(const MItem* pi, const std::vector<double>& var);

protected:
	MITEM	m_item;
};
