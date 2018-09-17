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
	MathObject() {}
	virtual ~MathObject() {}

	int Dim() { return (int)m_Var.size(); }

	void AddVariable(MVariable* pv);
	MVariable* FindVariable(const std::string& s);
	int Variables() { return (int)m_Var.size(); }

	virtual MathObject* copy() = 0;

protected:
	void BuildVarList(MItem* pi);

protected:
	MVarList	m_Var;		// list of variables
};

//-----------------------------------------------------------------------------
// This class defines a simple epxression that can be evaluated by
// setting the values of the variables.
class MSimpleExpression : public MathObject
{
public:
	void SetExpression(MITEM& e) { m_item = e; m_Var.clear(); BuildVarList(m_item.ItemPtr()); }
	MITEM& GetExpression() { return m_item; }

	MathObject* copy()
	{
		MSimpleExpression* po = new MSimpleExpression();
		po->SetExpression(m_item);
		return po;
	}

	double value(MItem* pi);

	double value() { return value(m_item.ItemPtr());  }

	int Items();

protected:
	MITEM	m_item;
};

//-----------------------------------------------------------------------------
// User definition
class MDefinition : public MathObject
{
public:
	MDefinition(){}
	MDefinition(std::string& s, MITEM& e) : m_name(s), m_item(e){}

	void SetName(std::string& s) { m_name = s; }
	void SetExpression(MITEM& e) { m_item = e; }

	std::string GetName() { return m_name; }
	MITEM& GetExpression() { return m_item; }

	MathObject* copy()
	{
		MDefinition* pd = new MDefinition();
		pd->SetName(m_name);
		pd->SetExpression(m_item);
		return pd;
	}

protected:
	std::string	m_name;
	MITEM		m_item;
};

//-----------------------------------------------------------------------------
// Class for defining user functions of N-variables
class MFuncDef : public MathObject
{
public:
	MFuncDef(){}
	~MFuncDef() {}
	MFuncDef(std::string& s) : m_name(s) {}

public:
	void SetName(std::string& s) { m_name = s; }
	std::string GetName() { return m_name; }

	void SetExpression(MITEM& e) { m_item = e; m_Var.clear(); BuildVarList(m_item.ItemPtr()); }
	MITEM& GetExpression() { return m_item; }

	void SetVariables(int n) { m_pv.resize(n); }
	int GetVariables() { return (int) m_pv.size(); }
	void SetVariable(int i, MVariable* pv) { m_pv[i] = pv; }
	MVariable* GetVariable(int i) { return m_pv[i]; }

	MathObject* copy()
	{
		MFuncDef* pd = new MFuncDef(m_name);
		pd->SetExpression(m_item);
		pd->m_pv = m_pv;
		return pd;
	}

	MSFuncND* CreateFunction(MSequence& s);

protected:
	std::string				m_name;
	MITEM					m_item;
	std::vector<MVariable*>	m_pv;
};

//-----------------------------------------------------------------------------
// class that holds a list of mathobjects
class MObjList
{
public:
	MObjList(){}
	~MObjList() { for (int i=0; i<(int)m_obj.size(); ++i) delete m_obj[i]; m_obj.clear(); }

	void Add(MathObject* po) { m_obj.push_back(po); }
	int Objects() { return (int) m_obj.size(); }
	MathObject* Object(int i) { return m_obj[i]; }

	MObjList* copy()
	{
		MObjList* po = new MObjList();
		for (int i=0; i<(int)m_obj.size();++i) po->Add(m_obj[i]->copy());
		return po;
	}

protected:
	std::vector<MathObject*>	m_obj;
};
