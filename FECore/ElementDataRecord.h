#pragma once
#include "DataStore.h"

//-----------------------------------------------------------------------------
//! Base class for element log data
class FELogElemData
{
public:
	FELogElemData(FEModel* pfem) : m_pfem(pfem){}
	virtual ~FELogElemData(){}
	virtual double value(FEElement& el) = 0;
protected:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
class ElementDataRecord : public DataRecord
{
	struct ELEMREF
	{
		int	ndom;
		int	nid;
	};

public:
	ElementDataRecord(FEModel* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	int Size() { return (int) m_Data.size(); }

protected:
	void BuildELT();

protected:
	vector<ELEMREF>	m_ELT;
	vector<FELogElemData*>	m_Data;
};
