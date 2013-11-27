#pragma once
#include "DataStore.h"
#include "FEObject.h"

//-----------------------------------------------------------------------------
//! Base class for object log data (e.g. rigid bodies)
class FELogObjectData
{
public:
	FELogObjectData(FEModel* pfem) : m_pfem(pfem) {}
	virtual ~FELogObjectData(){}
	virtual double value(FEObject& rb) = 0;
private:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
class ObjectDataRecord : public DataRecord
{
public:
	ObjectDataRecord(FEModel* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	int Size() { return (int) m_Data.size(); }

private:
	vector<FELogObjectData*>	m_Data;
};
