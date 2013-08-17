#pragma once
#include "DataStore.h"

//-----------------------------------------------------------------------------
//! This is the base class for a node data value.
//! \todo I'd like to modify this so I can pass the FENode class instead of the node number
class FENodeLogData
{ 
public:
	FENodeLogData(FEModel* pfem) : m_pfem(pfem) {}
	virtual ~FENodeLogData(){}
	virtual double value(int node) = 0; 
protected:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! This class records nodal data
//! \todo should I create a different class for each data record? Like for the plot file?
class NodeDataRecord : public DataRecord
{
public:
	NodeDataRecord(FEModel* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	void SetItemList(FENodeSet* pns);
	int Size() { return (int) m_Data.size(); }

private:
	vector<FENodeLogData*>	m_Data;
};
