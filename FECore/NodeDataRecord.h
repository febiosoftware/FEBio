#pragma once
#include "FECoreBase.h"
#include "DataRecord.h"

class FENodeSet;

//-----------------------------------------------------------------------------
//! This is the base class for a node data value.
//! \todo I'd like to modify this so I can pass the FENode class instead of the node number
class FECORE_API FENodeLogData : public FECoreBase
{ 
	FECORE_SUPER_CLASS

public:
	FENodeLogData(FEModel* fem);
	virtual ~FENodeLogData();
	virtual double value(int node) = 0; 
};

//-----------------------------------------------------------------------------
//! This class records nodal data
//! \todo should I create a different class for each data record? Like for the plot file?
class FECORE_API NodeDataRecord : public DataRecord
{
public:
	NodeDataRecord(FEModel* pfem, const char* szfile);
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	void SetItemList(FENodeSet* pns);
	int Size() const;

private:
	vector<FENodeLogData*>	m_Data;
};

//-----------------------------------------------------------------------------
// Special class for outputting nodal variables
class FECORE_API FENodeVarData : public FENodeLogData
{
public:
	FENodeVarData(FEModel* pfem, int ndof);
	double value(int node);
private:
	int	m_ndof;
};
