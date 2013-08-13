#pragma once
#include "FECore/DataStore.h"

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
class FENodeXPos : public FENodeLogData
{ 
public: 
	FENodeXPos(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYPos : public FENodeLogData 
{ 
public: 
	FENodeYPos(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZPos : public FENodeLogData
{ 
public: 
	FENodeZPos(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeXDisp : public FENodeLogData
{ 
public: 
	FENodeXDisp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYDisp : public FENodeLogData
{ 
public: 
	FENodeYDisp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZDisp : public FENodeLogData
{ 
public: 
	FENodeZDisp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeXVel : public FENodeLogData
{ 
public: 
	FENodeXVel(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYVel : public FENodeLogData
{ 
public: 
	FENodeYVel(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZVel : public FENodeLogData
{ 
public: 
	FENodeZVel(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeTemp : public FENodeLogData
{ 
public: 
	FENodeTemp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodePressure : public FENodeLogData
{ 
public: 
	FENodePressure(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeConcentration : public FENodeLogData
{ 
public: 
	FENodeConcentration(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeConcentration_ : public FENodeLogData
{ 
protected: 
	FENodeConcentration_(FEModel* pfem, int nsol) : FENodeLogData(pfem), m_nsol(nsol) {} 
	double value(int node); 
private:
	int	m_nsol;
};

//-----------------------------------------------------------------------------
class FENodeForceX: public FENodeLogData
{ 
public: 
	FENodeForceX(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceY: public FENodeLogData
{ 
public: 
	FENodeForceY(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceZ: public FENodeLogData
{ 
public: 
	FENodeForceZ(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
template <int N> class FENodeConcentration_T : public FENodeConcentration_
{ 
public: 
	FENodeConcentration_T(FEModel* pfem) : FENodeConcentration_(pfem, N) {} 
	double value(int node) { return FENodeConcentration_::value(node); }
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
