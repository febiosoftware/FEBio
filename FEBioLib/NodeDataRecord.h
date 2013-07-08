#pragma once
#include "DataStore.h"

//-----------------------------------------------------------------------------
//! This class records nodal data
//! \todo should I create a different class for each data record? Like for the plot file?
class NodeDataRecord : public DataRecord
{
	enum { X, Y, Z, UX, UY, UZ, VX, VY, VZ, RX, RY, RZ, T, P, C, C1, C2, C3, C4, C5, C6 };

public:
	NodeDataRecord(FEModel* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	void SetItemList(FENodeSet* pns);
};
