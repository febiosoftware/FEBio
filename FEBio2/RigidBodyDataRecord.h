#pragma once
#include "DataStore.h"

//-----------------------------------------------------------------------------
// TODO: should I create a different class for each data record? Like for the plot file?
class RigidBodyDataRecord : public DataRecord
{
	enum {X, Y, Z, QX, QY, QZ, QW, FX, FY, FZ, MX, MY, MZ};

public:
	RigidBodyDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
};
