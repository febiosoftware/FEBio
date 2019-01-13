#pragma once
#include <FECore/FECoreBase.h>
#include <FECore/DataRecord.h>
#include "FERigidBody.h"

//-----------------------------------------------------------------------------
//! Base class for object log data (e.g. rigid bodies)
class FEBIOMECH_API FELogObjectData : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FELogObjectData(FEModel* fem) : FECoreBase(fem, FEOBJLOGDATA_ID) {}
	virtual ~FELogObjectData(){}
	virtual double value(FERigidBody& rb) = 0;
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API ObjectDataRecord : public DataRecord
{
public:
	ObjectDataRecord(FEModel* pfem, const char* szfile);
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	int Size() const;

private:
	vector<FELogObjectData*>	m_Data;
};
