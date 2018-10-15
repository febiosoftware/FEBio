#pragma once
#include "FECoreBase.h"
#include "DataRecord.h"
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints log data (e.g. rigid connectors)
class FECORE_API FELogNLConstraintData : public FECoreBase
{
	DECLARE_SUPER_CLASS(FENLCLOGDATA_ID);

public:
    FELogNLConstraintData(FEModel* fem) : FECoreBase(fem, FENLCLOGDATA_ID) {}
    virtual ~FELogNLConstraintData(){}
    virtual double value(FENLConstraint& rc) = 0;
};

//-----------------------------------------------------------------------------
class FECORE_API NLConstraintDataRecord : public DataRecord
{
public:
	NLConstraintDataRecord(FEModel* pfem, const char* szfile) : DataRecord(pfem, szfile, FE_DATA_NLC){}
    double Evaluate(int item, int ndata);
    void Parse(const char* sz);
    void SelectAllItems();
    int Size() { return (int) m_Data.size(); }
    
private:
    vector<FELogNLConstraintData*>	m_Data;
};
