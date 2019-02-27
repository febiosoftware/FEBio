#pragma once
#include "FECoreBase.h"
#include "DataRecord.h"
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints log data (e.g. rigid connectors)
class FECORE_API FELogNLConstraintData : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
    FELogNLConstraintData(FEModel* fem) : FECoreBase(fem) {}
    virtual ~FELogNLConstraintData(){}
    virtual double value(FENLConstraint& rc) = 0;
};

//-----------------------------------------------------------------------------
class FECORE_API NLConstraintDataRecord : public DataRecord
{
public:
	NLConstraintDataRecord(FEModel* pfem, const char* szfile);
    double Evaluate(int item, int ndata);
    void Parse(const char* sz);
    void SelectAllItems();
	int Size() const;
    
private:
    vector<FELogNLConstraintData*>	m_Data;
};
