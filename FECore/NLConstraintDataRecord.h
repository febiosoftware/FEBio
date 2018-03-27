#pragma once
#include "FECoreBase.h"
#include "DataRecord.h"
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints log data (e.g. rigid connectors)
class FELogNLConstraintData : public FECoreBase
{
public:
    FELogNLConstraintData(FEModel* pfem) : FECoreBase(FENLCLOGDATA_ID), m_pfem(pfem) {}
    virtual ~FELogNLConstraintData(){}
    virtual double value(FENLConstraint& rc) = 0;
private:
    FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
class NLConstraintDataRecord : public DataRecord
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
