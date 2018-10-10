#include "stdafx.h"
#include "NLConstraintDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
void NLConstraintDataRecord::Parse(const char* szexpr)
{
    char szcopy[MAX_STRING] = {0};
    strcpy(szcopy, szexpr);
    char* sz = szcopy, *ch;
    m_Data.clear();
    strcpy(m_szdata, szexpr);
    do
    {
        ch = strchr(sz, ';');
        if (ch) *ch++ = 0;
        FELogNLConstraintData* pdata = fecore_new<FELogNLConstraintData>(sz, m_pfem);
        if (pdata) m_Data.push_back(pdata);
        else throw UnknownDataField(sz);
        sz = ch;
    }
    while (ch);
}

//-----------------------------------------------------------------------------
double NLConstraintDataRecord::Evaluate(int item, int ndata)
{
    int nc = item - 1;
    if ((nc < 0) || (nc >= m_pfem->NonlinearConstraints())) return 0;
    
	FENLConstraint& nlc = *m_pfem->NonlinearConstraint(nc);
	return m_Data[ndata]->value(nlc);
}

//-----------------------------------------------------------------------------
void NLConstraintDataRecord::SelectAllItems()
{
	int n = m_pfem->NonlinearConstraints();
	m_item.resize(n);
	for (int i = 0; i<n; ++i) m_item[i] = i + 1;
}
