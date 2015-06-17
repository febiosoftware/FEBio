#include "stdafx.h"
#include "NLConstraintDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"
#include "FEBioMech/FERigidConnector.h"

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
        FELogNLConstraintData* pdata = fecore_new<FELogNLConstraintData>(FENLCLOGDATA_ID, sz, m_pfem);
        if (pdata) m_Data.push_back(pdata);
        else throw UnknownDataField(sz);
        sz = ch;
    }
    while (ch);
}

//-----------------------------------------------------------------------------
double NLConstraintDataRecord::Evaluate(int item, int ndata)
{
    int nrc = item - 1;
    if ((nrc < 0) || (nrc >= m_pfem->NonlinearConstraints())) return 0;
    
    double val = 0;
    
    // find the nonlinear constraint that has this rigid connector
    int NLC = m_pfem->NonlinearConstraints();
    for (int i=0; i<NLC; ++i)
    {
        FENLConstraint& nlc = *m_pfem->NonlinearConstraint(i);
        FERigidConnector* rc = dynamic_cast<FERigidConnector*>(&nlc);
        if (rc && (rc->GetConnectorID() == nrc)) return m_Data[ndata]->value(nlc);
    }
    
    return val;
}

//-----------------------------------------------------------------------------
void NLConstraintDataRecord::SelectAllItems()
{
    int n = 0, i;
    for (i=0; i<m_pfem->NonlinearConstraints(); ++i)
    {
        FENLConstraint* pm = m_pfem->NonlinearConstraint(i);
        FERigidConnector* rc = dynamic_cast<FERigidConnector*>(pm);
        if (rc) ++n;
    }
    
    if (n > 0)
    {
        m_item.resize(n);
        n = 0;
        for (i=0; i<m_pfem->NonlinearConstraints(); ++i)
        {
            FENLConstraint* pm  = m_pfem->NonlinearConstraint(i);
            FERigidConnector* rc = dynamic_cast<FERigidConnector*>(pm);
            if (rc)
            {
                m_item[n++] = i+1;
            }
        }
    }
}
