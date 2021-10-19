/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "DomainDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"

REGISTER_SUPER_CLASS(FELogDomainData, FELOGDOMAINDATA_ID);

//-----------------------------------------------------------------------------
void FEDomainDataRecord::SetData(const char* szexpr)
{
    char szcopy[MAX_STRING] = { 0 };
    strcpy(szcopy, szexpr);
    char* sz = szcopy, * ch;
    m_Data.clear();
    strcpy(m_szdata, szexpr);
    do
    {
        ch = strchr(sz, ';');
        if (ch) *ch++ = 0;
        FELogDomainData* pdata = fecore_new<FELogDomainData>(sz, GetFEModel());
        if (pdata) m_Data.push_back(pdata);
        else throw UnknownDataField(sz);
        sz = ch;
    } while (ch);
}

//-----------------------------------------------------------------------------
FEDomainDataRecord::FEDomainDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_DOMAIN) {}

//-----------------------------------------------------------------------------
int FEDomainDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
double FEDomainDataRecord::Evaluate(int item, int ndata)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int nd = item - 1;
    if ((nd < 0) || (nd >= mesh.Domains())) return 0;

    FEDomain& dom = mesh.Domain(nd);
    return m_Data[ndata]->value(dom);
}

//-----------------------------------------------------------------------------
void FEDomainDataRecord::SetDomain(int domainIndex)
{
    m_item.clear();
    m_item.push_back(domainIndex + 1);
}

//-----------------------------------------------------------------------------
void FEDomainDataRecord::SelectAllItems()
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int n = mesh.Domains();
    m_item.resize(n);
    for (int i = 0; i < n; ++i) m_item[i] = i + 1;
}
