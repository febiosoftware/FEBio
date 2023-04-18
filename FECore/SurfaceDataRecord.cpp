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
#include "SurfaceDataRecord.h"
#include "FECoreKernel.h"
#include "FEModel.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
void FESurfaceDataRecord::SetData(const char* szexpr)
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
        FELogSurfaceData* pdata = fecore_new<FELogSurfaceData>(sz, GetFEModel());
        if (pdata) m_Data.push_back(pdata);
        else throw UnknownDataField(sz);
        sz = ch;
    } while (ch);
}

//-----------------------------------------------------------------------------
FESurfaceDataRecord::FESurfaceDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_SURFACE) {}

//-----------------------------------------------------------------------------
int FESurfaceDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
double FESurfaceDataRecord::Evaluate(int item, int ndata)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int nd = item - 1;
    if ((nd < 0) || (nd >= mesh.Surfaces())) return 0;

    FESurface& surf = mesh.Surface(nd);
    return m_Data[ndata]->value(surf);
}

//-----------------------------------------------------------------------------
void FESurfaceDataRecord::SetSurface(int surfIndex)
{
    m_item.clear();
    m_item.push_back(surfIndex + 1);
}

//-----------------------------------------------------------------------------
void FESurfaceDataRecord::SelectAllItems()
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int n = mesh.Surfaces();
    m_item.resize(n);
    for (int i = 0; i < n; ++i) m_item[i] = i + 1;
}
