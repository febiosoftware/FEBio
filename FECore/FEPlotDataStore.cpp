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
#include "FEPlotDataStore.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEPlotVariable::FEPlotVariable() {}

//-----------------------------------------------------------------------------
FEPlotVariable::FEPlotVariable(const FEPlotVariable& pv)
{
    m_svar = pv.m_svar;
    m_sdom = pv.m_sdom;
    m_item = pv.m_item;
}

//-----------------------------------------------------------------------------
void FEPlotVariable::operator = (const FEPlotVariable& pv)
{
    m_svar = pv.m_svar;
    m_sdom = pv.m_sdom;
    m_item = pv.m_item;
}

FEPlotVariable::FEPlotVariable(const std::string& var, std::vector<int>& item, const char* szdom)
{
    m_svar = var;
    if (szdom) m_sdom = szdom;
    m_item = item;
}

void FEPlotVariable::Serialize(DumpStream& ar)
{
    ar & m_svar;
    ar & m_sdom;
    ar & m_item;
}

//=======================================================================================
FEPlotDataStore::FEPlotDataStore()
{
    m_plot.clear();
    m_nplot_compression = 0;
}

//-----------------------------------------------------------------------------
FEPlotDataStore::FEPlotDataStore(const FEPlotDataStore& plt)
{
    m_splot_type = plt.m_splot_type;
    m_nplot_compression = plt.m_nplot_compression;
    m_plot = plt.m_plot;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::operator = (const FEPlotDataStore& plt)
{
    m_splot_type = plt.m_splot_type;
    m_nplot_compression = plt.m_nplot_compression;
    m_plot = plt.m_plot;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::AddPlotVariable(const char* szvar, std::vector<int>& item, const char* szdom)
{
    FEPlotVariable var(szvar, item, szdom);
    m_plot.push_back(var);
}

//-----------------------------------------------------------------------------
int FEPlotDataStore::GetPlotCompression() const
{
    return m_nplot_compression;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::SetPlotCompression(int n)
{
    m_nplot_compression = n;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::SetPlotFileType(const std::string& fileType)
{
    m_splot_type = fileType;
}

void FEPlotDataStore::Serialize(DumpStream& ar)
{
    ar & m_nplot_compression;
    ar & m_splot_type;
    ar & m_plot;
}
