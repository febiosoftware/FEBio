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
#include "DataStore.h"
#include "log.h"
#include "FEModel.h"
#include "FEAnalysis.h"
#include "FELogElemMath.h"
#include "FECoreKernel.h"

//-----------------------------------------------------------------------------
DataStore::DataStore(FEModel* fem) : m_fem(fem)
{
}

//-----------------------------------------------------------------------------
DataStore::~DataStore()
{
}

//-----------------------------------------------------------------------------
void DataStore::Clear()
{
	for (size_t i=0; i<m_data.size(); ++i) delete m_data[i];
	m_data.clear();
}

bool DataStore::Init()
{
	for (size_t i = 0; i < m_elemDefs.size(); ++i)
	{
		auto& def = *m_elemDefs[i];
		if (def.Init() == false) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------

void DataStore::Write()
{
	for (size_t i=0; i<m_data.size(); ++i)
	{
		DataRecord& DR = *m_data[i];
		DR.Write();
	}
}

//-----------------------------------------------------------------------------

void DataStore::AddRecord(DataRecord* prec)
{
	prec->m_nid = (int) m_data.size() + 1;
	m_data.push_back(prec);
}

void DataStore::AddElementDataDefinition(FELogElemDefinition* def)
{
	m_elemDefs.push_back(def);
}

FELogElemDefinition* DataStore::FindElementDataDefinition(const std::string& name)
{
	for (size_t i = 0; i < m_elemDefs.size(); ++i)
	{
		if (m_elemDefs[i]->GetName() == name) return m_elemDefs[i];
	}
	return nullptr;
}

FELogElemSource* DataStore::GetElementDataSource(const std::string& name)
{
	if (name.empty()) return nullptr;

	FELogElemSource* pdata = nullptr;
	if (name[0] == '=')
	{
		FELogElemMath* logMath = fecore_alloc(FELogElemMath, m_fem);
		if (logMath)
		{
			string smath(name.c_str() + 1);
			if (logMath->SetExpression(smath))
			{
				pdata = logMath;
			}
		}
	}
	else
	{
		// try to allocate a built-in element data source
		pdata = fecore_new<FELogElemData>(name.c_str(), m_fem);
		if (pdata == nullptr)
		{
			// see if it's a user-defined element data source
			pdata = FindElementDataDefinition(name);
		}
	}

	return pdata;
}
