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
#include "PlotFile.h"
#include <FECore/FEPlotDataStore.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
PlotFile::PlotFile(FEModel* fem) : m_pfem(fem)
{

}

//-----------------------------------------------------------------------------
PlotFile::~PlotFile()
{
	Close();

	// clear all arrays
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Glob.begin();
	for (int i = 0; i < (int)m_dic.m_Glob.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Mat.begin();
	for (int i = 0; i < (int)m_dic.m_Mat.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Node.begin();
	for (int i = 0; i < (int)m_dic.m_Node.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Elem.begin();
	for (int i = 0; i < (int)m_dic.m_Elem.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Face.begin();
	for (int i = 0; i < (int)m_dic.m_Face.size(); ++i, ++it) delete it->m_psave;
}

//-----------------------------------------------------------------------------
void PlotFile::Close()
{
	
}

//-----------------------------------------------------------------------------
bool PlotFile::AddVariable(FEPlotData* ps, const char* szname)
{
	vector<int> dummy;
	switch (ps->RegionType())
	{
	case FE_REGION_NODE: return m_dic.AddNodalVariable(ps, szname, dummy);
	case FE_REGION_DOMAIN: return m_dic.AddDomainVariable(ps, szname, dummy);
	case FE_REGION_SURFACE: return m_dic.AddSurfaceVariable(ps, szname, dummy);
	default:
		assert(false);
		return false;
	}
}

//-----------------------------------------------------------------------------
bool PlotFile::AddVariable(const char* sz)
{
	vector<int> dummy;
	return AddVariable(sz, dummy);
}

//-----------------------------------------------------------------------------
bool PlotFile::AddVariable(const char* sz, vector<int>& item, const char* szdom)
{
	return m_dic.AddVariable(GetFEModel(), sz, item, szdom);
}

//-----------------------------------------------------------------------------
// build the dictionary
void PlotFile::BuildDictionary()
{
	FEPlotDataStore& pltData = GetFEModel()->GetPlotDataStore();
	for (int n = 0; n < pltData.PlotVariables(); ++n)
	{
		FEPlotVariable& vi = pltData.GetPlotVariable(n);
		const std::string& varName = vi.Name();
		const std::string& domName = vi.DomainName();

		// add the plot output variable
		if (AddVariable(varName.c_str(), vi.m_item, domName.c_str()) == false)
		{
			feLog("FATAL ERROR: Output variable \"%s\" is not defined\n", varName.c_str());
			throw "FATAL ERROR";
		}
	}
}
