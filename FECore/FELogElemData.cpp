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
#include "FELogElemData.h"
#include "FEModel.h"

FELogElemData::FELogElemData(FEModel* fem) : FELogElemSource(fem) {}

FELogElemData::~FELogElemData() {}

BEGIN_FECORE_CLASS(FELogElemAlias, FELogElemDefinition)
	ADD_PARAMETER(m_var, "var");
END_FECORE_CLASS();

FELogElemAlias::FELogElemAlias(FEModel* fem) : FELogElemDefinition(fem) {}

bool FELogElemAlias::Init()
{
	if (m_var.empty()) return false;
	m_pdata = fecore_new<FELogElemData>(m_var.c_str(), GetFEModel());
	return (m_pdata != nullptr);
}

double FELogElemAlias::value(FEElement& el)
{
	if (m_pdata) return m_pdata->value(el);
	return 0.0;
}

BEGIN_FECORE_CLASS(FELogElemFunction, FELogElemDefinition)
	ADD_PARAMETER(m_var, "var");
	ADD_PROPERTY(m_func, "function");
END_FECORE_CLASS();

FELogElemFunction::FELogElemFunction(FEModel* fem) : FELogElemDefinition(fem) {}

bool FELogElemFunction::Init()
{
	DataStore& DS = GetFEModel()->GetDataStore();
	m_pdata = DS.GetElementDataSource(m_var);
	if (m_pdata == nullptr) return false;

	return FELogElemDefinition::Init();
}

double FELogElemFunction::value(FEElement& el)
{
	if (m_pdata == nullptr) return 0.0;
	double val = m_pdata->value(el);
	if (m_func) return m_func->value(val);
	return val;
}
