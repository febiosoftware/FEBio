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
#include "FELogElemMath.h"
#include "MObjBuilder.h"

FELogElemMath::FELogElemMath(FEModel* pfem) : FELogElemData(pfem) 
{

}

FELogElemMath::~FELogElemMath()
{
	Clear();
}

void FELogElemMath::Clear()
{
	for (FELogElemData* d : m_data) delete d;
	m_data.clear();
}

double FELogElemMath::value(FEElement& el)
{
	std::vector<double> val(m_data.size());
	for (size_t i = 0; i < m_data.size(); ++i) val[i] = m_data[i]->value(el);
	return m.value_s(val);
}

bool FELogElemMath::SetExpression(const std::string& smath)
{
	Clear();

	MObjBuilder o;
	o.setAutoVars(true);
	if (!o.Create(&m, smath, false)) return false;

	int nvar = m.Variables();
	for (int i = 0; i < nvar; ++i)
	{
		MVariable& var = *m.Variable(i);
		string varName = var.Name();
		FELogElemData* pd = fecore_new<FELogElemData>(varName.c_str(), GetFEModel());
		if (pd == nullptr) return false;
		m_data.push_back(pd);
	}
	return true;
}
