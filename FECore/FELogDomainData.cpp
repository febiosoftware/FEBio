/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FELogDomainData.h"
#include "FELogElemData.h"
#include "FEDomain.h"
#include "FEMesh.h"

FELogAvgDomainData::FELogAvgDomainData(FEModel* pfem) : FELogDomainData(pfem)
{
	m_elemData = nullptr;
}

FELogAvgDomainData::~FELogAvgDomainData()
{
	if (m_elemData) delete m_elemData;
	m_elemData = nullptr;
}

bool FELogAvgDomainData::SetParameters(const std::vector<std::string>& params)
{
	if (params.size() != 1) return false;
	const std::string& v1 = params[0];
	if (v1.empty()) return false;

	m_elemData = fecore_new<FELogElemData>(v1.c_str(), GetFEModel());
	if (m_elemData == nullptr) return false;

	return true;
}

double FELogAvgDomainData::value(FEDomain& dom)
{
	if (m_elemData == nullptr) return 0.0;

	double avg = 0.0;
	const int NE = dom.Elements();
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double eval = m_elemData->value(el);
		avg += eval;
	}
	avg /= (double)NE;
	return avg;
}

//============================================================================
FELogPctDomainData::FELogPctDomainData(FEModel* pfem) : FELogDomainData(pfem)
{
	m_pct = 0.0;
	m_elemData = nullptr;
}

FELogPctDomainData::~FELogPctDomainData()
{
	if (m_elemData) delete m_elemData;
	m_elemData = nullptr;
}

bool FELogPctDomainData::SetParameters(const std::vector<std::string>& params)
{
	if (params.size() != 2) return false;
	const std::string& v1 = params[0];
	const std::string& v2 = params[1];
	if (v1.empty() || v2.empty()) return false;

	m_elemData = fecore_new<FELogElemData>(v1.c_str(), GetFEModel());
	if (m_elemData == nullptr) return false;

	m_pct = atof(v2.c_str());
	if ((m_pct < 0.0) || (m_pct > 1.0)) return false;

	return true;
}

double FELogPctDomainData::value(FEDomain& dom)
{
	if (m_elemData == nullptr) return 0.0;

	const int NE = dom.Elements();
	vector<double> val(NE, 0.0);
	for (int i = 0; i < NE; ++i)
	{
		FEElement& el = dom.ElementRef(i);
		val[i] = m_elemData->value(el);
	}

	std::sort(val.begin(), val.end());

	int n = (int)(m_pct * ((double)val.size() - 1.0));
	return val[n];
}


//============================================================================
FELogIntegralDomainData::FELogIntegralDomainData(FEModel* pfem) : FELogDomainData(pfem)
{
	m_elemData = nullptr;
}

FELogIntegralDomainData::~FELogIntegralDomainData()
{
	if (m_elemData) delete m_elemData;
	m_elemData = nullptr;
}

bool FELogIntegralDomainData::SetParameters(const std::vector<std::string>& params)
{
	if (params.size() != 1) return false;
	const std::string& v1 = params[0];
	if (v1.empty()) return false;

	m_elemData = fecore_new<FELogElemData>(v1.c_str(), GetFEModel());
	if (m_elemData == nullptr) return false;

	return true;
}

double FELogIntegralDomainData::value(FEDomain& dom)
{
	if (m_elemData == nullptr) return 0.0;

	FEMesh* mesh = dom.GetMesh();

	double sum = 0.0;
	const int NE = dom.Elements();
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double eval = m_elemData->value(el);
		double vol = mesh->ElementVolume(el);
		sum += eval * vol;
	}
	return sum;
}
