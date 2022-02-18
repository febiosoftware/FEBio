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
#pragma once
#include "fecore_api.h"
#include <vector>
#include <string>

class DumpStream;

class FECORE_API FEPlotVariable
{
public:
	FEPlotVariable();
	FEPlotVariable(const FEPlotVariable& pv);
	void operator = (const FEPlotVariable& pv);

	FEPlotVariable(const std::string& var, std::vector<int>& item, const char* szdom = "");

	void Serialize(DumpStream& ar);

	const std::string& Name() const { return m_svar; }
	const std::string& DomainName() const { return m_sdom; }

public:
	std::string			m_svar;		//!< name of output variable
	std::string			m_sdom;		//!< (optional) name of domain
	std::vector<int>	m_item;		//!< (optional) list of items
};

class FECORE_API FEPlotDataStore
{
public:
	FEPlotDataStore();
	FEPlotDataStore(const FEPlotDataStore&);
	void operator = (const FEPlotDataStore&);

	void AddPlotVariable(const char* szvar, std::vector<int>& item, const char* szdom = "");

	int GetPlotCompression() const;
	void SetPlotCompression(int n);

	void SetPlotFileType(const std::string& fileType);

	void Serialize(DumpStream& ar);

	int PlotVariables() const { return (int)m_plot.size(); }
	FEPlotVariable& GetPlotVariable(int n) { return m_plot[n]; }

private:
	std::string					m_splot_type;
	std::vector<FEPlotVariable>	m_plot;
	int							m_nplot_compression;
};
