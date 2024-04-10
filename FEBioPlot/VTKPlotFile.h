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
#include "PlotFile.h"
#include <stdio.h>

//! This class stores the FEBio results to a family of VTK files. 
class VTKPlotFile : public PlotFile
{
public:
	VTKPlotFile(FEModel* fem);

	//! Open the plot database
	bool Open(const char* szfile) override;

	//! Open for appending
	bool Append(const char* szfile) override;

	//! Write current FE state to plot database
	bool Write(float ftime, int flag = 0) override;

	//! see if the plot file is valid
	bool IsValid() const override;

	void Serialize(DumpStream& ar) override;

private:
	void WriteHeader();
	void WritePoints();
	void WriteCells();
	void WritePointData();
	void WriteCellData();

	void WriteScalarData(std::vector<float>& val, const std::string& szname);
	void WriteVectorData(std::vector<float>& val, const std::string& szname);
	void WriteMat3FSData(std::vector<float>& val, const std::string& szname);
	void WriteMat3FDData(std::vector<float>& val, const std::string& szname);
	void WriteArrayData (std::vector<float>& val, const std::string& name, FEPlotData* pd);
	void WriteArrayVec3fData(std::vector<float>& val, const std::string& name, FEPlotData* pd);

private:
	FILE*	m_fp;
	int		m_count;
	bool	m_valid;
	std::string	m_filename;
};
