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
#include "DataRecord.h"
#include "DumpStream.h"
#include "FEModel.h"
#include "FEAnalysis.h"
#include "log.h"
#include <sstream>

UnknownDataField::UnknownDataField(const std::string& msg) : std::runtime_error(msg)
{
}

DataRecord::DataRecord(FEModel* pfem, int ntype) : FECoreBase(pfem), m_type(ntype)
{
	m_nid = 0;
	m_delim = " ";
	m_bcomm = true;
	m_fp = nullptr;
}

//-----------------------------------------------------------------------------
bool DataRecord::SetFileName(const char* szfile)
{
	if (szfile == nullptr) return false;

	m_filename = szfile;
	m_fp = fopen(szfile, "wt");
	if (m_fp == 0)
	{
		feLogError("FAILED CREATING DATA FILE %s\n\n", szfile);
		return false;
	}

	return true;
}

DataRecord::~DataRecord()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = nullptr;
	}
}

void DataRecord::SetName(const char* sz)
{
	m_name = sz;
}

void DataRecord::SetDelim(const char* sz)
{
	m_delim = sz;
}

void DataRecord::SetFormat(const char* sz)
{
	m_fmt = sz;
}

//-----------------------------------------------------------------------------
bool DataRecord::Initialize()
{
	if (m_item.empty()) SelectAllItems();
	return true;
}

//-----------------------------------------------------------------------------
std::string DataRecord::printToString(int i)
{
	std::stringstream ss;
	ss.precision(12);

	ss << m_item[i] << m_delim;
	int nd = Size();
	for (int j = 0; j<nd; ++j)
	{
		double val = Evaluate(m_item[i], j);
		ss << val;
		if (j != nd - 1) ss << m_delim;
		else ss << "\n";
	}

	return ss.str();
}

//-----------------------------------------------------------------------------
std::string DataRecord::printToFormatString(int i)
{
	int ndata = Size();
	string fmt = m_fmt;

	std::stringstream ss;

	int nitem = m_item[i];
	char* sz = fmt.data(), * ch = 0;
	int j = 0;
	do
	{
		ch = strchr(sz, '%');
		if (ch)
		{
			if (ch[1] == 'i')
			{
				*ch = 0;
				ss << sz;
				*ch = '%'; sz = ch + 2;
				ss << nitem;
			}
			else if (ch[1] == 'l')
			{
				*ch = 0;
				ss << sz;
				*ch = '%'; sz = ch + 2;
				ss << (int)i + 1;
			}
			else if (ch[1] == 'g')
			{
				*ch = 0;
				ss << sz;
				*ch = '%'; sz = ch + 2;
				if (j<ndata)
				{
					double val = Evaluate(nitem, j++);
					ss << val;
				}
			}
			else if (ch[1] == 't')
			{
				*ch = 0;
				ss << sz;
				*ch = '%'; sz = ch + 2;
				ss << "\t";
			}
			else if (ch[1] == 'n')
			{
				*ch = 0;
				ss << "%s";
				*ch = '%'; sz = ch + 2;
				ss << "\n";
			}
			else
			{
				*ch = 0;
				ss << sz;
				*ch = '%'; sz = ch + 1;
			}
		}
		else { ss << sz; break; }
	} while (*sz);
	ss << "\n";

	return ss.str();
}

bool DataRecord::Write()
{
	FEModel* fem = GetFEModel();
	int nstep = fem->GetCurrentStep()->m_ntimesteps;
	double ftime = fem->GetCurrentTime();

	// make a note in the log file
	feLog("\nData Record #%d\n", m_nid);
	feLog("===========================================================================\n");
	feLog("Step = %d\n", nstep);
	feLog("Time = %.9lg\n", ftime);
	feLog("Data = %s\n", m_name.c_str());

	// write some comments
	FILE* fp = m_fp;
	if (fp && m_bcomm)
	{
		// we save the data in a seperate file
		feLog("File = %s\n", m_filename.c_str());

		// make a note in the data file
		fprintf(fp,"*Step  = %d\n", nstep);
		fprintf(fp,"*Time  = %.9lg\n", ftime);
		fprintf(fp,"*Data  = %s\n", m_name.c_str());
	}

	// save the data
	if (m_fmt.empty())
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			std::string out = printToString((int)i);

			if (fp) fprintf(fp, "%s", out.c_str());
			else feLog(out.c_str(),"");
		}
	}
	else
	{
		// print using the format string
		for (size_t i=0; i<m_item.size(); ++i)
		{
			std::string out = printToFormatString((int)i);

			if (fp) fprintf(fp, "%s", out.c_str());
			else feLog(out.c_str(),"");
		}
	}

	if (fp) fflush(fp);

	return true;
}

void DataRecord::SetItemList(const std::vector<int>& items)
{
	m_item = items;
}

void DataRecord::SetItemList(FEItemList* items, const std::vector<int>& selection)
{
	// derived classes should override this
	assert(false);
}

void DataRecord::Serialize(DumpStream &ar)
{
	if (ar.IsShallow()) return;

	// serialize data
	ar & m_nid;
	ar & m_name;
	ar & m_delim;
	ar & m_filename;
	ar & m_bcomm;
	ar & m_item;
	ar & m_data;

	// when we're loading we need to reinitialize the file
	if (ar.IsLoading())
	{
		SetData(m_data.c_str());

		if (m_fp) fclose(m_fp);
		m_fp = 0;
		if (m_filename[0] != 0)
		{
			// reopen data file for appending
			m_fp = fopen(m_filename.c_str(), "a+");
		}
	}
}

std::vector<DataRecordItem> ProcessDataString(const char* szdata)
{
	std::vector<DataRecordItem> data;
	if ((szdata == nullptr) || (szdata[0] == 0)) return data;

	std::string s = szdata;
	while (!s.empty())
	{
		size_t pos = s.find(";");
		std::string name;
		if (pos != std::string::npos)
		{
			name = s.substr(0, pos);
			s.erase(0, pos + strlen(";"));
		}
		else
		{
			name = s;
			s.clear();
		}

		DataRecordItem item;

		// see if parameters are defined
		// TODO: This only processes one parameter. We need to implement a more robust parser that can handle multiple parameters.
		size_t cl = name.find("(");
		if (cl != std::string::npos)
		{
			size_t cr = name.rfind(")");
			if (cr == std::string::npos) throw UnknownDataField(name);

			string params = name.substr(cl + 1, cr - cl - 1);
			name = name.substr(0, cl);

			cl = params.find("'"); if (cl == std::string::npos) throw UnknownDataField(name);
			cr = params.rfind("'"); if (cr == std::string::npos) throw UnknownDataField(name);
			params = params.substr(cl + 1, cr - cl - 1);

			item.params.push_back(params);
		}
		item.name = name;

		data.push_back(item);
		
	}
	return data;
}
