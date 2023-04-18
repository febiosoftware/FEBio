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
#include <stdio.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "FECoreBase.h"
#include "fecore_api.h"

//-----------------------------------------------------------------------------
// forward declaration
class FEModel;
class DumpStream;
class FEItemList;

//-----------------------------------------------------------------------------
enum FEDataRecordType {
	FE_DATA_NODE = 0x01,
	FE_DATA_FACE,
	FE_DATA_ELEM,
	FE_DATA_RB,
	FE_DATA_NLC,
	FE_DATA_SURFACE,
	FE_DATA_DOMAIN,
	FE_DATA_MODEL
};

//-----------------------------------------------------------------------------
// Exception thrown when parsing fails
class FECORE_API UnknownDataField : public std::runtime_error
{
public:
	UnknownDataField(const char* sz);
};

//-----------------------------------------------------------------------------

class FECORE_API DataRecord : public FECoreBase
{
	FECORE_SUPER_CLASS(FEDATARECORD_ID)
	FECORE_BASE_CLASS(DataRecord)

public:
	enum {MAX_DELIM=16, MAX_STRING=1024};
public:
	DataRecord(FEModel* pfem, int ntype);
	virtual ~DataRecord();

	bool SetFileName(const char* szfile);

	bool Write();

	void SetItemList(const std::vector<int>& items);
	virtual void SetItemList(FEItemList* itemList, const std::vector<int>& selection);

	void SetName(const char* sz);
	void SetDelim(const char* sz);
	void SetFormat(const char* sz);
	void SetComments(bool b) { m_bcomm = b; }

public:
	virtual bool Initialize();
	virtual double Evaluate(int item, int ndata) = 0;
	virtual void SelectAllItems() = 0;
	virtual void Serialize(DumpStream& ar);
	virtual void SetData(const char* sz) = 0;
	virtual int Size() const = 0;

private:
	std::string printToString(int i);
	std::string printToFormatString(int i);

public:
	int					m_nid;		//!< ID of data record
	std::vector<int>	m_item;		//!< item list
	int					m_type;		//!< type of data record

protected:
	bool	m_bcomm;				//!< export comments or not
	char	m_szname[MAX_STRING];	//!< name of expression
	char	m_szdelim[MAX_DELIM];	//!< data delimitor
	char	m_szdata[MAX_STRING];	//!< data expression
	char	m_szfmt[MAX_STRING];	//!< max format string

protected:
	char	m_szfile[MAX_STRING];	//!< file name of data record
	FILE*		m_fp;
};

//=========================================================================
// Super class for log data classes. 
class FECORE_API FELogData : public FECoreBase
{
public:
	FELogData(FEModel* fem);
};
