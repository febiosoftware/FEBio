/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FECoreBase.h"
#include "DataRecord.h"

class FENodeSet;

//-----------------------------------------------------------------------------
//! This is the base class for a node data value.
//! \todo I'd like to modify this so I can pass the FENode class instead of the node number
class FECORE_API FENodeLogData : public FECoreBase
{ 
	FECORE_SUPER_CLASS

public:
	FENodeLogData(FEModel* fem);
	virtual ~FENodeLogData();
	virtual double value(int node) = 0; 
};

//-----------------------------------------------------------------------------
//! This class records nodal data
//! \todo should I create a different class for each data record? Like for the plot file?
class FECORE_API NodeDataRecord : public DataRecord
{
public:
	NodeDataRecord(FEModel* pfem, const char* szfile);
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	void SetItemList(FENodeSet* pns);
	int Size() const;

private:
	vector<FENodeLogData*>	m_Data;
};

//-----------------------------------------------------------------------------
// Special class for outputting nodal variables
class FECORE_API FENodeVarData : public FENodeLogData
{
public:
	FENodeVarData(FEModel* pfem, int ndof);
	double value(int node);
private:
	int	m_ndof;
};
