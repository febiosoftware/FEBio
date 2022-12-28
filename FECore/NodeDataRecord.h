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
#include "FECoreBase.h"
#include "DataRecord.h"

class FENodeSet;
class FENode;

//-----------------------------------------------------------------------------
//! This is the base class for a node data value.
class FECORE_API FELogNodeData : public FELogData
{ 
	FECORE_SUPER_CLASS(FELOGNODEDATA_ID)
	FECORE_BASE_CLASS(FELogNodeData)

public:
	FELogNodeData(FEModel* fem);
	virtual ~FELogNodeData();
	virtual double value(const FENode& node) = 0; 
};

//-----------------------------------------------------------------------------
//! This class records nodal data
//! \todo should I create a different class for each data record? Like for the plot file?
class FECORE_API NodeDataRecord : public DataRecord
{
public:
	NodeDataRecord(FEModel* pfem);
	double Evaluate(int item, int ndata);
	void SetData(const char* sz) override;
	void SelectAllItems();
	int Size() const;

	void SetItemList(FEItemList* items, const std::vector<int>& selection) override;

private:
	vector<FELogNodeData*>	m_Data;
};

//-----------------------------------------------------------------------------
// Special class for outputting nodal variables
class FECORE_API FENodeVarData : public FELogNodeData
{
public:
	FENodeVarData(FEModel* pfem, int ndof);
	double value(const FENode& node) override;

private:
	int	m_ndof;
};
