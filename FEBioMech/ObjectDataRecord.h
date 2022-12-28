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
#include <FECore/FECoreBase.h>
#include <FECore/DataRecord.h>
#include "FERigidBody.h"

//-----------------------------------------------------------------------------
//! Base class for object log data (e.g. rigid bodies)
class FEBIOMECH_API FELogObjectData : public FELogData
{
	FECORE_SUPER_CLASS(FELOGOBJECTDATA_ID)
	FECORE_BASE_CLASS(FELogObjectData)

public:
	FELogObjectData(FEModel* fem) : FELogData(fem) {}
	virtual ~FELogObjectData(){}
	virtual double value(FERigidBody& rb) = 0;
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API ObjectDataRecord : public DataRecord
{
public:
	ObjectDataRecord(FEModel* pfem);
	double Evaluate(int item, int ndata) override;
	void SetData(const char* sz) override;
	void SelectAllItems() override;
	int Size() const override;

private:
	vector<FELogObjectData*>	m_Data;
};
