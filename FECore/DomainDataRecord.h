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
#include "ElementDataRecord.h"

class FEDomain;

//-----------------------------------------------------------------------------
//! Base class for domain log data
class FECORE_API FELogDomainData : public FELogData
{
    FECORE_SUPER_CLASS(FELOGDOMAINDATA_ID)
    FECORE_BASE_CLASS(FELogDomainData)

public:
    FELogDomainData(FEModel* fem) : FELogData(fem) {}
    virtual ~FELogDomainData() {}
    virtual double value(FEDomain& rc) = 0;

    virtual bool SetParameters(std::vector<std::string>& params) { return false; }
};

//-----------------------------------------------------------------------------
class FECORE_API FEDomainDataRecord : public DataRecord
{
public:
    FEDomainDataRecord(FEModel* pfem);
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SelectAllItems();
    void SetDomain(int domainIndex);
    int Size() const;

private:
    vector<FELogDomainData*>	m_Data;
};

//-----------------------------------------------------------------------------
class FECORE_API FELogAvgDomainData : public FELogDomainData
{
public:
    FELogAvgDomainData(FEModel* pfem);
    ~FELogAvgDomainData();
    double value(FEDomain& rc) override;

    bool SetParameters(std::vector<std::string>& params);

private:
    FELogElemData* m_elemData;
};

//-----------------------------------------------------------------------------
class FECORE_API FELogPctDomainData : public FELogDomainData
{
public:
    FELogPctDomainData(FEModel* pfem);
    ~FELogPctDomainData();
    double value(FEDomain& rc) override;

    bool SetParameters(std::vector<std::string>& params);

private:
    double          m_pct;
    FELogElemData* m_elemData;
};

//-----------------------------------------------------------------------------
class FECORE_API FELogIntegralDomainData : public FELogDomainData
{
public:
	FELogIntegralDomainData(FEModel* pfem);
	~FELogIntegralDomainData();
	double value(FEDomain& rc) override;

	bool SetParameters(std::vector<std::string>& params);

private:
	FELogElemData* m_elemData;
};
