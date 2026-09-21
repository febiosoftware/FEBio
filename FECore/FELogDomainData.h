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
#pragma once
#include "FELogData.h"

class FEDomain;
class FELogElemData;

//! Base class for domain log data
class FECORE_API FELogDomainData : public FELogData
{
	FECORE_SUPER_CLASS(FELOGDOMAINDATA_ID)
		FECORE_BASE_CLASS(FELogDomainData)

public:
	FELogDomainData(FEModel* fem) : FELogData(fem) {}
	virtual ~FELogDomainData() {}
	virtual double value(FEDomain& rc) = 0;
};

using FELogDomainVec3dData  = FELogComponentData<FELogDomainData, FEDomain, vec3d , Vec3dLogTraits>;
using FELogDomainMat3dsData = FELogComponentData<FELogDomainData, FEDomain, mat3ds, Mat3dsLogTraits>;
using FELogDomainMat3dData  = FELogComponentData<FELogDomainData, FEDomain, mat3d , Mat3dLogTraits>;

class FECORE_API FELogAvgDomainData : public FELogDomainData
{
public:
	FELogAvgDomainData(FEModel* pfem);
	~FELogAvgDomainData();
	double value(FEDomain& rc) override;

	bool SetParameters(const std::vector<std::string>& params) override;

private:
	FELogElemData* m_elemData;
};

class FECORE_API FELogPctDomainData : public FELogDomainData
{
public:
	FELogPctDomainData(FEModel* pfem);
	~FELogPctDomainData();
	double value(FEDomain& rc) override;

	bool SetParameters(const std::vector<std::string>& params) override;

private:
	double          m_pct;
	FELogElemData* m_elemData;
};

class FECORE_API FELogIntegralDomainData : public FELogDomainData
{
public:
	FELogIntegralDomainData(FEModel* pfem);
	~FELogIntegralDomainData();
	double value(FEDomain& rc) override;

	bool SetParameters(const std::vector<std::string>& params) override;

private:
	FELogElemData* m_elemData;
};
