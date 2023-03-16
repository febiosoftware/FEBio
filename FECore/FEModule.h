/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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

class FEModel;
class FECoreKernel;

// The FEModule class defines the physics variables for a particular module
class FECORE_API FEModule
{
	class Impl;

public:
	enum Status {
		EXPERIMENTAL,
		RELEASED
	};

public:
	FEModule();
	FEModule(const char* szname, const char* szdescription = nullptr);

	virtual ~FEModule();

	// this function must be overridden by derived classes
	virtual void InitModel(FEModel* fem);

	void AddDependency(FEModule& mod);

	void ClearDependencies();

	std::vector<int> GetDependencies() const;

	int GetModuleID() const;

	const char* GetName() const;

	const char* GetDescription() const;

	int GetStatus() const;

	bool HasDependent(int modId) const;

protected:
	void SetStatus(FEModule::Status status);

private:
	void SetID(int newId);
	void SetName(const char* szname);
	void SetDescription(const char* szdesc);

public:
	Impl* im;

	friend class FECoreKernel;
};
