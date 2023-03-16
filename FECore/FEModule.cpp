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
#include "FEModule.h"
#include <vector>
#include <cstring>

class FEModule::Impl
{
public:
	const char*			szname = 0;		// name of module
	const char*			szdesc = 0;		// description of module (optional, can be null)
	unsigned int		id = 0;			// unqiue ID (starting at one)
	int					alloc_id = -1;	// ID of allocator
	int					m_status = FEModule::RELEASED;	// Status of module
	std::vector<int>	depMods;	// module dependencies

public:
	void AddDependency(int mid)
	{
		for (size_t i = 0; i < depMods.size(); ++i)
		{
			if (depMods[i] == mid) return;
		}
		depMods.push_back(mid);
	}

	void AddDependencies(const std::vector<int>& mid)
	{
		for (size_t i = 0; i < mid.size(); ++i)
		{
			AddDependency(mid[i]);
		}
	}
};

FEModule::FEModule() : im(new FEModule::Impl)
{

}

FEModule::FEModule(const char* szname, const char* szdescription) : im(new FEModule::Impl)
{
	SetName(szname);
	SetDescription(szdescription);
}

FEModule::~FEModule()
{
	delete im;
}

void FEModule::InitModel(FEModel* fem)
{

}

int FEModule::GetModuleID() const
{
	return im->id;
}

const char* FEModule::GetName() const
{
	return im->szname;
}

const char* FEModule::GetDescription() const
{
	return im->szdesc;
}

void FEModule::SetStatus(FEModule::Status status)
{
	im->m_status = status;
}

int FEModule::GetStatus() const
{
	return im->m_status;
}

bool FEModule::HasDependent(int modId) const
{
	if (modId == im->id) return true;
	for (int i : im->depMods)
	{
		if (i == modId) return true;
	}
	return false;
}

void FEModule::AddDependency(FEModule& mod)
{
	im->AddDependency(mod.GetModuleID());
	im->AddDependencies(mod.im->depMods);
}

void FEModule::ClearDependencies()
{
	im->depMods.clear();
}

std::vector<int> FEModule::GetDependencies() const
{
	return im->depMods;
}

void FEModule::SetID(int newId)
{
	im->id = newId;
}

void FEModule::SetName(const char* szname)
{
	im->szname = szname;
}

void FEModule::SetDescription(const char* szdesc)
{
	im->szdesc = szdesc;
}
