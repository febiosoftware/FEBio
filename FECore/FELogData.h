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
#include "fecore_api.h"

// Super class for log data classes. 
class FECORE_API FELogData : public FECoreBase
{
public:
	FELogData(FEModel* fem);

	virtual bool SetParameters(const std::vector<std::string>& params) { return params.empty(); }

	virtual bool SetComponent(const std::string& comp) { return false; }

	virtual bool SetIndex(int n) { return false; }
};

template <class Base, class Item, class T, class Traits>
class FELogComponentData : public Base
{
public:
	using Base::Base;

	bool SetComponent(const std::string& comp) override
	{
		m_comp = Traits::ComponentIndex(comp);
		return m_comp >= 0;
	}

	bool Init() override
	{
		return (m_comp >= 0) && Base::Init();
	}

	double value(Item& item) final
	{
		T v = typedValue(item);
		return Traits::Component(v, m_comp);
	}

protected:
	virtual T typedValue(Item& item) = 0;

private:
	int m_comp = -1;
};

struct FECORE_API Vec3dLogTraits
{
	static int ComponentIndex(const std::string& c);
	static double Component(const vec3d& v, int i);
};

struct FECORE_API Mat3dsLogTraits
{
	static int ComponentIndex(const std::string& c);
	static double Component(const mat3ds& v, int i);
};

struct FECORE_API Mat3dLogTraits
{
	static int ComponentIndex(const std::string& c);
	static double Component(const mat3d& v, int i);
};
