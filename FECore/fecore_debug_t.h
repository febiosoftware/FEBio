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
#include <list>
#include <vector>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>
#include "matrix.h"
#include "mat3d.h"
#include "vec3d.h"
#include "tens4d.h"
#include "fecore_api.h"

// This file defines template constructions that are used by the FECore debugger.
// Don't use anything in here or include this file directly. 
// Instead include fecore_debug.h and use the user functions and macros defined there.

template <typename T> void fecore_print_T(T* pd) { std::cout << (*pd); }
template <> void fecore_print_T<matrix>(matrix* pd);
template <> void fecore_print_T<mat3d>(mat3d* pd);
template <> void fecore_print_T<mat3ds>(mat3ds* pd);
template <> void fecore_print_T<mat3da>(mat3da* pd);
template <> void fecore_print_T<mat3dd>(mat3dd* pd);
template <> void fecore_print_T<vec3d>(vec3d* pd);
template <> void fecore_print_T<tens4ds>(tens4ds* pd);
template <> void fecore_print_T<std::vector<double> >(std::vector<double>* pv);
template <> void fecore_print_T<std::vector<int> >(std::vector<int>* pv);

class FECoreBreakPoint;

class FECORE_API FECoreDebugger
{
public:
	class Variable
	{
	public:
		Variable(void* pd, const char* sz) { m_pd = pd; m_szname = sz; }
		virtual ~Variable() {}

		virtual void print() = 0;

	public:
		void* m_pd;
		const char* m_szname;
	};

	template <typename T> class Variable_T : public Variable
	{
	public:
		Variable_T(void* pd, const char* sz) : Variable(pd, sz) {}

		void print()
		{
			std::cout << typeid(T).name() << std::endl;
			fecore_print_T<T>((T*)m_pd);
		}
	};

public:
	static void Break(FECoreBreakPoint* pbr);
	static void Clear();
	static void Add(Variable* pvar);
	static void Remove(Variable* pvar);

	static void Print(const char* szformat, ...);

private:
	FECoreDebugger() {}

	static std::list<Variable*>	m_var;
	static FILE* m_fp;
};

class FECoreWatchVariable
{
public:
	FECoreWatchVariable(FECoreDebugger::Variable* pvar) : m_pvar(pvar)
	{
		if (pvar) FECoreDebugger::Add(pvar);
	}

	~FECoreWatchVariable()
	{
		if (m_pvar) FECoreDebugger::Remove(m_pvar);
	}

protected:
	FECoreDebugger::Variable* m_pvar;
};

class FECoreBreakPoint
{
public:
	FECoreBreakPoint() { static int n = 1; m_bactive = true; m_nid = n++; }
	void Deactivate() { m_bactive = false; }
	bool IsActive() { return m_bactive; }

	void Break()
	{
		if (m_bactive) { FECoreDebugger::Break(this); }
	}

	int GetID() { return m_nid; }

private:
	int		m_nid;
	bool	m_bactive;
};

template <typename T> FECoreDebugger::Variable* create_watch_variable(T* pd, const char* sz) { return new FECoreDebugger::Variable_T<T>(pd, sz); }
