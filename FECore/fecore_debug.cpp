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



#include "stdafx.h"
#include "fecore_debug.h"
#include <string>
#include <iostream>
#include <stdarg.h>
#include "version.h"
using namespace std;

std::list<FECoreDebugger::Variable*>	FECoreDebugger::m_var;

void FECoreDebugger::Break(FECoreBreakPoint* pbr)
{
	string s;
	do
	{
		if (pbr) cout << "#" << pbr->GetID();
		cout << ">>";
		cin >> s;

		if      (s == "cont" ) break;
		else if (s == "print")
		{
			cin >> s;
			list<Variable*>::iterator it;
			bool bfound = false;
			for (it = m_var.begin(); it != m_var.end(); ++it)
			{
				if (s == (*it)->m_szname)
				{
					(*it)->print();
					cout << endl;
					bfound = true;
					break;
				}
			}
			if (bfound == false) cout << "Error: unknown variable\n";
		}
		else if (s == "list")
		{
			list<Variable*>::iterator it;
			for (it = m_var.begin(); it != m_var.end(); ++it)
			{
				cout << (*it)->m_szname << endl;
			}
		}
		else if (s == "remove")
		{
			if (pbr) { pbr->Deactivate(); break; }
			else cout << "No active breakpoint\n";
		}
		else if (s == "help")
		{
			cout << "cont      = continue\n";
			cout << "list      = list all watch variables\n";
			cout << "print var = print variable\n";
			cout << "remove    = remove breakpoint and continue\n";
			cout << "help      = show this information\n";
		}
		else cout << "Error: Unknown command\n";
	}
	while (1);
}

void FECoreDebugger::Clear()
{
	list<Variable*>::iterator it;
	for (it = m_var.begin(); it != m_var.end(); ++it) delete (*it);
	m_var.clear();
}

void FECoreDebugger::Add(FECoreDebugger::Variable* pvar)
{
	m_var.push_back(pvar);
}

void FECoreDebugger::Remove(FECoreDebugger::Variable* pvar)
{
	list<Variable*>::iterator it;
	for (it = m_var.begin(); it != m_var.end(); ++it)
	{
		if ((*it)==pvar)
		{
			delete *it;
			m_var.erase(it);
			break;
		}
	}
}

FILE* FECoreDebugger::m_fp = nullptr;

void FECoreDebugger::Print(const char* szformat, ...)
{
	if (m_fp == nullptr)
	{
		char fileName[256] = { 0 };
		sprintf(fileName, "febio_%d.%d.%d_debug.log", FE_SDK_MAJOR_VERSION, FE_SDK_SUB_VERSION, FE_SDK_SUBSUB_VERSION);
		m_fp = fopen(fileName, "wt"); assert(m_fp);
	}

	if (m_fp)
	{
		va_list	args;
		va_start(args, szformat);
		vfprintf(m_fp, szformat, args);
		va_end(args);

		fflush(m_fp);
	}
}

template <> void fecore_print_T<matrix>(matrix* pd)
{
	matrix& m = *pd;
	int nr = m.rows();
	int nc = m.columns();
	for (int i=0; i<nr; ++i)
	{
		for (int j=0; j<nc; ++j)
		{
			cout << m(i,j) << " ";
		}
		cout << endl;
	}
}

template <> void fecore_print_T<mat3d>(mat3d* pd)
{
	mat3d& m = *pd;
	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j) cout << m(i,j) << " ";
		cout << endl;
	}
}

template <> void fecore_print_T<mat3ds>(mat3ds* pd)
{
	mat3ds& m = *pd;
	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j) cout << m(i,j) << " ";
		cout << endl;
	}
}

template <> void fecore_print_T<mat3da>(mat3da* pd)
{
	mat3da& m = *pd;
	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j) cout << m(i,j) << " ";
		cout << endl;
	}
}

template <> void fecore_print_T<mat3dd>(mat3dd* pd)
{
	mat3dd& m = *pd;
	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j) cout << m(i,j) << " ";
		cout << endl;
	}
}

template <> void fecore_print_T<vec3d>(vec3d* pd)
{
	vec3d& v = *pd;
	cout << v.x << endl;
	cout << v.y << endl;
	cout << v.z << endl;
}

template <> void fecore_print_T<tens4ds>(tens4ds* pd)
{
	tens4ds& m = *pd;
	for (int i=0; i<6; ++i)
	{
		for (int j=0; j<6; ++j)
		{
			cout << m(i,j) <<" ";
		}
		cout << endl;
	}
}

template <> void fecore_print_T<std::vector<double> >(std::vector<double>* pv)
{
	vector<double>& v = *pv;
	int n = (int)v.size();
	for (int i=0; i<n; ++i) cout << v[i] << endl;
}

template <> void fecore_print_T<std::vector<int> >(std::vector<int>* pv)
{
	vector<int>& v = *pv;
	int n = (int)v.size();
	for (int i=0; i<n; ++i) cout << v[i] << endl;
}
