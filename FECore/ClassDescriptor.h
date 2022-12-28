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
#include <string>
#include <vector>
#include "fecore_api.h"

class FECORE_API FEClassDescriptor
{
public:
	class Variable
	{
	public:
		Variable(const std::string& name) : m_name(name) {}
		virtual ~Variable() {}

		virtual size_t Count() const = 0;

	public:
		std::string	m_name;
	};

	class SimpleVariable : public Variable
	{
	public:
		SimpleVariable(const std::string& name, const std::string& value) : Variable(name), m_val(value) {}

		size_t Count() const override { return 1; }

	public:
		std::string	m_val;
	};

	class ClassVariable : public Variable
	{
	public:
		ClassVariable(const std::string& name, const std::string& type) : Variable(name), m_type(type) {}
		size_t Count() const override { return m_var.size(); }

		void AddVariable(Variable* v) { m_var.push_back(v); }

		Variable* GetVariable(int i) { return m_var[i]; }
		const Variable* GetVariable(int i) const { return m_var[i]; }

	public:
		std::vector<Variable*> m_var;
		std::string	m_type;
	};

public:
	FEClassDescriptor(const std::string& classType) : m_var("root", classType) {}

	void AddVariable(Variable* v) { m_var.AddVariable(v); }

	ClassVariable* Root() { return &m_var; }
	const ClassVariable* Root() const { return &m_var; }

	SimpleVariable* FindParameter(const char* szparam);

	const std::string& ClassType() const { return m_var.m_type; }

public:
	ClassVariable	m_var;
};
