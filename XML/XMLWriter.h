/*This file is part of the FEBio Studio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio-Studio.txt for details.

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
#include <stdio.h>
#include <string.h>
#include <vector>
#include <string>
#include <assert.h>
class XMLWriter;

class XMLElement
{
public:
	class XMLAtt
	{
	public:
		XMLAtt() {}
		XMLAtt(const std::string& name, const std::string& val) : m_name(name), m_val(val) {}
		XMLAtt(const XMLAtt& att) : m_name(att.m_name), m_val(att.m_val) {}
		void operator = (const XMLAtt& att) { m_name = att.m_name; m_val = att.m_val; }

		const char* name() const { return m_name.c_str(); }
		const char* value() const { return m_val.c_str(); }

	public:
		std::string m_name;
		std::string m_val;
	};

public:
	XMLElement(const char* szname = 0)
	{
		clear();
		if (szname) m_tag = szname;
	}

	void clear()
	{
		m_att.clear();
		m_tag.clear();
		m_val.clear();
	}

	void name(const char* sz) { if (sz) m_tag = sz; else m_tag.clear(); }
	const char* name() const { return m_tag.c_str(); }

	void value(const char* sz) { if (sz) m_val = sz; else m_val.clear(); }
	void value(int    n);
	void value(int* pi, int n);
	void value(bool   b);
	void value(double g);
	void value(double* pg, int n);
	void value(const std::vector<int>& v);
	void value(const std::vector<double>& v);

	template <class T> void value(const T& v);

	int add_attribute(const char* szn, const char* szv);
	int add_attribute(const char* szn, int n);
	int add_attribute(const char* szn, bool b);
	int add_attribute(const char* szn, double g);
	int add_attribute(const char* szn, const std::string& s);

	void set_attribute(int nid, const char* szv);
	void set_attribute(int nid, int n);
	void set_attribute(int nid, bool b);
	void set_attribute(int nid, double g);

	int attributes() const { return (int) m_att.size(); }
	const XMLAtt& attribute(int i) const { return m_att[i]; }

protected:
	std::string	m_tag;		// element name
	std::string	m_val;		// element value

	std::vector<XMLAtt> m_att; // attributes

public:
	static void setDefaultFormats();
	static const char* intFormat;

	friend class XMLWriter;
};

class XMLWriter  
{
	enum {MAX_TAGS = 32};

public:
	enum XMLFloatFormat {
		ScientificFormat,
		FixedFormat
	};

public:
	XMLWriter();
	virtual ~XMLWriter();
	
	bool open(const char* szfile);
    bool setStringstream(std::ostringstream* stream);
    void init();

	void close();

	void add_branch(XMLElement& el, bool bclear = true);
	void add_branch(const char* szname);

	void add_empty(XMLElement& el, bool bclear = true);

	void add_leaf  (XMLElement& el, bool bclear = true);

	void add_leaf(const char* szn, const char* szv);
	void add_leaf(const char* szn, const std::string& s);

	void add_leaf(const char* szn, int    n){ char szv[256]; sprintf(szv, "%d" , n); write_leaf(szn, szv); }
	void add_leaf(const char* szn, bool   b){ char szv[256]; sprintf(szv, "%d" , b); write_leaf(szn, szv); }
	void add_leaf(const char* szn, double g){ char szv[256]; sprintf(szv, "%lg", g); write_leaf(szn, szv); }
	void add_leaf(const char* szn, float  g){ char szv[256]; sprintf(szv, "%g" , g); write_leaf(szn, szv); }
	void add_leaf(const char* szn, int *pi, int n);
	void add_leaf(const char* szn, float* pg, int n);
	void add_leaf(const char* szn, double* pg, int n);
	void add_leaf(XMLElement& el, const std::vector<int>& A);

	template <class T> void add_leaf(const char* szname, const T& v);

	void close_branch();

	void add_comment(const std::string& s, bool singleLine = false);

public:
	static void SetFloatFormat(XMLFloatFormat fmt);
	static XMLFloatFormat GetFloatFormat();

protected:
	void inc_level();
	void dec_level();

	void write_leaf(const char* sztag, const char* szval);

protected:
    std::ostream* m_stream;
	int		m_level;

	char	m_tag[MAX_TAGS][256];
	char	m_sztab[256];

	static XMLFloatFormat	m_floatFormat;
};

template <class T> std::string type_to_string(const T& v)
{
	return std::string(v);
}

template <class T> void XMLElement::value(const T& v)
{
	std::string s = type_to_string<T>(v);
	m_val = s;
}

template <class T> void XMLWriter::add_leaf(const char* szname, const T& v)
{
	std::string s = type_to_string<T>(v);
	add_leaf(szname, s);
}
