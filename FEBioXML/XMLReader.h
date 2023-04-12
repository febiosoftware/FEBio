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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "febioxml_api.h"
using namespace std;

//-------------------------------------------------------------------------
// forward declaration
class XMLReader;

//-------------------------------------------------------------------------
//! This class represents a xml-attribute
class FEBIOXML_API XMLAtt
{
	//! max buffer size for attribute name and value
	enum { MAX_TAG = 128 };

public:
	//! constructor
	XMLAtt();

	//! assignment operator
	bool operator == (const char* sz);

	//! Get the attribute name
	const char* name() { return m_szatt; }

	//! Get the attribute value
	const char* cvalue() { return m_szatv; }

public:
	char	m_szatt[MAX_TAG];	//!< attribute name
	char	m_szatv[MAX_TAG];	//!< attribute value
	bool	m_bvisited;			//!< was the attribute processed or not?
};

//-------------------------------------------------------------------------
//! This class implements a xml-tag. The value and attributes of this tag
//! can be queried.
//! \todo I would like to get rid of the m_szroot element and replace it with a 
//!       parent tag. The root element can then be identified by the tag that 
//!       does not have a parent
class FEBIOXML_API XMLTag
{
public:
	enum {MAX_TAG   = 128};
	enum {MAX_ATT   =   8};
	enum {MAX_LEVEL =  16};

public:
	char		m_sztag[MAX_TAG];		// tag name
	std::string m_szval;				// tag value

	XMLAtt	m_att[MAX_ATT];				// attribute list
	int		m_natt;						// nr of attributes

	int		m_nlevel;						// depth level
	char	m_szroot[MAX_LEVEL][MAX_TAG];	// name tag of parent's

	XMLReader*	m_preader;			// pointer to reader
    int64_t		m_fpos;				// file position of next tag
	int		m_nstart_line;		// line number at beginning of tag
	int		m_ncurrent_line;	// current line number

	bool	m_bend;		// end tag flag
	bool	m_bleaf;	// this is a leaf (i.e. has no child elements)
	bool	m_bempty;	// empty tag (i.e. no value)

public:
	XMLTag();

	void clear();

	bool operator == (const char* sztag) { return (strcmp(sztag, m_sztag) == 0); }
	bool operator != (const char* sztag) { return (strcmp(sztag, m_sztag) != 0); }
	void operator ++ ();

	bool isend() { return m_bend; }
	bool isleaf() { return m_bleaf; }
	bool isempty() { return m_bempty; }

	// count the number of children
	int children();

	const char* AttributeValue(const char* szat, bool bopt = false);
	XMLAtt* Attribute(const char* szat, bool bopt);
	XMLAtt& Attribute(const char* szat);

	bool AttributeValue(const char* szat, int&    n, bool bopt = false);
	bool AttributeValue(const char* szat, double& d, bool bopt = false);
		
	const char* Name() { return m_sztag; }

	void value(double& val) { val = atof(m_szval.c_str()); } 
	void value(float& val)  { val = (float) atof(m_szval.c_str()); }
	void value(int& val) { val = atoi(m_szval.c_str()); }
	void value(long& val) { val = (long) atoi(m_szval.c_str()); }
	void value(short& val) { val = (short) atoi(m_szval.c_str()); }
	int value(double* pf, int n);
	int value(float* pf, int n);
	int value(int* pi, int n);
	int value(std::vector<string>& stringList, int n);
	void value(bool& val);
	void value(char* szstr);
	void value(std::string& val);
	void value(vector<int>& l);

	const char* szvalue() { return m_szval.c_str(); }
};

//-----------------------------------------------------------------------------
//! This class implements a reader for XML files
class FEBIOXML_API XMLReader
{
public:
	enum {MAX_TAG   = 128};

	enum {BUF_SIZE = 32768};

public:
	// Base class for Exceptions
	class FEBIOXML_API Error : public std::runtime_error
	{
	public:
		Error(const std::string& err) : std::runtime_error(err) {}
		Error(XMLTag& tag, const std::string& err);
	};

	// End of file was discovered 
	class FEBIOXML_API EndOfFile : public Error {
	public: 
		EndOfFile() : Error("End of file") {}
	};

	// the end of file was detected unexpectedly.
	class FEBIOXML_API UnexpectedEOF : public Error {
	public:
		UnexpectedEOF() : Error("Unexpected end of file") {}
	};

	// A syntax error was found
	class FEBIOXML_API XMLSyntaxError : public Error
	{
	public:
		XMLSyntaxError(int line_number = -1);
	};

	// an end tag was not matched
	class FEBIOXML_API UnmatchedEndTag : public Error
	{
	public:
		UnmatchedEndTag(XMLTag& t);
	};

	// an unknown tag was encountered 
	class FEBIOXML_API InvalidTag : public Error
	{
	public:
		InvalidTag(XMLTag& t);
	};

	// The value of a tag was invald 
	class FEBIOXML_API InvalidValue : public Error
	{
	public:
		InvalidValue(XMLTag& t);
	};

	// the value of an attribute was invalid 
	class FEBIOXML_API InvalidAttributeValue : public Error
	{
	public:
		InvalidAttributeValue(XMLTag& t, const char* sza, const char* szv = 0);
	};

	// an attribute is invalid
	class FEBIOXML_API InvalidAttribute : public Error
	{
	public:
		InvalidAttribute(XMLTag& t, const char* sza);
	};

	// an attribute was missing
	class FEBIOXML_API MissingAttribute : public Error
	{
	public:
		MissingAttribute(XMLTag& t, const char* sza);
	};
	//------------------------

public:
	//! constructor/destructor
	XMLReader();
	virtual ~XMLReader();

	//! Open the xml file
	bool Open(const char* szfile);

	//! Close the xml file
	void Close();

	//! Find a tag
	bool FindTag(const char* xpath, XMLTag& tag);

	//! Get the next tag
	void NextTag(XMLTag& tag);

	//! return the current line
	int GetCurrentLine() { return m_nline; }

	//! Skip a tag
	void SkipTag(XMLTag& tag);

protected: // helper functions

	//! Get the next character in the file
	char GetChar();

	//! Read a tag
	void ReadTag(XMLTag& tag);

	//! Read the value of a tag
	void ReadValue(XMLTag& tag);

	//! process end tag
	void ReadEndTag(XMLTag& tag);

	//! read the next character of the buffer
	char readNextChar();

	//! get the current position
    int64_t currentPos();

	//! move the file pointer
    void rewind(int64_t nstep);

protected:
	FILE*	m_fp;			//!< the file pointer
	int		m_nline;		//!< current line (used only as temp storage)
    int64_t	m_currentPos;	//!< current file position

	char*	m_buf;
    int64_t    m_bufIndex, m_bufSize;
	bool	m_eof;
};

//-----------------------------------------------------------------------------
// some inline functions
inline void XMLTag::operator ++ () { m_preader->NextTag(*this); }
