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
#include <assert.h>

//-------------------------------------------------------------------------
// forward declaration
class XMLReader;

//-------------------------------------------------------------------------
//! This class represents a xml-attribute
class XMLAtt
{
public:
	//! constructor
	XMLAtt();

	//! comparison operators
	bool operator == (const char* sz);
	bool operator != (const char* szval);

	//! Get the attribute name
	const char* name() { return m_name.c_str(); }

	//! Get the attribute value
	const char* cvalue() { return m_val.c_str(); }

public:
	void value(int& val) { val = atoi(m_val.c_str()); }
	int value(double* v, int n);
	template <typename T> T value() { return T(0); }

public:
	std::string	m_name;
	std::string m_val;
	bool	m_bvisited;			//!< was the attribute processed or not?
};

template <> inline int XMLAtt::value<int>() { return atoi(m_val.c_str()); }
template <> inline double XMLAtt::value<double>() { return atof(m_val.c_str()); }
template <> inline std::string XMLAtt::value<std::string>() { return m_val; }

//-------------------------------------------------------------------------
//! This class implements a xml-tag. The value and attributes of this tag
//! can be queried.
//! \todo I would like to get rid of the m_szroot element and replace it with a 
//!       parent tag. The root element can then be identified by the tag that 
//!       does not have a parent
class XMLTag
{
public:
	std::string	m_sztag;			// tag name
	std::string m_szval;				// tag value

	std::vector<XMLAtt>	m_att;				// attribute list

	std::vector<std::string> m_path;	// current path (name tags of parents)

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

	int currentLine() const { return m_ncurrent_line; }

	bool operator == (const char* sztag) { return (strcmp(sztag, m_sztag.c_str()) == 0); }
	bool operator != (const char* sztag) { return (strcmp(sztag, m_sztag.c_str()) != 0); }
	void operator ++ ();

	void skip();

	bool isend() { return m_bend; }
	bool isleaf() { return m_bleaf; }
	bool isempty() { return m_bempty; }

	// count the number of children
	int children();

	const char* AttributeValue(const char* szat, bool bopt = false);
	XMLAtt* AttributePtr(const char* szat);
	XMLAtt* Attribute(const char* szat, bool bopt);
	XMLAtt& Attribute(const char* szat);

	template <typename T> T AttributeValue(const char* szatt, const T& def_val) { return def_val; }

	bool AttributeValue(const char* szat, int&    n, bool bopt = false);
	bool AttributeValue(const char* szat, double& d, bool bopt = false);
		
	const char* Name() { return m_sztag.c_str(); }

	void value(double& val) { val = atof(m_szval.c_str()); } 
	void value(float& val)  { val = (float) atof(m_szval.c_str()); }
	void value(int& val) { val = atoi(m_szval.c_str()); }
	void value(long& val) { val = (long) atoi(m_szval.c_str()); }
	void value(short& val) { val = (short) atoi(m_szval.c_str()); }
	int value(double* pf, int n);
	int value(float* pf, int n);
	int value(int* pi, int n);
	int value(std::vector<std::string>& stringList, int n);
	void value(bool& val);
	void value(char* szstr);
	void value(std::string& val);
	void value(std::vector<int>& l);
	void value2(std::vector<int>& l);
	void value(std::vector<double>& l);

	template <class T> void value(T& v);

	const char* szvalue() { return m_szval.c_str(); }

	std::string relpath(const char* szroot) const;

	const std::string& comment();
    const void clearComment();
};


template <> inline int XMLTag::AttributeValue<int>(const char* szatt, const int& def_val)
{
	XMLAtt* pa = AttributePtr(szatt);
	if (pa) return pa->value<int>();
	else return def_val;
}

template <> inline double XMLTag::AttributeValue<double >(const char* szatt, const double& def_val)
{
	XMLAtt* pa = AttributePtr(szatt);
	if (pa) return pa->value<double>();
	else return def_val;
}

template <> inline std::string XMLTag::AttributeValue<std::string>(const char* szatt, const std::string& def_val)
{
	XMLAtt* pa = AttributePtr(szatt);
	if (pa) return pa->value<std::string>();
	else return def_val;
}

//-----------------------------------------------------------------------------
//! This class implements a reader for XML files
class XMLReader
{
public:
	enum {BUF_SIZE = 32768};

public:
	// Base class for Exceptions
	class Error : public std::runtime_error
	{
	public:
		Error(const std::string& err) : std::runtime_error(err) {}
		Error(XMLTag& tag, const std::string& err);
	};

	// End of file was discovered 
	class EndOfFile : public Error {
	public: 
		EndOfFile() : Error("End of file") {}
	};

	// the end of file was detected unexpectedly.
	class UnexpectedEOF : public Error {
	public:
		UnexpectedEOF() : Error("Unexpected end of file") {}
	};

	// A syntax error was found
	class XMLSyntaxError : public Error
	{
	public:
		XMLSyntaxError(int line_number = -1);
	};

	// an end tag was not matched
	class UnmatchedEndTag : public Error
	{
	public:
		UnmatchedEndTag(XMLTag& t);
	};

	// an unknown tag was encountered 
	class InvalidTag : public Error
	{
	public:
		InvalidTag(XMLTag& t);
	};

	// The value of a tag was invald 
	class InvalidValue : public Error
	{
	public:
		InvalidValue(XMLTag& t);
	};

	// the value of an attribute was invalid 
	class InvalidAttributeValue : public Error
	{
	public:
		InvalidAttributeValue(XMLTag& t, const char* sza, const char* szv = 0);
		InvalidAttributeValue(XMLTag& t, XMLAtt& att);
	};

	// an attribute is invalid
	class InvalidAttribute : public Error
	{
	public:
		InvalidAttribute(XMLTag& t, const char* sza);
	};

	// an attribute was missing
	class MissingAttribute : public Error
	{
	public:
		MissingAttribute(XMLTag& t, const char* sza);
	};

	class MissingTag : public Error
	{
	public:
		MissingTag(XMLTag& t, const char* szt);
	};
	//------------------------

public:
	//! constructor/destructor
	XMLReader();
	virtual ~XMLReader();

    std::ifstream* GetFileStream();

	//! Open the xml file
	bool Open(const char* szfile, bool checkForXMLTag = true);

    //! Pass xml formatted string to reader
    bool OpenString(std::string& xml, bool checkForXMLTag = true);

	//! Close the xml file
	void Close();

	//! Find a tag
	bool FindTag(const char* xpath, XMLTag& tag);

	//! Get the next tag
	void NextTag(XMLTag& tag);

	//! return the current line
	int GetCurrentLine();

	//! Skip a tag
	void SkipTag(XMLTag& tag);

	const std::string& GetLastComment();

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

	// only used for processing comments
	char GetNextChar();
	
protected:
    std::istream* m_stream;

	int		m_nline;		//!< current line (used only as temp storage)
    int64_t	m_currentPos;	//!< current file position

	std::string	m_comment;	//!< last comment that was read

	char*		m_buf;
    int64_t    m_bufIndex, m_bufSize;
	bool		m_eof;
};

//-----------------------------------------------------------------------------
// some inline functions
inline void XMLTag::operator ++ () { m_preader->NextTag(*this); }

inline void XMLTag::skip() { m_preader->SkipTag(*this); }

inline const std::string& XMLTag::comment() { return m_preader->GetLastComment(); }

//-----------------------------------------------------------------------------
// mechanism for using custom types with XMLReader. 
template <class T> void string_to_type(const std::string& s, T& v) { assert(false); }
template <class T> void XMLTag::value(T& v) { string_to_type(m_szval, v); }
