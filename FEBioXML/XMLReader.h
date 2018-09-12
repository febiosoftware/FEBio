// XMLReader.h: interface for the XMLReader class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_XMLREADER_H__667494D8_4C95_4342_BC31_6C1C097A4C81__INCLUDED_)
#define AFX_XMLREADER_H__667494D8_4C95_4342_BC31_6C1C097A4C81__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <FEBioLib/febiolib_api.h>
using namespace std;

//-------------------------------------------------------------------------
// forward declaration
class XMLReader;

//-------------------------------------------------------------------------
//! This class represents a xml-attribute
class FEBIOLIB_API XMLAtt
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
};

//-------------------------------------------------------------------------
//! This class implements a xml-tag. The value and attributes of this tag
//! can be queried.
//! \todo I would like to get rid of the m_szroot element and replace it with a 
//!       parent tag. The root element can then be identified by the tag that 
//!       does not have a parent
class FEBIOLIB_API XMLTag
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

	XMLReader*	m_preader;		// pointer to reader
	fpos_t	m_fpos;				// file position of next tag
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
	void value(bool& val);
	void value(char* szstr);
	void value(vector<int>& l);

	const char* szvalue() { return m_szval.c_str(); }
};

//-----------------------------------------------------------------------------
//! This class implements a reader for XML files
class FEBIOLIB_API XMLReader
{
public:
	enum {MAX_TAG   = 128};

	enum {BUF_SIZE = 1024};

public:
	// Base class for Exceptions
	class FEBIOLIB_API Error
	{
	public:
		enum { MAX_ERROR_STRING = 128 };

	public:
		Error() { m_szerr[0] = 0; }
		virtual ~Error(){}

		// retrieve the error string
		const char* GetErrorString() const { return m_szerr; }

	protected:
		// derived classes use this function to set the error string
		void SetErrorString(const char* sz, ...);

	protected:
		char	m_szerr[MAX_ERROR_STRING];
	};

	// End of file was discovered 
	class EndOfFile : public Error {};

	// the end of file was detected unexpectedly.
	class UnexpectedEOF : public Error {};

	// A syntax error was found
	class XMLSyntaxError : public Error
	{
	public:
		XMLSyntaxError();
	};

	// an end tag was not matched
	class UnmatchedEndTag : public Error
	{
	public:
		UnmatchedEndTag(XMLTag& t);
	};

	// an unknown tag was encountered 
	class FEBIOLIB_API InvalidTag : public Error
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
	fpos_t currentPos();

	//! move the file pointer
	void rewind(fpos_t nstep);

protected:
	FILE*	m_fp;		//!< the file pointer
	int		m_nline;	//!< current line (used only as temp storage)
	fpos_t	m_currentPos;	//!< current file position

	char	m_buf[BUF_SIZE];
#ifdef WIN32
    __int64    m_bufIndex, m_bufSize;
#else
    __int64_t    m_bufIndex, m_bufSize;
#endif
	bool	m_eof;
};

//-----------------------------------------------------------------------------
// some inline functions
inline void XMLTag::operator ++ () { m_preader->NextTag(*this); }

#endif // !defined(AFX_XMLREADER_H__667494D8_4C95_4342_BC31_6C1C097A4C81__INCLUDED_)
