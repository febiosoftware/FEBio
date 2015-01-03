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
using namespace std;

//-------------------------------------------------------------------------
// forward declaration
class XMLReader;

//-------------------------------------------------------------------------
//! This class represents a xml-attribute
class XMLAtt
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
class XMLTag
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
class XMLReader  
{
public:
	enum {MAX_TAG   = 128};

public:
	// exceptions -----------
	class Error{};

	class EndOfFile{};

	class UnexpectedEOF{};

	class XMLSyntaxError{};

	class UnmatchedEndTag
	{
	public:
		XMLTag tag;
		UnmatchedEndTag(XMLTag& t) : tag(t) {}
	};

	class InvalidTag
	{
	public:
		XMLTag tag;
		InvalidTag(XMLTag& t) : tag(t) {}
	};

	class InvalidValue
	{
	public:
		XMLTag tag;
		InvalidValue(XMLTag& t) : tag(t) {}
	};

	class InvalidAttributeValue
	{
	public:
		XMLTag tag;
		char szatt[MAX_TAG];
		char szval[MAX_TAG];
		InvalidAttributeValue(XMLTag& t, const char* sza, const char* szv = 0);
	};

	class MissingAttribute
	{
	public:
		XMLTag tag;
		char szatt[MAX_TAG];
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
	bool FindTag(const char* sztag, XMLTag& tag);

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

protected:
	FILE*	m_fp;		//!< the file pointer
	int		m_nline;	//!< current line (used only as temp storage)
};

//-----------------------------------------------------------------------------
// some inline functions
inline void XMLTag::operator ++ () { m_preader->NextTag(*this); }

#endif // !defined(AFX_XMLREADER_H__667494D8_4C95_4342_BC31_6C1C097A4C81__INCLUDED_)
