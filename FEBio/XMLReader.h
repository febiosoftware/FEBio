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
#include "FECore/vec3d.h"

class XMLReader  
{
public:
	enum {MAX_TAG   = 128};
	enum {MAX_ATT   =   8};
	enum {MAX_LEVEL = 16};

	class XMLTag
	{
	// TODO: I would like to get rid of the m_szroot element and replace it with a 
	// parent tag. The root element can then be identified by the tag that 
	// does not have a parent
	public:
		char	m_sztag[MAX_TAG];			// tag name
		char	m_szval[MAX_TAG];			// tag value
		char	m_szatt[MAX_ATT][MAX_TAG];	// tag attributes
		char	m_szatv[MAX_ATT][MAX_TAG];	// value of attributes

		int		m_nlevel;	// depth level
		char	m_szroot[MAX_LEVEL][MAX_TAG];	// name tag of parent's

		int		m_natt;	// nr of attributes

		XMLReader*	m_preader;		// pointer to reader
		fpos_t	m_fpos;				// file position of next tag
		int		m_nstart_line;		// line number at beginning of tag
		int		m_ncurrent_line;	// current line number

		bool	m_bend;		// end tag flag
		bool	m_bleaf;	// this is a leaf (i.e. has no child elements)
		bool	m_bempty;	// empty tag (i.e. no value)

		void operator ++ () { m_preader->NextTag(*this); }


		XMLTag();

		void clear()
		{
			m_sztag[0] = 0;
			m_szval[0] = 0;
			m_natt = 0;
			m_bend = false;
			m_bleaf = true;
			m_bempty = false;
		}

	public:
		bool operator == (const char* sztag) { return (strcmp(sztag, m_sztag) == 0); }
		bool operator != (const char* sztag) { return (strcmp(sztag, m_sztag) != 0); }

		bool isend() { return m_bend; }
		bool isleaf() { return m_bleaf; }
		bool isempty() { return m_bempty; }

	public:
		const char* AttributeValue(const char* szat, bool bopt = false);

		bool AttributeValue(const char* szat, int&    n, bool bopt = false);
		bool AttributeValue(const char* szat, double& d, bool bopt = false);
		
		const char* Name() { return m_sztag; }

		void value(char* szstr) { strcpy(szstr, m_szval); }
		void value(double& val) { val = atof(m_szval); } 
		void value(float& val)  { val = (float) atof(m_szval); }
		void value(bool& val) { int n=0; sscanf(m_szval, "%d", &n); val = (n != 0); }
		void value(int& val) { val = atoi(m_szval); }
		void value(double* pf, int n);
		void value(float* pf, int n);
		void value(int* pi, int n);
		void value(vec3d& v);

		const char* szvalue() { return m_szval; }
	};

	// exceptions -----------

	class Error{};

	class UnexpectedEOF{};

	class XMLSyntaxError{};

	class UnmatchedEndTag
	{
	public:
		XMLTag tag;
		UnmatchedEndTag(XMLTag& t) : tag(t) {}
	};

	class EndOfBuffer
	{
	public:
		XMLTag tag;
		EndOfBuffer(XMLTag& t) : tag(t) {}
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
		InvalidAttributeValue(XMLTag& t, const char* sza, const char* szv) : tag(t)
		{ 
			strcpy(szatt, sza);
			strcpy(szval, szv);
		}
	};

	class MissingAttribute
	{
	public:
		XMLTag tag;
		char szatt[MAX_TAG];
		MissingAttribute(XMLTag& t, const char* sza) : tag(t) { strcpy(szatt, sza); }
	};

	//------------------------


public:
	XMLReader();
	virtual ~XMLReader();

	FILE* GetFilePtr() { return m_fp; } 

	bool Open(const char* szfile);

	void Close();

	bool FindTag(const char* sztag, XMLTag& tag);

	void NextTag(XMLTag& tag);

	int GetCurrentLine() { return m_nline; }

protected:
	char GetChar() 
	{
		char ch;
		while ((ch=fgetc(m_fp))=='\n') ++m_nline;
		if (feof(m_fp)) throw UnexpectedEOF();
		return ch;
	}

	void ReadTag(XMLTag& tag);
	void ReadValue(XMLTag& tag);
	void ReadEndTag(XMLTag& tag);

protected:
	FILE*	m_fp;		// the file pointer
	int		m_nline;	// current line (used only as temp storage)

	friend class XMLTag;
};

typedef XMLReader::XMLTag XMLTag;

#endif // !defined(AFX_XMLREADER_H__667494D8_4C95_4342_BC31_6C1C097A4C81__INCLUDED_)
