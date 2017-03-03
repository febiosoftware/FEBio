// XMLReader.cpp: implementation of the XMLReader class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "XMLReader.h"
#include <assert.h>
#include <stdarg.h>

//=============================================================================
// XMLAtt
//=============================================================================

//-----------------------------------------------------------------------------
//! The constructor sets the name and value to zero strings
XMLAtt::XMLAtt()
{
	m_szatt[0] = 0;
	m_szatv[0] = 0;
}

//-----------------------------------------------------------------------------
//! Returns true if the attribute's name is the same as the string that is passed.
//! Note that the comparison is case-sensitive 
bool XMLAtt::operator == (const char* sz)
{
	return (strcmp(m_szatv, sz) == 0);
}

//=============================================================================
// XMLTag
//=============================================================================

//-----------------------------------------------------------------------------
XMLTag::XMLTag()
{
	m_preader = 0;
	m_bend = false;

	m_sztag[0] = 0;
	m_nlevel = 0;

	m_natt = 0;
	for (int i=0; i<MAX_LEVEL; ++i) m_szroot[i][0] = 0;
}

//-----------------------------------------------------------------------------
//! Clear tag data
void XMLTag::clear()
{
	m_sztag[0] = 0;
	m_szval.clear();
	m_natt = 0;
	m_bend = false;
	m_bleaf = true;
	m_bempty = false;
}

//-----------------------------------------------------------------------------
//! This function reads in a comma delimited list of doubles. The function reads
//! in a maximum of n values. The actual number of values that are read is returned.
//!
int XMLTag::value(double* pf, int n)
{
	const char* sz = m_szval.c_str();
	int nr = 0;
	for (int i=0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pf[i] = atof(sz);
		nr++;

		if (sze) sz = sze+1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
//! This function reads in a comma delimited list of floats. The function reads
//! in a maximum of n values. The actual number of values that are read is returned.
//!
int XMLTag::value(float* pf, int n)
{
	const char* sz = m_szval.c_str();
	int nr = 0;
	for (int i=0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pf[i] = (float) atof(sz);
		nr++;

		if (sze) sz = sze+1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
//! This function reads in a comma delimited list of ints. The function reads
//! in a maximum of n values. The actual number of values that are read is returned.
//!
int XMLTag::value(int* pi, int n)
{
	const char* sz = m_szval.c_str();
	int nr = 0;
	for (int i=0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pi[i] = atoi(sz);
		nr++;

		if (sze) sz = sze+1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
void XMLTag::value(bool& val)
{ 
	int n=0; 
	sscanf(m_szval.c_str(), "%d", &n); 
	val = (n != 0); 
}

//-----------------------------------------------------------------------------
void XMLTag::value(char* szstr)
{
	strcpy(szstr, m_szval.c_str()); 
}

//-----------------------------------------------------------------------------
void XMLTag::value(vector<int>& l)
{
	int i, n = 0, n0, n1, nn;
	char* szval = strdup(m_szval.c_str());
	char* ch;
	char* sz = szval;
	int nread;
	do
	{
		ch = strchr(sz, ',');
		if (ch) *ch = 0;
		nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
			break;
		case 3:
			break;
		default:
			n0 = 0;
			n1 = -1;
			nn = 1;
		}

		for (i=n0; i<=n1; i += nn) ++n;

		if (ch) *ch = ',';
		sz = ch+1;
	}
	while (ch != 0);

	if (n != 0)
	{
		l.resize(n);

		sz = szval;
		n = 0;
		do
		{
			ch = strchr(sz, ',');
			if (ch) *ch = 0;
			nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
			switch (nread)
			{
			case 1:
				n1 = n0;
				nn = 1;
				break;
			case 2:
				nn = 1;
			}

			for (i=n0; i<=n1; i += nn) l[n++] = i;
			assert(n <= (int) l.size());

			if (ch) *ch = ',';
			sz = ch+1;
		}
		while (ch != 0);
	}

	free(szval);
}

//-----------------------------------------------------------------------------
//! Return the number of children of a tag
int XMLTag::children()
{
	XMLTag tag(*this); ++tag;
	int ncount = 0;
	while (!tag.isend()) { ncount++; m_preader->SkipTag(tag); }
	return ncount;
}

//-----------------------------------------------------------------------------
//! return the string value of an attribute
const char* XMLTag::AttributeValue(const char* szat, bool bopt)
{
	// find the attribute
	for (int i=0; i<m_natt; ++i)
		if (strcmp(m_att[i].m_szatt, szat) == 0) return m_att[i].m_szatv;

	// If the attribute was not optional, we throw a fit
	if (!bopt) throw XMLReader::MissingAttribute(*this, szat);

	// we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
//! return the attribute
XMLAtt* XMLTag::Attribute(const char* szat, bool bopt)
{
	// find the attribute
	for (int i=0; i<m_natt; ++i)
		if (strcmp(m_att[i].m_szatt, szat) == 0) return m_att+i;

	// If the attribute was not optional, we throw a fit
	if (!bopt) throw XMLReader::MissingAttribute(*this, szat);

	// we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
//! return the attribute
XMLAtt& XMLTag::Attribute(const char* szat)
{
	// find the attribute
	for (int i=0; i<m_natt; ++i)
		if (strcmp(m_att[i].m_szatt, szat) == 0) return m_att[i];

	// throw a fit
	throw XMLReader::MissingAttribute(*this, szat);
}

//-----------------------------------------------------------------------------
bool XMLTag::AttributeValue(const char* szat, double& d, bool bopt)
{
	const char* szv = AttributeValue(szat, bopt);
	if (szv == 0) return false;

	d = atof(szv);

	return true;
}

//-----------------------------------------------------------------------------
bool XMLTag::AttributeValue(const char* szat, int& n, bool bopt)
{
	const char* szv = AttributeValue(szat, bopt);
	if (szv == 0) return false;

	n = atoi(szv);

	return true;
}

//=============================================================================
// XMLReader - Exceptions
//=============================================================================

//-----------------------------------------------------------------------------
void XMLReader::Error::SetErrorString(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	va_start(args, sz);
	vsprintf(m_szerr, sz, args);
	va_end(args);
}

//-----------------------------------------------------------------------------
XMLReader::XMLSyntaxError::XMLSyntaxError()
{
	SetErrorString("syntax error");
}

//-----------------------------------------------------------------------------
XMLReader::UnmatchedEndTag::UnmatchedEndTag(XMLTag& tag)
{
	const char* sz = tag.m_szroot[tag.m_nlevel];
	SetErrorString("unmatched end tag for \"%s\"", sz);
}

//-----------------------------------------------------------------------------
XMLReader::InvalidTag::InvalidTag(XMLTag& tag)
{
	SetErrorString("unrecognized tag \"%s\"", tag.m_sztag);
}

//-----------------------------------------------------------------------------
XMLReader::InvalidValue::InvalidValue(XMLTag& tag)
{
	SetErrorString("the value for tag \"%s\" is invalid", tag.m_sztag);
}

//-----------------------------------------------------------------------------
XMLReader::InvalidAttributeValue::InvalidAttributeValue(XMLTag& tag, const char* sza, const char* szv)
{
	const char* szt = tag.m_sztag;
	if (szv) SetErrorString("invalid value \"%s\" for attribute \"%s.%s\"", szv, szt, sza);
	else SetErrorString("invalid value for attribute \"%s.%s\"", szt, sza);
}

//-----------------------------------------------------------------------------
XMLReader::MissingAttribute::MissingAttribute(XMLTag& tag, const char* sza)
{
	SetErrorString("Missing attribute \"%s\" of tag \"%s\"", sza, tag.m_sztag);
}

//=============================================================================
// XMLReader
//=============================================================================

//-----------------------------------------------------------------------------
XMLReader::XMLReader()
{
	m_fp = 0;
	m_nline = 0;
}

//-----------------------------------------------------------------------------
XMLReader::~XMLReader()
{
	Close();
}

//-----------------------------------------------------------------------------
void XMLReader::Close()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
		m_nline = 0;
	}
}

//-----------------------------------------------------------------------------
bool XMLReader::Open(const char* szfile)
{
	// make sure this reader has not been attached to a file yet
	if (m_fp != 0) return false;

	// open the file
	m_fp = fopen(szfile, "rb");
	if (m_fp == 0) return false;

	// read the first line
	char szline[256] = {0};
	fgets(szline, 255, m_fp);

	// make sure it is correct
	if (strncmp(szline,"<?xml", 5) != 0)
	{
		// This file is not an XML file
		return false;
	}

	// This file is ready to be processed
	return true;
}

//-----------------------------------------------------------------------------

class XMLPath
{
	struct TAG
	{
		char*	tag;
		char*	att;
		char*	atv;
	};

public:
	XMLPath(const char* xpath)
	{
		int l = strlen(xpath);
		m_path = new char[l+1];
		strncpy(m_path, xpath, l);
		m_path[l] = 0;
		m_index = 0;
		parse();
	}

	~XMLPath()
	{
		delete [] m_path;
	}

	void next() { m_index++; }

	bool match(XMLTag& tag)
	{
		// make sure index is valid
		if ((m_index < 0) || (m_index >= (int) m_tag.size())) return false;

		// get current tag
		TAG& t = m_tag[m_index];

		// do tag name compare?
		if (strcmp(t.tag, tag.m_sztag) != 0) return false;

		// if tag requires attribute, make sure it exist
		if (t.att)
		{
			const char* szatt = tag.AttributeValue(t.att, true);
			if (szatt)
			{
				// if the value is specified, make sure it matches too
				if (t.atv)
				{
					if (strcmp(t.atv, szatt) == 0) return true;
					else return false;
				}
				else return false;
			}
			else return false;
		}
		else return true;
	}

	bool valid() { return ((m_index >= 0) && (m_index < m_tag.size())); }

private:
	void parse()
	{
		char* sz = m_path;
		char* ch = strchr(m_path, '/');
		do
		{
			TAG tag;
			tag.tag = sz;
			tag.att = 0;
			tag.atv = 0;

			if (ch) *ch = 0;
			if (sz)
			{
				char* cl = strchr(sz, '[');
				if (cl)
				{
					char* cr = strchr(cl+1, ']');
					if (cr) 
					{
						*cl++ = 0;
						if (cl[0]=='@') cl++;
						tag.att = cl;

						*cr = 0;
						cr = strchr(cl, '=');
						if (cr)
						{
							*cr++ = 0;
							tag.atv = cr;
						}
					}
				}
			}

			sz = (ch ? ch + 1 : 0);
			if (ch) ch = strchr(sz, '/');

			m_tag.push_back(tag);
		}
		while (sz);
	}

private:
	char*	m_path;
	vector<TAG>	m_tag;
	int			m_index;
};

bool XMLReader::FindTag(const char* xpath, XMLTag& tag)
{
	// go to the beginning of the file
	fseek(m_fp, 0, SEEK_SET);

	// set the first tag
	tag.m_preader = this;
	tag.m_ncurrent_line = 1;
	fgetpos(m_fp, &m_currentPos);
	tag.m_fpos = m_currentPos;

	XMLPath path(xpath);

	// find the correct tag
	bool bfound = false;
	try
	{
		// get the next tag
		NextTag(tag);
		do
		{
			// check for match
			if (path.match(tag))
			{
				path.next();
				if (path.valid()) NextTag(tag);
				bfound = true;
			}
			else
			{
				SkipTag(tag);
				bfound = false;
			}
		}
		while (path.valid());
	}
	catch (...)
	{
		// an error has occured (or the end-of-file was reached)
		return false;
	}

	return bfound;
}

//-----------------------------------------------------------------------------
void XMLReader::NextTag(XMLTag& tag)
{
	assert( tag.m_preader == this);

	// set current line number
	m_nline = tag.m_ncurrent_line;

	// set the current file position
//	if (m_currentPos != tag.m_fpos)
	{
		fsetpos(m_fp, &tag.m_fpos);
		m_currentPos = tag.m_fpos;
	}

	// clear tag's content
	tag.clear();

	// read the start tag
	ReadTag(tag);

	try
	{
		// read value and end tag if tag is not empty
		if (!tag.isempty())
		{
			// read the value
			ReadValue(tag);

			// read the end tag
			ReadEndTag(tag);
		}
	}
	catch (EndOfFile)
	{
		if (!tag.isend()) throw UnexpectedEOF();
	}

	// store current line number
	tag.m_ncurrent_line = m_nline;

	// store start file pos for next element
	fgetpos(m_fp, &m_currentPos);
	tag.m_fpos = m_currentPos;
}

//-----------------------------------------------------------------------------
inline bool isvalid(char c)
{
	return (isalnum(c) || (c=='_') || (c=='.'));
}

//-----------------------------------------------------------------------------
void XMLReader::ReadTag(XMLTag& tag)
{
	// find the start token
	char ch, *sz;
	while (true)
	{
		while ((ch=GetChar())!='<') if (!isspace(ch)) throw XMLSyntaxError();

		ch = GetChar();
		if (ch == '!')
		{
			// parse the comment
			ch = GetChar(); if (ch != '-') throw XMLSyntaxError();
			ch = GetChar(); if (ch != '-') throw XMLSyntaxError();

			// find the end of the comment
			int n = 0;
			do
			{
				ch = GetChar();
				if (ch == '-') n++;
				else if ((ch == '>') && (n >= 2)) break;
				else n = 0;
			}
			while (1);
		}
		else if (ch == '?')
		{
			// parse the xml header tag
			while ((ch = GetChar()) != '?');
			ch = GetChar();
			if (ch != '>') throw XMLSyntaxError();
		}
		else break;
	}

	// record the startline
	tag.m_nstart_line = m_nline;

	if (ch == '/')
	{
		tag.m_bend = true;
		ch = GetChar();
	}

	// skip whitespace
	while (isspace(ch)) ch = GetChar();

	// read the tag name
	if (!isvalid(ch)) throw XMLSyntaxError();
	sz = tag.m_sztag;
	*sz++ = ch;
	while (isvalid(ch=GetChar())) *sz++ = ch;
	*sz = 0;

	// read attributes
	tag.m_natt = 0;
	int n = 0;
	sz = 0;
	while (true)
	{
		// skip whitespace
		while (isspace(ch)) ch = GetChar();
		if (ch == '/')
		{
			tag.m_bempty = true;
			ch = GetChar();
			if (ch != '>') throw XMLSyntaxError();
			break;
		}
		else if (ch == '>') break;

		// read the attribute's name
		sz = tag.m_att[n].m_szatt;
		if (!isvalid(ch)) throw XMLSyntaxError();
		*sz++ = ch;
		while (isvalid(ch=GetChar())) *sz++ = ch;
		*sz=0; sz=0;

		// skip whitespace
		while (isspace(ch)) ch=GetChar();
		if (ch != '=') throw XMLSyntaxError();

		// skip whitespace
		while (isspace(ch=GetChar()));

		// make sure the attribute's value starts with an apostrophe or quotation mark
		if ((ch != '"')&&(ch != '\'')) throw XMLSyntaxError();
		char quot = ch;

		sz = tag.m_att[n].m_szatv;
		while ((ch=GetChar())!=quot) *sz++ = ch;
		*sz=0; sz=0;
		ch=GetChar();

		++n;
		++tag.m_natt;
	}

	if (!tag.isend() && !tag.isempty())
	{
		// keep a copy of the name
		strcpy(tag.m_szroot[tag.m_nlevel], tag.m_sztag);
		tag.m_nlevel++;
	}
}

//-----------------------------------------------------------------------------
void XMLReader::ReadValue(XMLTag& tag)
{
	char ch;
	if (!tag.isend())
	{
		tag.m_szval.clear();
		tag.m_szval.reserve(256);
		while ((ch=GetChar())!='<') 
		{ 
			tag.m_szval.push_back(ch);
		}
		tag.m_szval.push_back(0);
	}
	else while ((ch=GetChar())!='<');
}

//-----------------------------------------------------------------------------
void XMLReader::ReadEndTag(XMLTag& tag)
{
	char ch, *sz = tag.m_sztag;
	if (!tag.isend())
	{
		ch = GetChar();
		if (ch == '/')
		{
			// this is the end tag
			// make sure it matches the tag name
			--tag.m_nlevel;

			// skip whitespace
			while (isspace(ch=GetChar()));

			int n = 0;
			do
			{
				if (ch != *sz++) throw UnmatchedEndTag(tag);
				ch = GetChar();
				++n;
			}
			while (!isspace(ch) && (ch!='>'));
			if (n != (int) strlen(tag.m_sztag)) throw UnmatchedEndTag(tag);

			// skip whitespace
			while (isspace(ch)) ch=GetChar();
			if (ch != '>') throw XMLSyntaxError();

			// find the start of the next tag
			if (tag.m_nlevel)
			{
				while (isspace(ch=GetChar()));
				if (ch != '<') throw XMLSyntaxError();
				fseek(m_fp, -1, SEEK_CUR);
			}
		}
		else
		{
			// this element has child elements
			// and therefor is not a leaf

			tag.m_bleaf = false;
			fseek(m_fp, -2, SEEK_CUR);
		}
	}
	else
	{
		fseek(m_fp, -1, SEEK_CUR);

		--tag.m_nlevel;

		// make sure the name is the same as the root
		if (strcmp(tag.m_sztag, tag.m_szroot[tag.m_nlevel]) != 0) throw UnmatchedEndTag(tag);
	}
}

//-----------------------------------------------------------------------------
//! Read the next character in the file.
char XMLReader::GetChar()
{
	char ch;
	while ((ch=fgetc(m_fp))=='\n') ++m_nline;

	// read entity references
	if (ch=='&')
	{
		char szbuf[16]={0};
		szbuf[0] = '&';
		int i = 1;
		do
		{
			ch = fgetc(m_fp);
			szbuf[i++]=ch;
		}
		while ((i<16)&&(ch!=';'));
		if (ch!=';') throw XMLSyntaxError();

		if      (strcmp(szbuf, "&lt;"  )==0) ch = '<';
		else if (strcmp(szbuf, "&gt;"  )==0) ch = '>';
		else if (strcmp(szbuf, "&amp;" )==0) ch = '&';
		else if (strcmp(szbuf, "&apos;")==0) ch = '\'';
		else if (strcmp(szbuf, "&quot;")==0) ch = '"';
		else throw XMLSyntaxError();
	}

	if (feof(m_fp)) 
	{
		int a = 0;
		throw EndOfFile();
	}
	return ch;
}

//-----------------------------------------------------------------------------
//! Skip a tag
void XMLReader::SkipTag(XMLTag& tag)
{
	// if this tag is a leaf we just return
	if (tag.isleaf()) { ++tag; return; }

	// if it is not a leaf we have to loop over all 
	// the children, skipping each child in turn
	NextTag(tag);
	do
	{
		SkipTag(tag);
	}
	while (!tag.isend());

	++tag;
}
