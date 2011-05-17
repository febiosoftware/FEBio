// XMLReader.cpp: implementation of the XMLReader class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "XMLReader.h"
#include <assert.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// XMLTag
//////////////////////////////////////////////////////////////////////

XMLTag::XMLTag()
{
	m_preader = 0;
	m_bend = false;

	m_sztag[0] = 0;
	m_szval[0] = 0;
	m_nlevel = 0;

	m_natt = 0;
	int i;
	for (i=0; i<XMLReader::MAX_ATT; ++i)
	{
		m_szatt[i][0] = 0;
		m_szatv[i][0] = 0;
	}

	for (i=0; i<XMLReader::MAX_LEVEL; ++i) m_szroot[i][0] = 0;
}

//////////////////////////////////////////////////////////////////////

void XMLTag::value(double* pf, int n)
{
	char* sz = m_szval;

	for (int i=0; i<n; ++i)
	{
		char* sze = strchr(sz, ',');
		if (sze) *sze = 0;

		pf[i] = atof(sz);

		if (sze)
		{
			*sze = ',';
			sz = sze+1;
		}
	}
}

//////////////////////////////////////////////////////////////////////

void XMLTag::value(float* pf, int n)
{
	char* sz = m_szval;

	for (int i=0; i<n; ++i)
	{
		char* sze = strchr(sz, ',');
		if (sze) *sze = 0;

		pf[i] = (float) atof(sz);

		if (sze)
		{
			*sze = ',';
			sz = sze+1;
		}
	}
}

//////////////////////////////////////////////////////////////////////

void XMLTag::value(int* pi, int n)
{
	char* sz = m_szval;

	for (int i=0; i<n; ++i)
	{
		char* sze = strchr(sz, ',');
		if (sze) *sze = 0;

		pi[i] = atoi(sz);

		if (sze)
		{
			*sze = ',';
			sz = sze+1;
		}
	}
}

//////////////////////////////////////////////////////////////////////

void XMLTag::value(vec3d& v)
{
	int n = sscanf(m_szval, "%lg,%lg,%lg", &v.x, &v.y, &v.z);
	if (n != 3) throw XMLReader::XMLSyntaxError();
}

//////////////////////////////////////////////////////////////////////

const char* XMLTag::AttributeValue(const char* szat, bool bopt)
{
	// find the attribute
	for (int i=0; i<m_natt; ++i)
		if (strcmp(m_szatt[i], szat) == 0) return m_szatv[i];

	// If the attribute was not optional, we throw a fit
	if (!bopt) throw XMLReader::MissingAttribute(*this, szat);

	// we didn't find it
	return 0;
}

//////////////////////////////////////////////////////////////////////

bool XMLTag::AttributeValue(const char* szat, double& d, bool bopt)
{
	const char* szv = AttributeValue(szat, bopt);
	if (szv == 0) return false;

	d = atof(szv);

	return true;
}

//////////////////////////////////////////////////////////////////////

bool XMLTag::AttributeValue(const char* szat, int& n, bool bopt)
{
	const char* szv = AttributeValue(szat, bopt);
	if (szv == 0) return false;

	n = atoi(szv);

	return true;
}


//////////////////////////////////////////////////////////////////////
// XMLReader
//////////////////////////////////////////////////////////////////////

XMLReader::XMLReader()
{
	m_fp = 0;
	m_nline = 0;
}

XMLReader::~XMLReader()
{
	Close();
}

//////////////////////////////////////////////////////////////////////

void XMLReader::Close()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
		m_nline = 0;
	}
}

//////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////

bool XMLReader::FindTag(const char* sztag, XMLTag& tag)
{
	// go to the beginning of the file
	fseek(m_fp, 0, SEEK_SET);

	// set the first tag
	tag.m_preader = this;
	tag.m_ncurrent_line = 1;
	fgetpos(m_fp, &tag.m_fpos);

	// find the correct tag
	bool bfound = false;
	do
	{
		NextTag(tag);
		if (strcmp(sztag, tag.m_sztag) == 0) bfound = true;
	}
	while (!bfound);

	return true;
}

//////////////////////////////////////////////////////////////////////

void XMLReader::NextTag(XMLTag& tag)
{
	assert( tag.m_preader == this);

	// set current line number
	m_nline = tag.m_ncurrent_line;

	// set the current file position
	fsetpos(m_fp, &tag.m_fpos);

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
	catch (UnexpectedEOF)
	{
		if (!tag.isend()) throw;
	}

	// store current line number
	tag.m_ncurrent_line = m_nline;

	// store start file pos for next element
	fgetpos(m_fp, &tag.m_fpos);
}

inline bool isvalid(char c)
{
	return (isalnum(c) || (c=='_') || (c=='.'));
}

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
			do
			{
				ch = GetChar();
				if (ch == '-')
				{
					ch = GetChar();
					if (ch == '-')
					{
						ch = GetChar();
						if (ch == '>') break;
					}
				}
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
		sz = tag.m_szatt[n];
		if (!isvalid(ch)) throw XMLSyntaxError();
		*sz++ = ch;
		while (isvalid(ch=GetChar())) *sz++ = ch;
		*sz=0; sz=0;

		// skip whitespace
		while (isspace(ch)) ch=GetChar();
		if (ch != '=') throw XMLSyntaxError();

		// skip whitespace
		while (isspace(ch=GetChar()));
		if (ch != '"') throw XMLSyntaxError();

		sz = tag.m_szatv[n];
		while ((ch=GetChar())!='"') *sz++ = ch;
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

void XMLReader::ReadValue(XMLTag& tag)
{
	char ch;
	int n = 0;
	if (!tag.isend())
	{
		char *sz = tag.m_szval;
		while ((ch=GetChar())!='<') 
		{ 
			*sz++ = ch; n++; 
			if (n > MAX_TAG) throw EndOfBuffer(tag);
		}
		*sz=0;
	}
	else while ((ch=GetChar())!='<');
}

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
