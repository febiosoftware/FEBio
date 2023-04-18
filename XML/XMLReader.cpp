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
#include "XMLReader.h"
#include <assert.h>
#include <stdarg.h>
#include <fstream>
#include <sstream>
using namespace std;

//=============================================================================
// XMLAtt
//=============================================================================

//-----------------------------------------------------------------------------
//! The constructor sets the name and value to zero strings
XMLAtt::XMLAtt()
{
	m_bvisited = false;
}

//-----------------------------------------------------------------------------
//! Returns true if the attribute's value is the same as the string that is passed.
//! Note that the comparison is case-sensitive 
bool XMLAtt::operator == (const char* sz)
{
	return (strcmp(m_val.c_str(), sz) == 0);
}

bool XMLAtt::operator != (const char* szval) 
{ 
	return (strcmp(szval, m_val.c_str()) != 0); 
}

int XMLAtt::value(double* pf, int n)
{
	const char* sz = m_val.c_str();
	int nr = 0;
	for (int i = 0; i < n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pf[i] = atof(sz);
		nr++;

		if (sze) sz = sze + 1;
		else break;
	}
	return nr;
}

//=============================================================================
// XMLTag
//=============================================================================

//-----------------------------------------------------------------------------
XMLTag::XMLTag()
{
	m_preader = 0;
	m_bend = false;
	m_bleaf = true;
	m_bempty = false;
}

//-----------------------------------------------------------------------------
//! Clear tag data
void XMLTag::clear()
{
	m_sztag.clear();
	m_szval.clear();
	m_att.clear();
//	m_path.clear(); // NOTE: Do not clear the path!
	m_bend = false;
	m_bleaf = true;
	m_bempty = false;
}

//-----------------------------------------------------------------------------
std::string XMLTag::relpath(const char* szroot) const
{
	std::string s;
	int m = -1;
	if (szroot)
	{
		for (int i = 0; i < m_path.size(); ++i)
		{
			if (strcmp(szroot, m_path[i].c_str()) == 0)
			{
				m = i;
				break;
			}
		}
	}

	if (m != -1)
	{
		for (int i = m+1; i < m_path.size(); i++)
		{
			if (i != m+1) s += ".";
			s += m_path[i];
		}
	}

	if (s.empty() == false) s += ".";
	s += m_sztag;

	return s;
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
int XMLTag::value(std::vector<string>& stringList, int n)
{
	stringList.clear();

	char tmp[256] = { 0 };

	const char* sz = m_szval.c_str();
	int nr = 0;
	for (int i = 0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');
		if (sze)
		{
			int l = (int)(sze - sz);
			if (l > 0) strncpy(tmp, sz, l);
			tmp[l] = 0;

			stringList.push_back(tmp);
		}
		else stringList.push_back(sz);
		nr++;

		if (sze) sz = sze + 1;
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
void XMLTag::value(std::string& val)
{
	val = m_szval;
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
				break;
			case 3:
				break;
			default:
				n0 = 0;
				n1 = -1;
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
void XMLTag::value(vector<double>& l)
{
	l.clear();
	const char *sz = m_szval.c_str();
	while (sz && *sz)
	{
		// skip space
		while (*sz == ' ') ++sz;

		// read the value
		if (sz && *sz)
		{
			double v = atof(sz);
			l.push_back(v);

			// find next space or comma
			while (*sz && (*sz != ' ') && (*sz != ',')) sz++;
			if (*sz == ',') sz++;
		}
	}
}

void XMLTag::value2(std::vector<int>& l)
{
	l.clear();
	const char* sz = m_szval.c_str();
	while (sz && *sz)
	{
		// skip space
		while (*sz == ' ') ++sz;

		// read the value
		if (sz && *sz)
		{
			int v = (int)atoi(sz);
			l.push_back(v);

			// find next space or comma
			while (*sz && (*sz != ' ') && (*sz != ',')) sz++;
			if (*sz == ',') sz++;
		}
	}
}

//-----------------------------------------------------------------------------
//! Return the number of children of a tag
int XMLTag::children()
{
	XMLTag tag(*this); ++tag;
	int ncount = 0;
	while (!tag.isend()) { ncount++; tag.skip(); ++tag; }
	return ncount;
}

//-----------------------------------------------------------------------------
//! return the string value of an attribute
const char* XMLTag::AttributeValue(const char* szat, bool bopt)
{
	// find the attribute
	for (XMLAtt& att : m_att)
	{
		if (strcmp(att.m_name.c_str(), szat) == 0)
		{
			att.m_bvisited = true;
			return att.m_val.c_str();
		}
	}

	// If the attribute was not optional, we throw a fit
	if (!bopt) throw XMLReader::MissingAttribute(*this, szat);

	// we didn't find it
	return 0;
}

//-----------------------------------------------------------------------------
XMLAtt* XMLTag::AttributePtr(const char* szat)
{
	return Attribute(szat, true);
}

//-----------------------------------------------------------------------------
//! return the attribute
XMLAtt* XMLTag::Attribute(const char* szat, bool bopt)
{
	// find the attribute
	for (XMLAtt& att : m_att)
	{
		if (strcmp(att.m_name.c_str(), szat) == 0)
		{
			att.m_bvisited = true;
			return &(att);
		}
	}

	// If the attribute was not optional, we throw a fit
	if (!bopt) throw XMLReader::MissingAttribute(*this, szat);

	// we didn't find it
	return nullptr;
}

//-----------------------------------------------------------------------------
//! return the attribute
XMLAtt& XMLTag::Attribute(const char* szat)
{
	// find the attribute
	for (XMLAtt& att : m_att)
	{
		if (strcmp(att.m_name.c_str(), szat) == 0)
		{
			att.m_bvisited = true;
			return att;
		}
	}

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

// helper function for formatting a string
string format_string(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	char szbuf[512];
	va_start(args, sz);
	vsprintf(szbuf, sz, args);
	va_end(args);

	return szbuf;
}

//-----------------------------------------------------------------------------
XMLReader::Error::Error(XMLTag& tag, const std::string& err) : \
std::runtime_error(format_string("tag \"%s\" (line %d) : ", tag.Name(), tag.m_nstart_line) + err) {}

//-----------------------------------------------------------------------------
XMLReader::XMLSyntaxError::XMLSyntaxError(int line_number) : \
XMLReader::Error(format_string("syntax error (line %d)", line_number)) {}

//-----------------------------------------------------------------------------
XMLReader::UnmatchedEndTag::UnmatchedEndTag(XMLTag& tag) : \
XMLReader::Error(tag, "unmatched end tag") {}

//-----------------------------------------------------------------------------
XMLReader::InvalidTag::InvalidTag(XMLTag& tag) : \
XMLReader::Error(tag, "unrecognized tag") {}

//-----------------------------------------------------------------------------
XMLReader::InvalidValue::InvalidValue(XMLTag& tag) : \
XMLReader::Error(tag, format_string("invalid value: %s", tag.isleaf() ? tag.szvalue() : "")) {}

//-----------------------------------------------------------------------------
XMLReader::InvalidAttributeValue::InvalidAttributeValue(XMLTag& tag, const char* sza, const char* szv) : \
XMLReader::Error(tag, format_string("invalid value for attribute \"%s\"", sza)) {}

XMLReader::InvalidAttributeValue::InvalidAttributeValue(XMLTag& tag, XMLAtt& att) : \
XMLReader::Error(tag, format_string("invalid value for attribute \"%s\"", att.cvalue())) {}

//-----------------------------------------------------------------------------
XMLReader::InvalidAttribute::InvalidAttribute(XMLTag& tag, const char* sza) :\
XMLReader::Error(tag, format_string("invalid attribute \"%s\"", sza)) {}

//-----------------------------------------------------------------------------
XMLReader::MissingAttribute::MissingAttribute(XMLTag& tag, const char* sza) : \
XMLReader::Error(tag, format_string("missing attribute \"%s\"", sza)) {}

//-----------------------------------------------------------------------------
XMLReader::MissingTag::MissingTag(XMLTag& tag, const char* sza) : \
XMLReader::Error(tag, format_string("missing tag \"%s\"", sza)) {}

//=============================================================================
// XMLReader
//=============================================================================

//-----------------------------------------------------------------------------
XMLReader::XMLReader()
{
    m_stream = nullptr;
	m_nline = 0;
	m_bufIndex = 0;
	m_bufSize = 0;
	m_eof = false;
	m_currentPos = 0;
	m_buf = new char[BUF_SIZE];
}

//-----------------------------------------------------------------------------
XMLReader::~XMLReader()
{
	Close();
	delete[] m_buf;
}

//-----------------------------------------------------------------------------
void XMLReader::Close()
{
    if(m_stream)
    {
        ifstream* fileStream = dynamic_cast<ifstream*>(m_stream);
        if(fileStream)
        {
            fileStream->close();
        }

        delete m_stream;
        m_stream = nullptr;
    }

	m_nline = 0;
	m_bufIndex = 0;
	m_bufSize = 0;
	m_eof = false;
	m_currentPos = 0;
}

//-----------------------------------------------------------------------------
bool XMLReader::Open(const char* szfile, bool checkForXMLTag)
{
	// make sure this reader has not been attached to a file yet
    if(m_stream) return false;

	// open the file
    m_stream = new ifstream;
    static_cast<ifstream*>(m_stream)->open(szfile, ifstream::in|ifstream::binary);
    if(m_stream->fail()) return false;


	// read the first line
	if (checkForXMLTag)
	{
		char szline[256] = { 0 };
		// fgets(szline, 255, m_fp);
        m_stream->get(szline, 255);

		// make sure it is correct
		if (strncmp(szline, "<?xml", 5) != 0)
		{
			// This file is not an XML file
			return false;
		}
	}

	m_currentPos = 0;

	// This file is ready to be processed
	return true;
}

bool XMLReader::OpenString(std::string& xml, bool checkForXMLTag)
{
    // make sure this we don't already have a stream
    if(m_stream) return false;

    // create string stream
    m_stream = new std::istringstream(xml);
    if(m_stream->fail()) return false;

    // read the first line
	if (checkForXMLTag)
	{
		char szline[256] = { 0 };
		// fgets(szline, 255, m_fp);
        m_stream->get(szline, 255);

		// make sure it is correct
		if (strncmp(szline, "<?xml", 5) != 0)
		{
			// This file is not an XML file
			return false;
		}
	}

	m_currentPos = 0;

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
		int l = (int)strlen(xpath);
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
		if (strcmp(t.tag, tag.Name()) != 0) return false;

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
    m_stream->seekg(0, ios_base::beg);
	m_bufIndex = m_bufSize = 0;
	m_currentPos = 0;
	m_eof = false;

	// set the first tag
	tag.m_preader = this;
	tag.m_ncurrent_line = 1;
	tag.m_fpos = currentPos();

	XMLPath path(xpath);

	// find the correct tag
	bool bfound = false;
	try
	{
		// get the next tag
		++tag;
		do
		{
			// check for match
			if (path.match(tag))
			{
				path.next();
				if (path.valid()) ++tag;
				bfound = true;
			}
			else
			{
				tag.skip();
				++tag;
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
	if (m_currentPos != tag.m_fpos)
	{
        m_stream->seekg(tag.m_fpos, ios_base::beg);
		m_currentPos = tag.m_fpos;
		m_bufSize = m_bufIndex = 0;
		m_eof = false;
	}

	// update the path
	if (!tag.isend() && !tag.isempty() && !tag.isleaf())
	{
		// keep a copy of the name
		tag.m_path.push_back(tag.m_sztag);
	}
	else if (tag.isend())
	{
		// make sure the name is the same as the root
		if (strcmp(tag.m_sztag.c_str(), tag.m_path.back().c_str()) != 0) throw UnmatchedEndTag(tag);

		tag.m_path.erase(--tag.m_path.end());
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
	tag.m_fpos = currentPos();
}

//-----------------------------------------------------------------------------
inline bool isvalid(char c)
{
	return (isalnum(c) || (c=='_') || (c=='.') || (c=='-') || (c==':'));
}

//-----------------------------------------------------------------------------
void XMLReader::ReadTag(XMLTag& tag)
{
	// find the start token
	char ch;
	while (true)
	{
		while ((ch=GetChar())!='<') 
			if (!isspace(ch)) 
			{
				throw XMLSyntaxError(m_nline);
			}

		ch = GetChar();
		if (ch == '!')
		{
            m_comment.clear();

			// parse the comment
			ch = GetChar(); if (ch != '-') throw XMLSyntaxError(m_nline);
			ch = GetChar(); if (ch != '-') throw XMLSyntaxError(m_nline);

			// find the end of the comment
			int n = 0;
			do
			{
				ch = GetChar();
				if (ch == '-') n++;
				else if ((ch == '>') && (n >= 2)) break;
				else
				{
					if (n > 0) m_comment += '-';
					if (n > 1) m_comment += '-';
					if (ch != '\r') m_comment += ch; // don't append \r
					n = 0;
				}
			}
			while (1);

			// eat whitespace at the start and end of the comment
			if (m_comment.empty() == false)
			{
				int n = m_comment.size();
				char* tmp = new char[n + 1];
				strncpy(tmp, m_comment.c_str(), n);
				tmp[n] = 0;

				char* cl = tmp;
				while (*cl && isspace(*cl)) cl++;
				n = strlen(cl);
				if (n > 0)
				{
					char* cr = &cl[n - 1];
					while ((cr > cl) && isspace(*cr)) *cr-- = 0;
				}

				m_comment = cl;
				delete[] tmp;
			}
		}
		else if (ch == '?')
		{
			// parse the xml header tag
			while ((ch = GetChar()) != '?');
			ch = GetChar();
			if (ch != '>') throw XMLSyntaxError(m_nline);
		}
		else break;
	}

	// record the startline
	tag.m_nstart_line = m_nline;

	if (ch == '/')
	{
		tag.m_bend = true;
        m_comment.clear();
		ch = GetChar();
	}

	// skip whitespace
	while (isspace(ch)) ch = GetChar();

	// read the tag name
	if (!isvalid(ch)) throw XMLSyntaxError(m_nline);
	tag.m_sztag.clear();
	tag.m_sztag.push_back(ch);
	while (isvalid(ch=GetChar())) tag.m_sztag.push_back(ch);

	// read attributes
	tag.m_att.clear();
	int n = 0;
	while (true)
	{
		// skip whitespace
		while (isspace(ch)) ch = GetChar();
		if (ch == '/')
		{
			tag.m_bempty = true;
			ch = GetChar();
			if (ch != '>') throw XMLSyntaxError(m_nline);
			break;
		}
		else if (ch == '>') break;

		// read the attribute's name
		XMLAtt att;
		std::string& name = att.m_name;
		name.clear();
		if (!isvalid(ch)) throw XMLSyntaxError(m_nline);
		name.push_back(ch);
		while (isvalid(ch=GetChar())) name.push_back(ch);

		// skip whitespace
		while (isspace(ch)) ch=GetChar();
		if (ch != '=') throw XMLSyntaxError(m_nline);

		// skip whitespace
		while (isspace(ch=GetChar()));

		// make sure the attribute's value starts with an apostrophe or quotation mark
		if ((ch != '"')&&(ch != '\'')) throw XMLSyntaxError(m_nline);
		char quot = ch;

		// read the value
		std::string& val = att.m_val;
		while ((ch=GetChar())!=quot) val.push_back(ch);
		ch=GetChar();

		// mark tag as unvisited
		att.m_bvisited = false;

		++n;
		tag.m_att.push_back(att);
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
	}
	else while ((ch=GetChar())!='<');
}

//-----------------------------------------------------------------------------
void XMLReader::ReadEndTag(XMLTag& tag)
{
	char ch;
	const char* sz = tag.m_sztag.c_str();
	if (!tag.isend())
	{
		ch = GetChar();
		if (ch == '/')
		{
			// this is the end tag

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
			if (n != (int) tag.m_sztag.size()) throw UnmatchedEndTag(tag);

			// skip whitespace
			while (isspace(ch)) ch=GetChar();
			if (ch != '>') throw XMLSyntaxError(m_nline);

			// find the start of the next tag
			if (tag.m_path.empty() == false)
			{
				while (isspace(ch=GetChar()));
				if (ch != '<') throw XMLSyntaxError(m_nline);
				rewind(1);
			}
		}
		else
		{
			// this element has child elements
			// and therefor is not a leaf

			tag.m_bleaf = false;
			rewind(2);
		}
	}
	else
	{
		rewind(1);
	}
}

//-----------------------------------------------------------------------------
char XMLReader::readNextChar()
{
	if (m_bufIndex >= m_bufSize)
	{
		if (m_eof) 
        {
            throw EndOfFile();
        }

        m_stream->read(m_buf, BUF_SIZE);
        
        // clear eofbit and failbit if eof is hit. Unless we do this, seekg doesn't work
        // we handle eof manually, so we can clear it without issue. 
        m_stream->clear();
        
        m_bufSize = m_stream->gcount();
        m_bufIndex = 0;
        m_eof = (m_bufSize != BUF_SIZE);
	}
	m_currentPos++;
	return m_buf[m_bufIndex++];
}

//-----------------------------------------------------------------------------
int64_t XMLReader::currentPos()
{
	return m_currentPos;
}

//-----------------------------------------------------------------------------
//! move the file pointer
void XMLReader::rewind(int64_t nstep)
{
	m_bufIndex -= nstep;
	m_currentPos -= nstep;

	if (m_bufIndex < 0)
	{
        m_stream->seekg(m_bufIndex - m_bufSize, ios_base::cur);
		m_bufIndex = m_bufSize = 0;
		m_eof = false;
	}
}

//-----------------------------------------------------------------------------
//! Read the next character in the file.
char XMLReader::GetChar()
{
	char ch;
    if ((ch=readNextChar())=='\n') ++m_nline;

	// read entity references
	if (ch=='&')
	{
		char szbuf[16]={0};
		szbuf[0] = '&';
		int i = 1;
		do
		{
			ch = readNextChar();
			szbuf[i++]=ch;
		}
		while ((i<16)&&(ch!=';'));
		if (ch!=';') throw XMLSyntaxError(m_nline);

		if      (strcmp(szbuf, "&lt;"  )==0) ch = '<';
		else if (strcmp(szbuf, "&gt;"  )==0) ch = '>';
		else if (strcmp(szbuf, "&amp;" )==0) ch = '&';
		else if (strcmp(szbuf, "&apos;")==0) ch = '\'';
		else if (strcmp(szbuf, "&quot;")==0) ch = '"';
		else throw XMLSyntaxError(m_nline);
	}
	return ch;
}

//-----------------------------------------------------------------------------
//! Skip a tag
void XMLReader::SkipTag(XMLTag& tag)
{
	// if this tag is a leaf we just return
	if (tag.isleaf()) { return; }

	// if it is not a leaf we have to loop over all 
	// the children, skipping each child in turn
	++tag;
	do
	{
		SkipTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
char XMLReader::GetNextChar()
{
	char ch;
	do
	{
		ch = readNextChar();
		if (ch == '\n') ++m_nline;
	} while (ch == '\r');
	return ch;
}

ifstream* XMLReader::GetFileStream()
{
    return dynamic_cast<ifstream*>(m_stream);
}

//! return the current line
int XMLReader::GetCurrentLine() { return m_nline; }


const std::string& XMLReader::GetLastComment()
{
    // Get rid of starting and trailing new lines in comments
	if (m_comment.empty() == false)
	{
		if (*m_comment.begin() == '\n')
		{
			m_comment.erase(m_comment.begin());
		}
		if (*m_comment.rbegin() == '\n')
		{
			m_comment.erase(m_comment.end());
		}
	}

    return m_comment;
}