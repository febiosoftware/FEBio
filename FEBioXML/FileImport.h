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
#include <XML/XMLReader.h>
#include <FECore/vec3d.h>
#include <FECore/mat3d.h>
#include <FECore/tens3d.h>
#include <FECore/FECoreBase.h>
#include "FEModelBuilder.h"
#include "febioxml_api.h"
#include <map>

//-----------------------------------------------------------------------------
// Forward declarations
class FEModel;
class FEFileImport;

//-----------------------------------------------------------------------------
// Base class for FEBio import exceptions
// Derived classes should set the error string in their constructor
class FEBIOXML_API FEFileException
{
public:
	enum { MAX_ERR_STRING = 1024 };

	FEFileException();

	FEFileException(const char* sz, ...);

public:
	// retrieve the error string
	const char* GetErrorString() const { return m_szerr; }

protected:
	// set the error string (used by derived classes)
	void SetErrorString(const char* sz, ...);

protected:
	char	m_szerr[MAX_ERR_STRING];
};

//-----------------------------------------------------------------------------
// Class for handling unrecognized tags. 
// Use FEFileSection::SetInvalidTagHandler to set the handler for unrecognized parameters. 
class FEInvalidTagHandler
{
public:
	FEInvalidTagHandler() {}
	virtual ~FEInvalidTagHandler() {}

	virtual bool ProcessTag(XMLTag& tag) { return false; }
};

//-----------------------------------------------------------------------------
// The FEObsoleteParamHandler class tries to map an unrecognized tag to 
// a model parameter. Obsolete parameters are added with the AddParam function. 
class FEObsoleteParamHandler : public FEInvalidTagHandler
{
	struct FEObsoleteParam
	{
		const char* oldName = nullptr;
		const char* newName = nullptr;
		int paramType = FE_PARAM_INVALID;
		bool readIn = false;
		union {
			bool		bVal;
			int			iVal;
			double		gVal;
		};
	};

public:
	FEObsoleteParamHandler(XMLTag& tag, FECoreBase* pc);

	// Add an obsolete parameter.
	// The oldname is a relative path w.r.t. the XMLTag passed in the constructor.
	// the newName is a ParamString, relative to the pc parameter. 
	// To mark a parameter as ignored, set newName to nullptr, and paramType to FE_PARAM_INVALID.
	void AddParam(const char* oldName, const char* newName, FEParamType paramType);

	// This will try to find an entry in the m_param list. If a match is not find, this function returns false.
	bool ProcessTag(XMLTag& tag) override;

	// This function will try to map the obsolete parameters to model parameters. 
	// Or it will print a warning if the obsolete parameter will be ignored. 
	virtual void MapParameters();

private:
	std::string	m_root;
	FECoreBase* m_pc;
	vector<FEObsoleteParam>	m_param;
};

//-----------------------------------------------------------------------------
// Base class for XML sections parsers
class FEBIOXML_API FEFileSection
{
public:
	FEFileSection(FEFileImport* pim) { m_pim = pim; }
    virtual ~FEFileSection() {}

	virtual void Parse(XMLTag& tag) = 0;

	FEFileImport* GetFileReader() { return m_pim; }

	FEModel* GetFEModel();

	FEModelBuilder* GetBuilder();

	// Set the handler for unrecognized tags
	void SetInvalidTagHandler(FEInvalidTagHandler* ith);

public:
	//! read a nodal ID
	//! This assumes the node ID is defined via the "id" attribute
	int ReadNodeID(XMLTag& tag);

public:
	bool ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam = 0, FECoreBase* pc = 0, bool parseAttributes = true);
	bool ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam = 0, bool parseAttributes = true);
	void ReadParameterList(XMLTag& tag, FEParameterList& pl);
	void ReadParameterList(XMLTag& tag, FECoreBase* pc);
	void ReadAttributes(XMLTag& tag, FECoreBase* pc);

public:
	void value(XMLTag& tag, int&    n);
	void value(XMLTag& tag, double& g);
	void value(XMLTag& tag, bool&   b);
	void value(XMLTag& tag, vec3d&  v);
	void value(XMLTag& tag, mat3d&  m);
	void value(XMLTag& tag, mat3ds& m);
	void value(XMLTag& tag, tens3drs& m);
	void value(XMLTag& tag, char* szstr);
	int value(XMLTag& tag, int* pi, int n);
	int value(XMLTag& tag, double* pf, int n);
	void value(XMLTag& tag, std::string& v);
	void value(XMLTag& tag, std::vector<int>& v);
	void value(XMLTag& tag, std::vector<double>& v);

protected:
	bool parseEnumParam(FEParam* pp, const char* val);

private:
	FEFileImport*	m_pim;
	FEInvalidTagHandler* m_ith = nullptr;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FEBIOXML_API FEFileSectionMap : public std::map<string, FEFileSection*>
{
public:
	~FEFileSectionMap();
	void Clear();

	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
//! Base class for file import classes. 
//! FEBio import files are XML formatted files, where each major section (children of root) is represented
//! by an FEFileSection. 
//! This class also offers a simple error reporting mechanism and manages the FILE* pointer. 
//! This class also manages "xml parameters". This is a feature of FEBio files that allow users to use parameters
//! as values for xml tag. A parameter is defined by a name-value pair and referenced in the input file using the $(parameter_name) syntax.

class FEBIOXML_API FEFileImport
{
public:
	//! constructor
	FEFileImport();

	//! destructor
	virtual ~FEFileImport();

	//! get the error message
	void GetErrorMessage(char* szerr);

	//! Get the current FE model that is being processed
	FEModel* GetFEModel();

	//! Get the model builder
	FEModelBuilder* GetBuilder();

	// return the file path
	const char* GetFilePath();

	// set file version
	void SetFileVerion(int nversion);

	// get file version
	int GetFileVersion() const;

	// throw exception if an unknown attribute is found
	void SetStopOnUnknownAttribute(bool b);
	bool StopOnUnknownAttribute() const;

protected:
	//! open a file
	bool Open(const char* szfile, const char* szmode);

	//! close the file
	void Close();

	//! helper function for reporting errors
	bool errf(const char* szerr, ...);

	//! parse the file
	bool ParseFile(XMLTag& tag);

protected:
	FILE*	m_fp;			//!< file pointer
	char	m_szfile[256];	//!< file name
	char	m_szerr[256];	//!< error message
	char	m_szpath[512];	//!< file path

protected:
	FEFileSectionMap	m_map;
	FEModelBuilder*		m_builder;
	bool				m_stopOnUnknownAttribute;

private:
	int	m_nversion;	// version of file
};
