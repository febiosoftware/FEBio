#pragma once
#include <stdio.h>
#include "XMLReader.h"
#include <FECore/vec3d.h>
#include <FECore/mat3d.h>
#include <FECore/tens3d.h>
#include <FECore/FECoreBase.h>
#include "FEModelBuilder.h"
#include "febioxml_api.h"
#include <map>
using namespace std;

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

//-------------------------------------------------------------------------
class FEBIOXML_API FEFileParam
{
	enum { MAX_TAG = 128 };

public:
	FEFileParam() { m_szname[0] = 0; m_szval[0] = 0; }

public:
	char	m_szname[MAX_TAG];	// parameter name
	char	m_szval[MAX_TAG];	// parameter value
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
	const char* get_value_string(XMLTag& tag);
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

private:
	FEFileImport*	m_pim;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FEBIOXML_API FEFileSectionMap : public map<string, FEFileSection*>
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

public:
	// Add a file parameter
	void AddFileParameter(const char* szname, const char* szval);

	// find a file parameter
	FEFileParam* FindFileParameter(const char* sz);

	// clear all file parameters
	void ClearFileParams();

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
	FEModel*			m_fem;
	FEModelBuilder*		m_builder;
	vector<FEFileParam>	m_Param;	// parameter list
	bool				m_stopOnUnknownAttribute;

private:
	int	m_nversion;	// version of file
};
