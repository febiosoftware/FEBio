#pragma once
#include "FileImport.h"
#include "XMLReader.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FESolver.h"
#include <map>
#include <string>
using namespace std;

class FEFEBioImport;

//-----------------------------------------------------------------------------
// Base class for XML sections parsers
class FEBioFileSection
{
public:
	FEBioFileSection(FEFEBioImport* pim) { m_pim = pim; }

	virtual void Parse(XMLTag& tag) = 0;

	FEModel* GetFEModel();
	FECore::FEAnalysis* GetStep();

protected:
	FEFEBioImport*	m_pim;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FEBioFileSectionMap : public map<string, FEBioFileSection*>
{
public:
	~FEBioFileSectionMap();
	void Clear();
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEFEBioImport : public FEFileImport
{
public:
	class InvalidVersion{};
	class InvalidMaterial
	{ 
	public: 
		InvalidMaterial(int nel) : m_nel(nel){}
		int m_nel; 
	};
	class InvalidDomainType{};
	class FailedCreatingDomain{};
	class InvalidElementType{};
	class FailedLoadingPlugin
	{
	public:
		FailedLoadingPlugin(const char* sz) : m_szfile(sz) {}
		const char* FileName() { return m_szfile.c_str(); }
	public:
		string	m_szfile;
	};
	class DuplicateMaterialSection {};
	class InvalidDomainMaterial
	{ 
	public: 
		InvalidDomainMaterial(int ndom) : m_ndom(ndom){}
		int m_ndom; 
	};

public:
	//-------------------------------------------------------------------------
	class PlotVariable
	{
	public:
		PlotVariable(const char* szvar, vector<int>& item);
		PlotVariable(const PlotVariable& pv);

	public:
		char		m_szvar[64];	//!< name of output variable
		vector<int>	m_item;			//!< (optional) list of items
	};

	//-------------------------------------------------------------------------
	class XMLParam
	{
		enum {MAX_TAG = 128};
	public:
		XMLParam() { m_szname[0] = 0; m_szval[0] = 0; }
	public:
		char	m_szname[MAX_TAG];	// parameter name
		char	m_szval[MAX_TAG];	// parameter value
	};


public:
	//! constructor
	FEFEBioImport();

	//! Load the model data from file.
	bool Load(FEModel& fem, const char* szfile);

	//! read the contents of a file
	bool ReadFile(const char* szfile);

public:
	FEModel* GetFEModel() { return m_pfem; }
	FEMesh* GetFEMesh() { return m_pMesh; }
	FECore::FEAnalysis*	GetStep();

	int Version() { return m_nversion; }

	bool ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam = 0);
	bool ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam = 0);

	void ReadList(XMLTag& tag, vector<int>& l);

	FECore::FEAnalysis* CreateNewStep();

public:
	void SetDumpfileName(const char* sz) { sprintf(m_szdmp, sz); }
	void SetLogfileName (const char* sz) { sprintf(m_szlog, sz); }
	void SetPlotfileName(const char* sz) { sprintf(m_szplt, sz); }

	void AddPlotVariable(const char* szvar, vector<int>& item);

	void SetPlotCompression(int n);

public:
	void ClearParams();

	XMLParam* FindParameter(const char* sz);

	void AddParameter(const char* szname, const char* szval);

public:
	const char* get_value_string(XMLTag& tag);
	void value(XMLTag& tag, int&    n);
	void value(XMLTag& tag, double& g);
	void value(XMLTag& tag, bool&   b);
	void value(XMLTag& tag, vec3d&  v);
	void value(XMLTag& tag, mat3d&  m);
	void value(XMLTag& tag, mat3ds& m);
	void value(XMLTag& tag, char* szstr);
	int value(XMLTag& tag, int* pi, int n);
	int value(XMLTag& tag, double* pf, int n);

protected:
	void ParseVersion			(XMLTag& tag);

public:
	FEModel*			m_pfem;		//!< pointer to the fem class
	FECore::FEAnalysis*	m_pStep;	//!< pointer to current analysis step
	FEMesh*				m_pMesh;	//!< pointer to the mesh class

public:
	char	m_szpath[512];
	char	m_szdmp[512];
	char	m_szlog[512];
	char	m_szplt[512];

public:
	char					m_szplot_type[256];
	vector<PlotVariable>	m_plot;
	int						m_nplot_compression;

	vector<XMLParam>	m_Param;	// parameter list

public:
	int m_nsteps;		//!< nr of step sections read
	int	m_nstep_type;	//!< step type
	int	m_maxid;		//!< max element ID

	bool	m_b3field;	//!< three-field element flag
	bool	m_but4;		//!< use UT4 formulation flag
	FE_Element_Type		m_nhex8;	//!< hex integration rule
	FE_Element_Type		m_ntet4;	//!< tet4 integration rule
	FE_Element_Type		m_ntet10;	//!< tet10 integration rule
	FE_Element_Type		m_ntet15;	//!< tet15 integration rule
	FE_Element_Type		m_ntri6;	//!< tri6 integration rule
	FE_Element_Type		m_ntri7;	//!< tri7 integration rule
	FE_Element_Type		m_ntri3;	//!< tri3 integration rule

protected:
	FEBioFileSectionMap	m_map;

protected:
	int	m_nversion;	// version of file
};
